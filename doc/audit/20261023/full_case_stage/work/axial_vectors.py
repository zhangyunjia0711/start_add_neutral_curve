#!/usr/bin/env python3
"""Vector-tracked axial seeds; the original axial script/output remain unchanged.

Execution requires an explicitly supplied immutable solver path and SHA-256.
Only a reviewed v4-or-later build may be used for the first diagnostic blocks.
"""
import argparse
import concurrent.futures
import configparser
import hashlib
import inspect
import json
from pathlib import Path
import re
import sys

import axial_seeds as axial
import full_driver as driver

ROOT = driver.ROOT
OUTPUT = ROOT / "outputs" / "axial_vectors.json"
OVERLAP = 0.95
OPTIONS = ["-pep_tol", "1e-12", "-pep_max_it", "1000",
           "-st_pc_factor_mat_solver_type", "superlu", "-pep_ncv", "8"]
ARGS = None
DRIVER_KEYS = None
BLOCK_CONFIGURATION = {}


def require_driver_keys(adaptive_key, vector_key, overlap_key):
    """Fail before solving when this driver's explicit spec mapping is absent."""
    source = inspect.getsource(driver.run)
    available = set(re.findall(r"spec\.get\(\s*['\"]([^'\"]+)['\"]", source))
    required = (adaptive_key, vector_key, overlap_key, "minimum_step")
    missing = [key for key in required if key not in available]
    if missing:
        raise RuntimeError("Unconfirmed full_driver spec mapping: " + repr(missing)
                           + "; available keys: " + repr(sorted(available)))
    return dict(adaptive=adaptive_key, vector=vector_key, overlap=overlap_key)


def inspect_config(path):
    cfg = configparser.ConfigParser()
    cfg.read(path)
    neutral = dict(cfg["neutral"])
    adaptive = {key: value for key, value in neutral.items() if "adaptive" in key}
    vectors = {key: value for key, value in neutral.items()
               if "vector" in key or "overlap" in key}
    true = lambda value: value.strip().casefold() in ("t", "true", "1", "yes")
    if not any(true(value) for value in adaptive.values()):
        raise RuntimeError("adaptive stepping was not enabled in actual config: " + repr(adaptive))
    if not any(true(value) for key, value in vectors.items() if "track" in key):
        raise RuntimeError("vector tracking was not enabled in actual config: " + repr(vectors))
    thresholds = []
    for key, value in vectors.items():
        if "overlap" in key:
            try:
                thresholds.append(float(value))
            except ValueError:
                pass
    if not thresholds or not any(abs(value - OVERLAP) < 1e-12 for value in thresholds):
        raise RuntimeError("Actual vector-overlap threshold is not 0.95: " + repr(vectors))
    return dict(adaptive=adaptive, vectors=vectors)


def execute_vector_block(family, lo, hi, anchor, direction, omega_center=None):
    tag = ARGS.solver_sha[:8]
    name = f"axial_vectors_{family}_{direction}_i{lo:03d}_{hi:03d}_{tag}"
    center = axial.NOMINAL[family][1] if omega_center is None else omega_center
    width = axial.EPS if omega_center is None else 1e-12
    spec = dict(name=name, solver=str(ARGS.solver), coordinates=ARGS.coordinates, mode="M1", family=family,
                station=lo, seed_station=anchor["iloc"], task=33,
                outer_parameter="iloc", inner_parameter="omega",
                outer_min=lo, outer_max=hi, outer_points=hi - lo + 1,
                inner_min=center - width, inner_max=center + width, inner_points=2,
                step=.2, jump=2, minimum_step=1e-5, seed=axial.seed_dict(anchor),
                options=OPTIONS, timeout=900)
    spec[DRIVER_KEYS["adaptive"]] = True
    spec[DRIVER_KEYS["vector"]] = True
    spec[DRIVER_KEYS["overlap"]] = OVERLAP
    summary = driver.run(spec)
    configuration = inspect_config(ROOT / "runs" / name / "config.ini")
    result, records, endpoint = axial.audit_run(summary, family, lo, hi, anchor, direction)
    result["actual_tracking_configuration"] = configuration
    result["solver_sha256"] = ARGS.solver_sha
    result["vector_overlap_metric"] = "normalized Euclidean solver-vector inner product"
    BLOCK_CONFIGURATION[name] = configuration
    return result, records, endpoint


def save_vector_notes():
    axial.STATE["tracking"] = dict(
        adaptive_steps=True, vector_tracking=True, minimum_euclidean_overlap=OVERLAP,
        driver_spec_keys=DRIVER_KEYS, solver=str(ARGS.solver),
        solver_sha256=ARGS.solver_sha, coordinates=ARGS.coordinates,
        qualification="Per-step Euclidean vector overlap is a continuation check; it is not a physical-energy norm or a full-spectrum certificate.")
    axial.STATE["notes"] = list(dict.fromkeys(axial.STATE.get("notes", []) + [
        "This output is independent of outputs/axial_seeds.json; the original axial files are unchanged.",
        "The immutable solver hash is explicitly checked before execution.",
        "M3 i120 positive/negative tiny-frequency side checks are independent actual solves, not interpolated across a failure.",
    ]))
    with axial.LOCK:
        temporary = OUTPUT.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(axial.STATE, indent=2, allow_nan=False) + "\n")
        temporary.replace(OUTPUT)


def i120_side_checks():
    family = axial.STATE["families"].get("M3", {})
    anchor = family.get("seed")
    if not anchor or anchor["status"] != 0:
        axial.STATE["i120_side_checks"] = dict(status="not run", reason="M3 initial main-grid seed not accepted")
        save_vector_notes()
        return
    checks = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        jobs = {}
        for label, omega in (("negative", -1e-8), ("positive", 1e-8)):
            job = executor.submit(execute_vector_block, "M3", 120, 120, anchor,
                                  "i120_side_" + label, omega)
            jobs[job] = (label, omega)
        for job in concurrent.futures.as_completed(jobs):
            label, omega = jobs[job]
            block, records, endpoint = job.result()
            checks.append(dict(side=label, requested_omega=omega,
                               input_anchor=anchor, block=block, records=records))
            print("AXIAL VECTOR I120 " + json.dumps(dict(
                side=label, omega=omega, exit=block["exit"], terminal=block["terminal"],
                records=records)), flush=True)
    # Side repeats remain separate evidence; they do not silently overwrite the
    # nominal-window continuation record or restart a stopped axial direction.
    axial.STATE["i120_side_checks"] = sorted(checks, key=lambda value: value["side"])
    save_vector_notes()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--solver-sha", required=True)
    parser.add_argument("--coordinates", choices=("global",), default="global",
                        help="Frozen v4 diagnostic uses the global-coordinate model")
    parser.add_argument("--phase", choices=("diagnostic", "full"), default="diagnostic")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--block-size", type=int, default=20)
    parser.add_argument("--adaptive-key", default="adaptive")
    parser.add_argument("--vector-key", default="vector_tracking")
    parser.add_argument("--overlap-key", default="minimum_overlap")
    parser.add_argument("--full-approved", action="store_true")
    global ARGS, DRIVER_KEYS
    ARGS = parser.parse_args()
    ARGS.solver = ARGS.solver.resolve()
    if not re.fullmatch(r"[0-9a-fA-F]{64}", ARGS.solver_sha):
        parser.error("--solver-sha must be the confirmed SHA-256")
    ARGS.solver_sha = ARGS.solver_sha.casefold()
    if "v3" in ARGS.solver.name.casefold():
        parser.error("v3 is not allowed for long scans; await the confirmed leak-fixed v4")
    if not ARGS.solver.is_file() or driver.sha(ARGS.solver) != ARGS.solver_sha:
        parser.error("solver is absent or differs from the confirmed immutable hash")
    if ARGS.phase == "full" and not ARGS.full_approved:
        parser.error("full continuation requires explicit post-diagnostic --full-approved")
    DRIVER_KEYS = require_driver_keys(ARGS.adaptive_key, ARGS.vector_key, ARGS.overlap_key)
    print("AXIAL VECTOR CONFIG " + json.dumps(dict(
        solver=str(ARGS.solver), sha256=ARGS.solver_sha, driver_spec_keys=DRIVER_KEYS,
        options=OPTIONS, phase=ARGS.phase, coordinates=ARGS.coordinates, overlap=OVERLAP)), flush=True)
    # Reuse the existing audit/continuity/atomic-writing functions in this
    # process only. No original source, output, or full_driver file is edited.
    axial.OUTPUT = OUTPUT
    axial.SOLVER = ARGS.solver
    axial.APPROVED_SHA = ARGS.solver_sha
    axial.TAG = "vectors_" + ARGS.solver_sha[:8]
    axial.OPTIONS = OPTIONS
    axial.execute_block = execute_vector_block
    previous_argv = sys.argv
    try:
        sys.argv = [str(Path(__file__).resolve()), "--phase", ARGS.phase,
                    "--workers", str(ARGS.workers), "--block-size", str(ARGS.block_size)]
        axial.main()
    finally:
        sys.argv = previous_argv
    save_vector_notes()
    if ARGS.phase == "diagnostic":
        i120_side_checks()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


