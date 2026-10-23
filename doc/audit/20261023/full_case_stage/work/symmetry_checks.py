#!/usr/bin/env python3
"""Fresh paired solves for the real-coefficient common-conjugation symmetry."""
import concurrent.futures
import configparser
import datetime
import json
import math
from pathlib import Path
import re
import subprocess
import sys
import time

import full_driver as driver

ROOT = driver.ROOT
OUTPUT = ROOT / "outputs" / "symmetry_checks.json"
POINTS = [
    dict(id=1, mode="M3", station=119, omega=0.0, beta=-7.0,
         alpha_r=51.76909439, alpha_i=-0.55946676),
    dict(id=2, mode="M1", station=317, omega=8.17586, beta=-0.98903,
         alpha_r=30.87103096, alpha_i=-0.8416318),
    dict(id=3, mode="M3", station=599, omega=-5.0, beta=-7.0,
         alpha_r=5.245889049789906, alpha_i=-0.04450596796116299),
    dict(id=4, mode="M1", station=321, omega=22.61994307483593,
         beta=10.57775824698773, alpha_r=44.02144804943075,
         alpha_i=5.812488106864006e-8),
]
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
PAIR_ERROR = re.compile(
    r"^\s*(\d+)\s+\(\s*(" + NUMBER + r")\s*,\s*(" + NUMBER
    + r")\s*\)\s+The error is\s+(" + NUMBER + r")", re.M)
TOLERANCE = 1e-8


def numeric(value):
    return float(value.replace("D", "E").replace("d", "e"))


def local_result(path):
    records = []
    if not path.exists():
        return None
    for raw in path.read_text(errors="replace").splitlines():
        try:
            values = [numeric(token) for token in raw.split()]
        except ValueError:
            continue
        if len(values) == 13:
            records.append(dict(station=int(round(values[0])), x=values[1],
                                omega_r=values[5], omega_i=values[6],
                                alpha_r=values[7], alpha_i=values[8],
                                beta_r=values[9], beta_i=values[10]))
    return records[-1] if records else None


def summarize(directory, exit_code, seconds):
    log_path = directory / "solver.log"
    log = log_path.read_text(errors="replace")
    local = local_result(directory / "out" / "LST_LOCAL_eigenvalues.plt")
    candidates = [dict(index=int(match[1]), alpha_r=numeric(match[2]),
                       alpha_i=numeric(match[3]), pep_error=numeric(match[4]))
                  for match in PAIR_ERROR.finditer(log)]
    chosen = None
    if local is not None and candidates:
        reference = complex(local["alpha_r"], local["alpha_i"])
        chosen = min(candidates, key=lambda item: abs(complex(item["alpha_r"], item["alpha_i"]) - reference))
        chosen = dict(chosen, table_match_distance=abs(complex(chosen["alpha_r"], chosen["alpha_i"]) - reference))
    petsc_errors = [line for line in log.splitlines()
                    if "PETSC ERROR" in line.upper() or "SLEPC ERROR" in line.upper()]
    return dict(directory=str(directory), exit=exit_code, seconds=round(seconds, 3),
                local_table_record=local, selected_raw_eigenpair=chosen,
                raw_pep_eigenpairs=candidates, petsc_errors=petsc_errors,
                solver_log=str(log_path), config=str(directory / "config.ini"))


def single_fallback(spec):
    name = spec["name"].replace("_neutral", "_single")
    directory = ROOT / "runs" / name
    directory.mkdir(parents=True, exist_ok=True)
    summary_path = directory / "symmetry_summary.json"
    if summary_path.exists():
        raise RuntimeError("Refusing to reuse a previous actual endpoint solve: " + str(directory))
    cfg = configparser.ConfigParser()
    cfg.read(ROOT / "runs" / spec["name"] / "config.ini")
    cfg["neutral"]["enabled"] = "F"
    with (directory / "config.ini").open("w") as stream:
        cfg.write(stream)
    command = [str(driver.SOLVER), "-f", "config.ini"] + list(driver.DEFAULT_OPTIONS)
    started = time.time()
    with (directory / "solver.log").open("w") as stream:
        stream.write("COMMAND " + json.dumps(command) + "\n")
        stream.write("STDIN: 0\\n (exit the interactive menu after the local solve)\n")
        stream.flush()
        try:
            process = subprocess.run(command, cwd=directory, env=driver.ENV,
                                     input="0\n", text=True, stdout=stream,
                                     stderr=subprocess.STDOUT, timeout=300)
            code = process.returncode
        except subprocess.TimeoutExpired:
            code = 124
    summary = summarize(directory, code, time.time() - started)
    summary["neutral_enabled"] = False
    summary["fresh_execution"] = True
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def endpoint(point, flip):
    sign = -1 if flip else 1
    side = "flip" if flip else "original"
    seed = dict(omega_r=sign * point["omega"], beta_r=sign * point["beta"],
                alpha_r=sign * point["alpha_r"], alpha_i=point["alpha_i"])
    prefix = "sym_m1" if point["id"] == 4 else "sym"
    spec = dict(name=f"{prefix}_{point['id']:02d}_{side}_neutral", mode=point["mode"],
                station=point["station"], task=33, omega_min=seed["omega_r"],
                omega_max=seed["omega_r"], omega_points=1, beta_min=seed["beta_r"],
                beta_max=seed["beta_r"], beta_points=1, seed=seed, timeout=300)
    directory = ROOT / "runs" / spec["name"]
    if (directory / "summary.json").exists():
        raise RuntimeError("Refusing cached neutral endpoint: " + str(directory))
    attempt = driver.run(spec)
    attempt_log = (directory / "out" / "neutral_run.log")
    attempt_text = attempt_log.read_text(errors="replace") if attempt_log.exists() else ""
    # A one-sample neutral scan is rejected by the current numerical scanner.
    # Retain that output but independently solve the local problem with neutral
    # disabled, rather than treating the supplied alpha seed as a measurement.
    needs_fallback = (attempt["exit"] != 0 or attempt["terminal"] == "MISSING"
                      or attempt["field_records"] == 0 or "invalid" in attempt["terminal"].casefold()
                      or "invalid" in attempt_text.casefold())
    if needs_fallback:
        actual = single_fallback(spec)
    else:
        actual = summarize(directory, attempt["exit"], attempt["seconds"])
        actual["neutral_enabled"] = True
        actual["fresh_execution"] = True
    actual["point_id"] = point["id"]
    actual["side"] = side
    actual["input_seed"] = seed
    actual["seed_real_shift"] = 1e-3
    actual["requested_station"] = point["station"]
    actual["neutral_attempt"] = attempt
    attach_target_candidate(actual)
    print(json.dumps(dict(point=point["id"], side=side, exit=actual["exit"],
                          measured=actual["selected_raw_eigenpair"],
                          petsc_error_lines=len(actual["petsc_errors"]),
                          fallback=needs_fallback)), flush=True)
    return actual


def attach_target_candidate(actual):
    seed = actual["input_seed"]
    target = complex(seed["alpha_r"], seed["alpha_i"])
    candidates = actual["raw_pep_eigenpairs"]
    closest = min(candidates, key=lambda item: abs(complex(item["alpha_r"], item["alpha_i"]) - target)) if candidates else None
    if closest is not None:
        closest = dict(closest, distance_to_requested_alpha=abs(
            complex(closest["alpha_r"], closest["alpha_i"]) - target))
    actual["nearest_requested_seed_candidate"] = closest
    actual["requested_seed_reproduced"] = (
        closest is not None and closest["distance_to_requested_alpha"] <= 1e-4)


def existing_endpoint(point, flip):
    side = "flip" if flip else "original"
    sign = -1 if flip else 1
    directory = ROOT / "runs" / f"sym_{point['id']:02d}_{side}_single"
    actual = json.loads((directory / "symmetry_summary.json").read_text())
    actual.update(point_id=point["id"], side=side,
                  input_seed=dict(omega_r=sign * point["omega"], beta_r=sign * point["beta"],
                                  alpha_r=sign * point["alpha_r"], alpha_i=point["alpha_i"]),
                  seed_real_shift=1e-3, requested_station=point["station"])
    actual["neutral_attempt"] = json.loads(
        (ROOT / "runs" / f"sym_{point['id']:02d}_{side}_neutral" / "summary.json").read_text())
    attach_target_candidate(actual)
    return actual


def main():
    results = []
    points_to_solve = POINTS
    if "--recheck-fourth" in sys.argv:
        results = [existing_endpoint(point, flip) for point in POINTS[:3] for flip in (False, True)]
        points_to_solve = POINTS[3:]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        jobs = [executor.submit(endpoint, point, flip) for point in points_to_solve for flip in (False, True)]
        for job in concurrent.futures.as_completed(jobs):
            results.append(job.result())
    checks = []
    for point in POINTS:
        original = next(item for item in results if item["point_id"] == point["id"] and item["side"] == "original")
        flipped = next(item for item in results if item["point_id"] == point["id"] and item["side"] == "flip")
        left, right = original["selected_raw_eigenpair"], flipped["selected_raw_eigenpair"]
        verified = all(item["exit"] == 0 and not item["petsc_errors"]
                       and item["selected_raw_eigenpair"] is not None
                       and item["selected_raw_eigenpair"]["table_match_distance"] <= 2e-7
                       and item["requested_seed_reproduced"]
                       for item in (original, flipped))
        residual_r = residual_i = residual_complex = None
        if left is not None and right is not None:
            residual_r = abs(left["alpha_r"] + right["alpha_r"])
            residual_i = abs(left["alpha_i"] - right["alpha_i"])
            residual_complex = abs(complex(right["alpha_r"], right["alpha_i"])
                                   + complex(left["alpha_r"], left["alpha_i"]).conjugate())
            verified = verified and residual_complex <= TOLERANCE
        checks.append(dict(point=point, original=original, flipped=flipped,
                           relation="alpha(-omega,-beta) = -conj(alpha(omega,beta))",
                           alpha_real_antisymmetry_error=residual_r,
                           alpha_imag_symmetry_error=residual_i,
                           complex_relation_error=residual_complex,
                           absolute_tolerance=TOLERANCE, passed=bool(verified)))
    summary = dict(
        generated=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat(),
        description="Both endpoints of every pair are fresh actual START solves; input alpha values are only seeds.",
        relation="(omega,beta,alpha_r,alpha_i) -> (-omega,-beta,-alpha_r,alpha_i)",
        options=driver.DEFAULT_OPTIONS, seed_shift=1e-3, maximum_parallel_processes=4,
        pep_error_note="Raw values printed by the unmodified solver as 'The error is'; no error-type reinterpretation.",
        checks=checks, all_passed=all(check["passed"] for check in checks))
    invalid_attempts = []
    for side in ("original", "flip"):
        directory = ROOT / "runs" / f"sym_04_{side}_single"
        if (directory / "symmetry_summary.json").exists():
            invalid_attempts.append(dict(
                reason="invalid case/grid: M2 661-station template used for a required M1 841-station point",
                excluded_from_checks=True,
                actual=json.loads((directory / "symmetry_summary.json").read_text())))
    summary["invalid_case_grid_attempts"] = invalid_attempts
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(dict(output=str(OUTPUT), all_passed=summary["all_passed"],
                          errors=[check["complex_relation_error"] for check in checks])), flush=True)
    return 0 if summary["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

