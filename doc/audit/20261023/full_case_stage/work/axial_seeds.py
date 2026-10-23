#!/usr/bin/env python3
"""Main-841-grid axial continuation seeds, with explicit gaps and atomic output."""
import argparse
import concurrent.futures
import copy
import datetime
import json
import math
from pathlib import Path
import re
import threading
import time

import full_driver as driver

ROOT = driver.ROOT
OUTPUT = ROOT / "outputs" / "axial_seeds.json"
SOLVER = ROOT / "work" / "solver_source" / "bin" / "START_LST"
APPROVED_SHA = "bcde0c5d8beaff8ca95f4c8ecc984272bf47ff6423c143fbb192387c491b8597"
OPTIONS = ["-pep_tol", "1e-12", "-pep_max_it", "1000",
           "-st_pc_factor_mat_solver_type", "superlu", "-pep_ncv", "8"]
NOMINAL = {"M1": (317, 8.17586, -0.98903),
           "M2": (312, 15.0, 15.0), "M3": (119, 0.0, -7.0)}
NUMBER = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"
PAIR = re.compile(r"^\s*(\d+)\s+\(\s*(" + NUMBER + r")\s*,\s*(" + NUMBER
                  + r")\s*\)\s+The error is\s+(" + NUMBER + r")", re.M)
LOCK = threading.Lock()
STATE = {}
X_GRID = []
TAG = APPROVED_SHA[:8]
EPS = 1e-8


def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=8))).isoformat()


def load_sources():
    catalogue = json.loads((ROOT / "outputs" / "seed_catalog.json").read_text())
    seeds = {value["mode"]: copy.deepcopy(value) for value in catalogue["baseline_converged_seeds"]}
    x_grid = catalogue["grids"]["main_841"]["x_nondim"]
    assert len(x_grid) == 841
    return seeds, x_grid


def persist(family, local, records):
    with LOCK:
        lookup = {(value["family"], value["iloc"]): value for value in STATE["records"]}
        for value in records:
            key = (value["family"], value["iloc"])
            # A later failed request must not overwrite an accepted actual seed.
            if key not in lookup or lookup[key]["status"] != 0 or value["status"] == 0:
                lookup[key] = value
        STATE["records"] = sorted(lookup.values(), key=lambda value: (value["family"], value["iloc"]))
        STATE["families"][family] = copy.deepcopy(local)
        STATE["updated"] = stamp()
        STATE["summary"] = {}
        for label in NOMINAL:
            valid = [value for value in STATE["records"] if value["family"] == label and value["status"] == 0]
            failed = [value["iloc"] for value in STATE["records"]
                      if value["family"] == label and value["status"] != 0]
            covered = {value["iloc"] for value in valid}
            STATE["summary"][label] = dict(accepted_stations=len(covered),
                failed_or_unverified_stations=failed,
                not_yet_covered=[station for station in range(1, 842) if station not in covered])
        temporary = OUTPUT.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(STATE, indent=2, allow_nan=False) + "\n")
        temporary.replace(OUTPUT)


def seed_dict(value):
    return {name: value[name] for name in ("omega_r", "beta_r", "alpha_r", "alpha_i")}


def audit_run(summary, family, lo, hi, previous, direction):
    directory = ROOT / "runs" / summary["spec"]["name"]
    source = directory / "out" / "neutral_field.plt"
    field = driver.table(source)
    log = (directory / "solver.log").read_text(errors="replace")
    candidates = [dict(alpha=complex(float(match[2].replace("D", "E")),
                                    float(match[3].replace("D", "E"))),
                       error=float(match[4].replace("D", "E")))
                  for match in PAIR.finditer(log)]
    error_lines = [line for line in log.splitlines()
                   if re.search(r"PETSC ERROR|SLEPC ERROR|Option left:|options you set that were not used", line, re.I)]
    fatal_lines = [line for line in log.splitlines()
                   if re.search(r"SIGSEGV|segmentation fault|Fortran runtime error|floating point exception", line, re.I)]
    global_bad = bool(summary.get("petsc_errors", 0) or summary.get("unused_options", 0)
                      or error_lines or fatal_lines or summary["exit"] == 124)
    by_station = {}
    for value in field:
        station = int(round(value["iloc"]))
        by_station.setdefault(station, []).append(value)
    provisional = {}
    for station in range(lo, hi + 1):
        values = by_station.get(station, [])
        reasons = []
        good = []
        for value in values:
            if value["status"] != 0 or not all(math.isfinite(value[key]) for key in
                    ("x", "omega_r", "beta_r", "alpha_r", "alpha_i", "growth")):
                continue
            if abs(value["omega_r"] - NOMINAL[family][1]) > EPS + 2e-12:
                continue
            if abs(value["beta_r"] - NOMINAL[family][2]) > 1e-10:
                continue
            good.append(value)
        if not values:
            reasons.append("no field record")
        if any(value["status"] != 0 for value in values):
            reasons.append("one or more failed field evaluations")
        if not good:
            reasons.append("no valid fixed-parameter candidate")
        if global_bad:
            reasons.append("global PETSc/options/runtime error in block")
        chosen = min(good, key=lambda value: abs(value["omega_r"] - NOMINAL[family][1])) if good else None
        selected_error = matching_distance = None
        if chosen is not None:
            actual = complex(chosen["alpha_r"], chosen["alpha_i"])
            if candidates:
                match = min(candidates, key=lambda value: abs(value["alpha"] - actual))
                matching_distance = abs(match["alpha"] - actual)
                selected_error = match["error"]
            if matching_distance is None or matching_distance > 1e-8:
                reasons.append("actual field eigenvalue has no matching printed PEP residual")
            elif selected_error > 1e-7:
                reasons.append("printed selected PEP residual exceeds 1e-7")
            if abs(chosen["x"] - X_GRID[station - 1]) > 1e-11 * max(1.0, abs(X_GRID[station - 1])):
                reasons.append("physical x differs from main-841 grid")
        record = dict(family=family, iloc=station, x=X_GRID[station - 1],
                      omega_r=chosen["omega_r"] if chosen else NOMINAL[family][1],
                      beta_r=chosen["beta_r"] if chosen else NOMINAL[family][2],
                      alpha_r=chosen["alpha_r"] if chosen else None,
                      alpha_i=chosen["alpha_i"] if chosen else None,
                      growth=chosen["growth"] if chosen else None,
                      status=0 if chosen is not None and not reasons else 2,
                      solver_sample_status=chosen["status"] if chosen else None,
                      source_path=str(source), source_row=int(chosen["row"]) if chosen else None,
                      printed_pep_error=selected_error, printed_pep_matching_distance=matching_distance,
                      nominal_omega=NOMINAL[family][1], omega_offset=(chosen["omega_r"] - NOMINAL[family][1]) if chosen else None,
                      reasons=reasons, case_template="M1 main_841", direction=direction)
        provisional[station] = record
    traversal = list(range(lo, hi + 1)) if direction != "up" else list(range(hi, lo - 1, -1))
    anchor = previous
    first_gap = None
    jumps = []
    for station in traversal:
        value = provisional[station]
        if first_gap is not None:
            if value["status"] == 0:
                value["status"] = 3
                value["reasons"].append("actual sample exists beyond unresolved gap; family continuity unverified")
            continue
        if value["status"] != 0:
            first_gap = station
            continue
        if anchor is not None and anchor["iloc"] != station:
            if abs(anchor["iloc"] - station) != 1:
                value["status"] = 3
                value["reasons"].append("not adjacent to the last accepted station")
                first_gap = station
                continue
            jump = abs(complex(value["alpha_r"] - anchor["alpha_r"],
                               value["alpha_i"] - anchor["alpha_i"]))
            value["alpha_jump_from_adjacent"] = jump
            jumps.append(jump)
            if jump > 2:
                value["status"] = 3
                value["reasons"].append("adjacent alpha jump exceeds configured 2.0")
                first_gap = station
                continue
        anchor = value
    result = dict(name=summary["spec"]["name"], direction=direction, lo=lo, hi=hi,
                  summary_path=str(directory / "summary.json"), exit=summary["exit"],
                  terminal=summary["terminal"], seconds=summary["seconds"],
                  evaluations=summary["evaluations"], solves=summary["solves"],
                  rejected=summary["rejected"], failed_field_records=summary["failed_field_records"],
                  petsc_errors=summary.get("petsc_errors"), unused_options=summary.get("unused_options"),
                  error_lines=error_lines, fatal_lines=fatal_lines,
                  raw_pep_error_max=max((value["error"] for value in candidates), default=None),
                  selected_pep_error_max=max((value["printed_pep_error"] for value in provisional.values()
                                              if value["printed_pep_error"] is not None), default=None),
                  max_adjacent_alpha_jump=max(jumps, default=0), first_unresolved_station=first_gap,
                  accepted_stations=[value["iloc"] for value in provisional.values() if value["status"] == 0],
                  failed_or_unverified=[value["iloc"] for value in provisional.values() if value["status"] != 0])
    return result, list(provisional.values()), anchor


def execute_block(family, lo, hi, anchor, direction):
    name = f"axial_{family}_{direction}_i{lo:03d}_{hi:03d}_{TAG}"
    spec = dict(name=name, solver=str(SOLVER), mode="M1", family=family,
                station=lo, seed_station=anchor["iloc"], task=33,
                outer_parameter="iloc", inner_parameter="omega",
                outer_min=lo, outer_max=hi, outer_points=hi - lo + 1,
                inner_min=NOMINAL[family][1] - EPS, inner_max=NOMINAL[family][1] + EPS,
                inner_points=2, step=.2, jump=2, seed=seed_dict(anchor),
                options=OPTIONS, timeout=900)
    return audit_run(driver.run(spec), family, lo, hi, anchor, direction)


def initial_seed(family, baseline):
    station, omega, beta = NOMINAL[family]
    anchor = dict(family=family, iloc=station, x=X_GRID[station - 1],
                  omega_r=omega, beta_r=beta,
                  alpha_r=baseline["alpha_r"], alpha_i=baseline["alpha_i"])
    block, records, endpoint = execute_block(family, station, station, anchor, "seed")
    actual = records[0]
    if actual["status"] == 0:
        distance = abs(complex(actual["alpha_r"] - baseline["alpha_r"],
                               actual["alpha_i"] - baseline["alpha_i"]))
        relative = distance / max(1.0, abs(complex(baseline["alpha_r"], baseline["alpha_i"])))
        actual["mapped_seed_alpha_distance"] = distance
        actual["mapped_seed_alpha_relative_distance"] = relative
        if family == "M2" and (distance > 1.0 or relative > .2):
            actual["status"] = 3
            actual["reasons"].append("mapped M2 seed differs substantially; family stopped for review")
    block["mapped_seed_confirmed"] = actual["status"] == 0
    return block, actual


def family_worker(family, baseline, phase, block_size):
    with LOCK:
        local = copy.deepcopy(STATE["families"].get(family, {}))
    seed_station = NOMINAL[family][0]
    if not local.get("seed"):
        block, seed = initial_seed(family, baseline)
        local = dict(seed=seed, seed_provenance=baseline, blocks=[block],
                     directions={}, diagnostic_complete=False)
        persist(family, local, [seed])
        print("AXIAL INITIAL " + json.dumps(dict(family=family, station=seed_station,
              alpha_r=seed["alpha_r"], alpha_i=seed["alpha_i"],
              status=seed["status"], mapped_distance=seed.get("mapped_seed_alpha_distance"),
              reasons=seed["reasons"])), flush=True)
        if seed["status"] != 0:
            local["stopped"] = "initial main-grid seed not verified"
            persist(family, local, [])
            return
    if local["seed"]["status"] != 0:
        return
    for direction in ("down", "up"):
        if direction not in local["directions"]:
            local["directions"][direction] = dict(anchor=copy.deepcopy(local["seed"]),
                next_station=seed_station + (1 if direction == "down" else -1),
                stopped=False, completed_blocks=0)
        state = local["directions"][direction]
        count = 0
        while not state["stopped"] and 1 <= state["next_station"] <= 841:
            if phase == "diagnostic" and state["completed_blocks"] >= 1:
                break
            start = state["next_station"]
            lo, hi = ((start, min(841, start + block_size - 1)) if direction == "down"
                      else (max(1, start - block_size + 1), start))
            block, records, anchor = execute_block(family, lo, hi, state["anchor"], direction)
            local["blocks"].append(block)
            state["completed_blocks"] += 1
            if anchor is not None:
                state["anchor"] = anchor
            if block["first_unresolved_station"] is not None:
                state["stopped"] = True
                state["stop_reason"] = "unresolved station; continuation beyond it is not assumed"
            else:
                state["next_station"] = hi + 1 if direction == "down" else lo - 1
            persist(family, local, records)
            print("AXIAL BLOCK " + json.dumps(dict(family=family, direction=direction,
                  range=[lo, hi], seconds=block["seconds"], exit=block["exit"],
                  terminal=block["terminal"], accepted=len(block["accepted_stations"]),
                  failures=block["failed_or_unverified"], alpha_jump=block["max_adjacent_alpha_jump"],
                  selected_pep_error_max=block["selected_pep_error_max"])), flush=True)
            count += 1
            if phase == "diagnostic":
                break
        if not 1 <= state["next_station"] <= 841:
            state["finished_native_boundary"] = True
    local["diagnostic_complete"] = all(
        local["directions"][direction]["completed_blocks"] >= 1
        or local["directions"][direction]["stopped"] for direction in ("down", "up"))
    persist(family, local, [])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase", choices=("diagnostic", "full"), default="diagnostic")
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--block-size", type=int, default=20)
    args = parser.parse_args()
    if not 1 <= args.workers <= 3 or args.block_size < 1:
        parser.error("workers must be 1..3 and block-size positive")
    assert SOLVER.exists(), "Validated independent solver is not yet built"
    assert driver.sha(SOLVER) == APPROVED_SHA, "Solver differs from the confirmed build"
    validation = json.loads((ROOT / "outputs" / "solver_build_validation.json").read_text())
    assert validation["all_benchmarks_accepted"] and validation["initial_target"]["accepted"]
    seeds, grid = load_sources()
    global X_GRID, STATE
    X_GRID = grid
    if OUTPUT.exists():
        STATE = json.loads(OUTPUT.read_text())
        assert STATE["solver_sha256"] == APPROVED_SHA
    else:
        STATE = dict(created=stamp(), solver=str(SOLVER), solver_sha256=APPROVED_SHA,
                     options=OPTIONS, grid="main_841 via M1 template for all three families",
                     records=[], families={}, notes=[
                         "Actual axial continuation seeds, not neutral roots or a complete spectrum.",
                         "Inner omega lies within nominal +/-1e-8; each record retains actual omega.",
                         "Only consecutive verified stations receive status=0; failed/gap-crossed samples remain nonzero.",
                         "M2 archived native-grid alpha is only a mapped initial guess until the new main-grid solve.",
                         "Only three original families are solved; common conjugates may be derived algebraically."])
    if args.phase == "full":
        assert all(STATE["families"].get(family, {}).get("diagnostic_complete")
                   for family in NOMINAL), "Complete and review all first-block diagnostics before full continuation"
    started = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        jobs = [executor.submit(family_worker, family, seeds[family], args.phase, args.block_size)
                for family in NOMINAL]
        for job in concurrent.futures.as_completed(jobs):
            job.result()
    print("AXIAL SUMMARY " + json.dumps(dict(phase=args.phase, seconds=round(time.time() - started, 3),
          summary={family:dict(accepted_stations=value["accepted_stations"],
                              failed_or_unverified_stations=value["failed_or_unverified_stations"],
                              not_yet_covered_count=len(value["not_yet_covered"]))
                   for family, value in STATE["summary"].items()})), flush=True)


if __name__ == "__main__":
    main()

