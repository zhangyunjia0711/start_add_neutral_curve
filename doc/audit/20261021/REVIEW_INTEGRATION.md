# Neutral integration read-only review (2026-10-21)

Reviewed the first integration draft of `src/lst/module_neutral_curve.f90` and
`src/lst/module_lst.f90`, plus the grid/baseflow ownership, LST solver and SLEPc
adapter code needed to interpret those interfaces. No target source was edited.

## Required correction reported to main agent

In the first draft, `module_neutral_curve.f90:207-217` discards the global status
from `scan_roots` and sets the overall result to success whenever any roots were
returned. The numerical core deliberately returns valid partial roots alongside
an evaluation/iteration failure. It also preserves such failures without
bracketing across invalid points. Therefore the orchestration can wrongly return
success when a row is partial, or when later rows fail after an earlier valid row.
The main agent was asked to preserve hard failures separately and give them
priority over an earlier success. Each row should retain its existing diagnostic.
The analogous initial scan before task 36 can return partial roots and failure;
if the selected valid branch is still traced, this needs an explicit diagnostic
and a deliberate partial-result status policy.

## Additional observations reported

- `refine_peak` exits silently on callback failure. It should report that the
  ridge fell back to a sampled maximum rather than imply successful refinement.
- The final restoration solve is not checked. If that solve fails, the run can
  return success while its interactive eigenvalue is NaN. An actual workflow
  check should run a curve, then output the eigenfunction or solve the same point.
  Grid and baseflow objects themselves use allocatable deep copies; the local
  cleanup does not leave the solver pointing to freed profile storage.
- The within-row eigenvalue seed follows the latest valid sample, and bracket
  correction uses the closer endpoint. Across rows the seed comes from the
  previous ridge while evaluation restarts at the lower scan bound. Broad scans
  can therefore challenge mode following on the first evaluation. This is a
  physical continuation limitation, not a proven code defect; a narrow real
  benchmark and explicit rejected-evaluation reporting are appropriate.

The draft's dependent-eigenvalue exclusions and growth signs are consistent:
spatial scans exclude alpha and use `-Im(alpha)`; temporal scans exclude omega
and use `Im(omega)`. Bounds are passed as `bounds(:,1)=lower` and
`bounds(:,2)=upper`, matching the core. Task 36 rejects discrete station axes.
No other mandatory fix was established in the reviewed source snapshot.
