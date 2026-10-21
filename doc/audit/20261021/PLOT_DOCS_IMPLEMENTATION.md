# Plotting and documentation implementation record

Prepared only inside `neutral_curve_audit_20261021/staging`; the destination program was not edited by this subtask.

Files for integration:
- `tools/plot_neutral_curve.py`
- `tools/test_plot_neutral_curve.py`
- `doc/neutral_curve.md`

The parser validates the fourteen-column schema, handles multi-line VARIABLES and Fortran D exponents, and preserves ZONE row numbering and failed samples. Root lines are split by branch, missing row number and explicit failed records; arclength point order is retained. Growth plots use a masked rectilinear grid when every row has the same inner coordinate set. Local refinement clouds are triangulated with duplicate failures masked, only adjacent original rows joined, same-row edges restricted to sorted coordinate neighbors, and excessive coordinate gaps removed. Failure records without finite coordinates conservatively mask their row. An empty file parses successfully, while the CLI reports exit status 2 if neither input contains usable data.

All executed validation commands were run through `run_logged.py` into `commands.log`. Eleven regression tests passed on Python 3.8, NumPy 1.24.4 and Matplotlib 3.7.5:
1. Multi-line variables / D exponents / multi-zone branches and missing rows.
2. Empty and header-only files, useful CLI failure.
3. Rejected wrong schema and truncated records.
4. Rectilinear failed samples and completely empty rows.
5. Irregular failed samples and missing rows.
6. Irregular excessive outer gaps.
7. Refine evaluation order and duplicate failed samples.
8. Arclength turning point order.
9. Duplicate evaluations cannot falsely force a refinement cloud into a rectilinear grid.
10. Rectilinear excessive inner gaps.
11. Successful PNG export with roots and failed field samples.

`make_plot_preview.py` and `plot_preview/` are audit-only synthetic fixtures, not a physical LST result. The exported `synthetic_preview.png` was opened and visually inspected: axis labels, legend, colorbar and roots are readable, and both the field and root lines break across the failed outer row.

A first attempt to serialize source text through a quoted WSL command failed with a JSON escaping error before writing the file. This failure is in `commands.log`; subsequent source files were written via the file patch tool and all validation completed normally.

## Real START result verification

After the baseline upload and target integration were completed by the parent task, `START_ALL_IN_ONE/doc/neutral_curve.md` was updated directly and synchronized with the staging copy. No Fortran source was changed by this subtask.

Documentation now specifies the case-insensitive parameter names, continuation step default and per-parameter limits, nearest-predicted-target tracking only during neutral tasks, preserved partial output with batch failure exit status, full root collection, local/sampled peak interpretation, and executable sample configuration path and commands.

The completed real temporal task 37 result was parsed and plotted using the integrated target script:
- Data: `real_cases/task37/out/neutral_roots.plt` and `neutral_field.plt`.
- beta rows: -3.25, -3.20, -3.15; inner alpha range: 2.8 to 4.5.
- Two roots in each row; 140 evaluations, 545 internal solves, 0 rejected points, final result converged.
- Maximum absolute root growth: less than 6.5e-9.
- Plot saved as `real_cases/task37/out/neutral_preview.png`.
- The PNG was opened and visually checked: both root curves enclose the positive growth region; axes, legend and colorbar are legible.

The single-row predictor result at beta=-3.2059 was also inspected: alpha roots 2.986692664284999 and 3.957146023034614, with growth 8.282616152864068e-10 and -2.515682509477935e-11. These are numerical checks of the configured eigenmode, not proof of mesh convergence.


