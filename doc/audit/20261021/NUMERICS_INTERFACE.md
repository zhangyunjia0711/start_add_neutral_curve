# Standalone neutral numerics (2026-10-21)

`module_neutral_numerics.f90` is PETSc/SLEPc independent Fortran 2008, using
`iso_fortran_env:real64`. It contains no file I/O and does not terminate the caller.
The caller supplies a continuously tracked eigenmode through this callback:

```fortran
subroutine evaluate(point,seed,eigen,growth,ok)
  real(dp), intent(in) :: point(2)
  complex(dp), intent(in) :: seed
  complex(dp), intent(out) :: eigen
  real(dp), intent(out) :: growth
  logical, intent(out) :: ok
end subroutine
```

`point=(outer coordinate,inner coordinate)`. `options%bounds(:,1)` contains
the lower bounds and `options%bounds(:,2)` contains the upper bounds.
Spatial growth is normally `-aimag(alpha)` and temporal growth is normally
`aimag(omega)`; the callback defines that physical convention.

```fortran
call scan_roots(evaluate,fixed_point,axis,lower,upper,seed,options,roots,status)
call trace_curve(evaluate,start_point,seed,direction,options,points,status,initial_tangent)
```

`roots` and `points` are allocatable arrays of `neutral_sample`. Fields are
`point(2),eigen,growth,residual,status,iterations`; `residual=abs(growth)`.
Every returned sample has finite eigenvalue/growth and `residual<=f_tol`.
The global status reports a partial scan or trace failure even if valid samples
were already returned. No failed point is written as a converged neutral point.
Callback failures or NaN/Inf values split scan intervals; no bracket spans them.

`scan_roots` changes one coordinate and scans all `scan_samples` locations over
the finite specified interval. An accepted sampled root or sign-changing interval
can yield a result. Each bracket is corrected with bounded secant/bisection,
or pure bisection if `options%use_bisection=.true.`. Width alone never certifies
convergence. All converged roots are retained in ascending coordinate order;
adjacent detections within `x_tol` are deduplicated. A finite scan cannot guarantee
detection of tangential roots or multiple roots between adjacent scan locations.
An unsampled tangency is therefore not reported as a proven root.

`trace_curve` first projects a nonneutral start with damped Newton updates.
It then estimates the growth gradient, predicts along the tangent, and corrects
along the normal through the predictor. This enforces the pseudo-arclength
constraint and passes coordinate turning points. Correction failure halves the
step down to `ds_min`. Normal correction stays within twice the current step
from the predictor. `direction=+1/-1` controls orientation. The optional
`initial_tangent(2)` selects the initial tangent orientation, then `direction`
applies its sign; a vector normal to the curve cannot select orientation.
Gradients use centered finite differences, with one-sided differences only at
the configured boundary. Callback failures are not silently used as one-sided
derivative data. The initial zero gradient or a singular normal derivative fails
explicitly. `ds`, finite differences and coordinates use the caller's scales.

Closure requires enough points, traveled distance, proximity to the initial
point and matching tangent orientation. A closed trace ends with an exact copy
of the already converged initial sample. Both endpoints therefore satisfy the
growth tolerance. The trace ends on closure, parameter boundary, evaluation
failure, inability to correct, or `max_points`; loops and corrections are bounded.

Public constants (also accepted by `neutral_status_name`):

| Constant | Value | Meaning |
|---|---:|---|
| NC_SUCCESS | 0 | Converged scan/point |
| NC_NO_ROOT | 1 | No sampled/bracketed root found |
| NC_EVALUATION_FAILED | 2 | Callback failed or returned nonfinite data |
| NC_MAX_ITERATIONS | 3 | Root/correction iteration or trace point limit |
| NC_BOUNDARY | 4 | Parameter bound reached |
| NC_CLOSED | 5 | Closed curve completed |
| NC_STEP_TOO_SMALL | 6 | Further trace step cannot progress |
| NC_INVALID_INPUT | 7 | Invalid controls, coordinates, seed, or tangent |
| NC_DEGENERATE | 8 | Zero gradient or singular correction direction |

Analytic validation uses the unstable interior of
`growth=1-x*x/4-y*y`: two roots, a complete closed ellipse and coordinate turning
points, both trace directions, nonneutral starting point and prescribed initial
tangent. Other cases cover NaN and failed scan gaps, unsampled tangency, four
disjoint roots, discontinuous sign jump with no actual zero, constant growth,
nonfinite eigenvalue, parameter bounds, point limits and trace failures.
These tests validate algorithms and interfaces, not physical LST eigenvalues.

```sh
gfortran -std=f2008 -Wall -Wextra -Wimplicit-interface -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -g module_neutral_numerics.f90 \
  test_neutral_numerics.f90 -o test_neutral_numerics
./test_neutral_numerics
```
