# SLEPc adapter safety staging changes

- Scope: only a staging copy of `src/public/module_slepc_adapter.f90`; original source remains untouched.
- Both EPS and PEP destroy all old PETSc Vec objects before reallocating result handles and on Finalize. This fixes repeated-solve object leaks without changing public method signatures.
- Solve now fills converged eigenvalues and vectors immediately. Previous GetEigenVal-only callers otherwise left vectors uninitialized, even though the next Solve copied one as its initial space.
- A repeated solve reuses a prior initial vector only when nonempty and dimensionally compatible. Its temporary duplicate is destroyed after Solve.
- GetEigenVal and GetEigenVector return a complex quiet NaN when nothing converged; the latter also returns PETSC_NULL_VEC. Array getters return allocated empty arrays, avoiding index-zero access and sorting of empty data.
- Added public type component `last_error`, initialized/reset to huge. EPS/PEP scalar and selected-vector getters record the selected eigenpair's absolute error using EPS_ERROR_ABSOLUTE / PEP_ERROR_ABSOLUTE. A zero-convergence result retains huge.
- Other solvers: every existing EPS/PEP client sees the same eigenpair selection and method signatures; all benefit from vector lifecycle cleanup. Eager extraction adds one eigenpair/vector read per converged result; previous repeated solves using uninitialized vectors may show different iteration paths. Returned vector handles remain borrowed from the solver adapter, valid until its next Solve or Finalize.

Validation delegated to the root agent's isolated integration build. No target source was changed by this subtask.

## Isolated verification

- Staging adapter compiled successfully with the existing MPICH/GFortran/PETSc/SLEPc headers and bounds checks; output object/module stayed in staging.
- `adapter_safety_test.f90` passed (exit 0): unallocated/zero-length EPS and PEP scalar/vector/array getters; five successive EPS solves and five successive PEP solves where only GetEigenVal runs between solves; final single-vector queries; selected absolute residual <= 1e-8; vector cleanup on Finalize. No PETSc object errors or leak warning appeared.
- All Linux compile/run commands and their output are appended to commands.log through run_logged.py. The test uses a standalone diagonal generalized eigenproblem and polynomial eigenproblem, without modifying or running the target working tree.


## Continuation branch selection follow-up

- Added public `track_target = .false.` in slepc_type. When explicitly enabled, single-value and single-vector EPS/PEP getters select the converged eigenpair with minimum `abs(kr - target)`. Their evIndex arrays place that selected pair last, so both query APIs agree.
- With track_target disabled, the existing growth-rate sorting remains unchanged. Multi-eigenpair array getter sorting also remains unchanged.
- Both Solve methods now reuse the last selected evIndex vector as the next initial space, rather than the last raw SLEPc convergence index. The helper checks allocation, nonzero size, matching array length, and index bounds, then falls back to vector 1 when no valid selection exists.
- Updated only target `src/public/module_slepc_adapter.f90` and its audit staging copy; the root agent manages LST wrappers and scan option lifecycle.

Tracking verification passed in both spatial (mode 2) and temporal (mode 1) diagonal test configurations. In temporal mode EPS default selects approximately 14+0.014i while target tracking at 2.01 selects 2+0.002i; PEP default selects approximately -1.010205-0.001021i while target tracking at -2.01 selects -2.041685-0.002085i. Single-vector queries agree with the tracked value; one further solve after selecting that vector retains the same branch. Existing empty-result, repeated-solve, and selected-residual checks also pass (exit 0). All compilation and execution logs are retained in commands.log.
