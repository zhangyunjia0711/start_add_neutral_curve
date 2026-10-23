"""Prepare a focused external regression test; production source tests unchanged."""
import pathlib
ROOT=pathlib.Path(__file__).resolve().parent.parent;S=ROOT/'work/solver_source'
s=(S/'tests/adapter_safety_test.f90').read_text()
s=s.replace('program adapter_safety_test','program raw_vector_copy_test')
s=s.replace(' use mod_slepc_adapter\n',' use mod_slepc_adapter\n use mod_lst_solver_abstract, only: normal_u\n')
s=s.replace(' Vec :: vector\n',' Vec :: vector, raw, raw_saved, normalization_vector\n logical :: raw_valid\n PetscScalar :: raw_value\n real(REAL64) :: raw_error, raw_norm, saved_norm\n PetscLogDouble :: memory_before, memory_after\n integer :: memory_pass\n')
s=s.replace(' ! Before Solve and after an empty result, scalar getters must be safe.', ''' call eps_adapter%CopySelectedEigenVector(raw_value, raw, raw_error, raw_valid)
 if (raw_valid) stop 101
 call pep_adapter%CopySelectedEigenVector(raw_value, raw, raw_error, raw_valid)
 if (raw_valid) stop 102
 ! Before Solve and after an empty result, scalar getters must be safe.''')
for kind in ['eps','pep']:
    marker=' call '+kind+'_adapter%GetEigenVal(default_value)\n '+kind+'_adapter%track_target = .true.'
    block=''' call {k}_adapter%GetEigenVal(default_value)
 call {k}_adapter%CopySelectedEigenVector(raw_value, raw_saved, raw_error, raw_valid)
 if (.not. raw_valid .or. abs(raw_value-default_value)>1.e-12_REAL64) stop 103
 if (raw_error>1.e-8_REAL64 .or. raw_error/={k}_adapter%last_error) stop 104
 call VecNorm(raw_saved,NORM_2,saved_norm,ierr)
 if (ierr/=0) stop 105
 {k}_adapter%track_target = .true.
 ! Copy must not reselect a different eigenpair after the tracking flag changes.
 call {k}_adapter%CopySelectedEigenVector(raw_value, raw, raw_error, raw_valid)
 if (.not. raw_valid .or. abs(raw_value-default_value)>1.e-12_REAL64) stop 106
 call VecScale(raw,cmplx(2.0_REAL64,0.0_REAL64,REAL64),ierr)
 if (ierr/=0) stop 107
 call VecNorm(raw_saved,NORM_2,raw_norm,ierr)
 if (ierr/=0 .or. abs(raw_norm-saved_norm)>1.e-12_REAL64) stop 108
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 109
 {k}_adapter%track_target = .true.'''.format(k=kind)
    assert marker in s;s=s.replace(marker,block)
    marker=" print *, '"+kind.upper()+" default growth selection:'"
    block=''' ! A copied vector must survive adapter replacement during the next Solve.
 call VecNorm(raw_saved,NORM_2,raw_norm,ierr)
 if (ierr/=0 .or. abs(raw_norm-saved_norm)>1.e-12_REAL64) stop 110
 call VecDestroy(raw_saved,ierr)
 if (ierr/=0) stop 111
'''
    assert marker in s;s=s.replace(marker,block+marker)
marker=' call SlepcFinalize(ierr)'
block=''' ! Repeated normalization must release temporary Vec and VecScatter objects.
 call VecCreateSeq(PETSC_COMM_SELF,2000,normalization_vector,ierr)
 if (ierr/=0) stop 112
 call VecSet(normalization_vector,cmplx(1._REAL64,0._REAL64,REAL64),ierr)
 if (ierr/=0) stop 113
 do memory_pass=1,10
   call normal_u(normalization_vector)
 end do
 call PetscMallocGetCurrentUsage(memory_before,ierr)
 if (ierr/=0) stop 114
 do memory_pass=1,5000
   call normal_u(normalization_vector)
 end do
 call PetscMallocGetCurrentUsage(memory_after,ierr)
 if (ierr/=0) stop 115
 print *, 'NORMALIZATION MEMORY BEFORE/AFTER/DELTA:',memory_before,memory_after,memory_after-memory_before
 if (memory_after-memory_before>2.e6_REAL64) stop 116
 call VecDestroy(normalization_vector,ierr)
 if (ierr/=0) stop 117
 print *, 'RAW COPY AND NORMALIZATION MEMORY PASS'
'''
assert marker in s;s=s.replace(marker,block+marker)
p=ROOT/'work/test_raw_vector_copy.f90';assert not p.exists();p.write_text(s)
print(str(p),flush=True)