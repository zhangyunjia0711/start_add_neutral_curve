"""Extend the focused external regression with explicit non-nearest MAC selection."""
import pathlib
from full_driver import ROOT
p=ROOT/'work/test_raw_vector_copy.f90';s=p.read_text()
assert 'EXPLICIT MAC SELECT' not in s
s=s.replace(' real(REAL64) :: raw_error, raw_norm, saved_norm',' real(REAL64) :: raw_error, raw_norm, saved_norm, selected_overlap')
for k in ['eps','pep']:
    marker=' call '+k+'_adapter%GetEigenVector(value, vector)\n if (abs(value-tracked_value)'
    insertion=''' ! EXPLICIT MAC SELECT: choose the reference mode even when the scalar
 ! nearest-target getter selected another pair. Gates are not relaxed.
 call {k}_adapter%SelectEigenVectorByOverlap(raw_saved,default_value, &
      2.0_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 121
 if (selected_overlap<.999999_REAL64 .or. raw_error/={k}_adapter%last_error) stop 122
 if ({k}_adapter%overlap_selected_index<=0) stop 123
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 124
 ! The only candidate inside this tight lambda gate has a different vector.
 call {k}_adapter%SelectEigenVectorByOverlap(raw_saved,tracked_value, &
      1.e-5_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (raw_valid) stop 125
 call {k}_adapter%CopySelectedEigenVector(raw_value,raw,raw_error,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 126
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 127
 ! No candidate exists inside a remote tight gate; selection must stay intact.
 call {k}_adapter%SelectEigenVectorByOverlap(raw_saved,cmplx(1000._REAL64,0._REAL64,REAL64), &
      1.e-5_REAL64,1.e-8_REAL64,.95_REAL64,raw_value,raw,raw_error,selected_overlap,raw_valid)
 if (raw_valid) stop 128
 call {k}_adapter%CopySelectedEigenVector(raw_value,raw,raw_error,raw_valid)
 if (.not.raw_valid .or. abs(raw_value-default_value)>1.e-10_REAL64) stop 129
 call VecDestroy(raw,ierr)
 if (ierr/=0) stop 130
 ! Restore the original scalar query to check its unchanged default policy.
 call {k}_adapter%GetEigenVal(tracked_value)
'''.format(k=k)
    assert marker in s;s=s.replace(marker,insertion+marker)
p.write_text(s);print('Focused external MAC selection regression prepared; not executed until shared source is ready.',flush=True)