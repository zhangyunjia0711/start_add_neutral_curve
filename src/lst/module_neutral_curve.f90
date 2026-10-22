! Neutral-curve orchestration for START's direct LST solver.
#include <petsc/finclude/petsc.h>
module mod_neutral_curve
  use global_parameters, mpi_size_unused => size
  use mod_lst_solver_abstract, only: type_lst_abstract
  use mod_lst_grid_profile, only: type_lst_grid_profile, assemble_lst_grid_profile
  use mod_lst_baseflow_profile, only: type_lst_bf_profile, assemble_lst_baseflow_profile
  use mod_lst_eqn, only: lst_eigen_error, lst_track_target, lst_is_tracking, lst_copy_selected_eigenvector, &
    lst_select_eigenvector_by_overlap, lst_prepare_candidate_retry
  use mod_neutral_numerics, only: neutral_options, neutral_sample, scan_roots, trace_curve, &
    neutral_status_name, NC_SUCCESS, NC_NO_ROOT, NC_EVALUATION_FAILED, NC_INVALID_INPUT, NC_CLOSED
  use petsc
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: run_neutral_curve, neutral_batch_enabled, neutral_exit_code
  integer :: neutral_exit_code=0
#include "../include/start.f90"
contains
  logical function neutral_batch_enabled()
    neutral_batch_enabled=.false.
    if (cfg%has_key('neutral','enabled')) call cfg%get('neutral','enabled',neutral_batch_enabled)
  end function

  subroutine run_neutral_curve(lst,ome,alf,bat,task_override,result_status)
    class(type_lst_abstract), intent(inout) :: lst
    complex(real64), intent(inout) :: ome,alf,bat
    integer, intent(in), optional :: task_override
    integer, intent(out), optional :: result_status
    type(neutral_options) :: opt
    type(neutral_sample), allocatable :: roots(:), points(:), history(:)
    Vec, allocatable :: history_vectors(:)
    logical, allocatable :: history_has_vector(:)
    type(neutral_sample) :: best
    type(type_lst_grid_profile) :: grid
    type(type_lst_bf_profile) :: flow
    complex(real64) :: saved_ome,saved_alf,saved_bat,seed,next_seed,row_seed,sample_eigen,row_anchor_eigen
    real(real64) :: sample_growth
    real(real64) :: saved_re,center(2),p(2),eigen_tol,jump_tol,theta,initial_tangent(2)
    real(real64) :: outer_lo,outer_hi,inner_lo,inner_hi,row_value,continuation_step,minimum_step,minimum_overlap
    integer :: saved_iloc,outer_id,inner_id,task,nrows,row,n,status,overall,branch,comm_size,ios
    integer :: roots_unit,field_unit,log_unit,ridge_unit,eval_count,failed_count,solve_count,nhistory,hard_failure
    integer :: retry_eigenpairs
    logical :: ok,best_valid,have_outputs,has_tangent,saved_tracking,adaptive_steps,vector_tracking,sample_only,sample_valid
    logical :: reverse_outer,reverse_inner,anchor_row_seed,have_row_anchor
    character(32) :: outer_name,inner_name
    character(256) :: prefix

    overall=NC_INVALID_INPUT
    if (present(result_status)) result_status=overall
    call MPI_Comm_size(PETSC_COMM_WORLD,comm_size,ierr)
    if (comm_size/=1) then
      print *, 'Neutral curves currently require a single MPI process.'
      return
    end if
    if (lst_IsAdjoint) then
      print *, 'Neutral curves require direct LST (adjoint=F).'
      return
    end if
    task=37
    call read_int('task',task)
    if (present(task_override)) task=task_override
    if (.not.any(task==[31,32,33,36,37])) then
      print *, '[neutral]task must be 31,32,33,36 or 37.'
      return
    end if
    saved_ome=ome; saved_alf=alf; saved_bat=bat
    saved_re=Re; saved_iloc=iloc
    if (lst_mode==LST_TEMPORAL_MODE) then
      inner_name='alpha'; seed=ome
      if (abs(aimag(alf))+abs(aimag(bat))>1.e-12_real64) then
        print *, 'Temporal neutral curves require real alpha and beta.'
        return
      end if
    else
      inner_name='omega'; seed=alf
      if (abs(aimag(ome))+abs(aimag(bat))>1.e-12_real64) then
        print *, 'Spatial neutral curves require real omega and beta.'
        return
      end if
    end if
    outer_name='beta'
    if (any(task==[31,32,33])) outer_name='iloc'
    call read_string('outer_parameter',outer_name)
    call read_string('inner_parameter',inner_name)
    outer_id=parameter_id(outer_name); inner_id=parameter_id(inner_name)
    if (outer_id==0 .or. inner_id==0 .or. outer_id==inner_id .or. inner_id==5) then
      print *, 'Neutral axes must be distinct: alpha,beta,omega,re,iloc (iloc only outer).'
      return
    end if
    if ((lst_mode==LST_TEMPORAL_MODE .and. (outer_id==3 .or. inner_id==3)) .or. &
        (lst_mode==LST_SPATIAL_MODE .and. (outer_id==1 .or. inner_id==1))) then
      print *, 'Do not scan the dependent eigenvalue (omega temporal / alpha spatial).'
      return
    end if
    if (task==36 .and. outer_id==5) then
      print *, 'Arclength continuation requires continuous axes; iloc is discrete.'
      return
    end if
    center=[parameter_value(outer_id),parameter_value(inner_id)]
    outer_lo=center(1)-.25_real64*max(1._real64,abs(center(1)))
    outer_hi=center(1)+.25_real64*max(1._real64,abs(center(1)))
    inner_lo=center(2)-.25_real64*max(1._real64,abs(center(2)))
    inner_hi=center(2)+.25_real64*max(1._real64,abs(center(2)))
    if (outer_id==5) then
      outer_lo=real(max(1,iloc-1),real64); outer_hi=real(min(in,iloc+1),real64)
    end if
    nrows=3; eigen_tol=1.e-7_real64; jump_tol=.5_real64; continuation_step=.02_real64
    prefix='out/neutral'
    call read_real('outer_min',outer_lo); call read_real('outer_max',outer_hi)
    call read_real('inner_min',inner_lo); call read_real('inner_max',inner_hi)
    call read_int('outer_points',nrows); call read_int('scan_points',opt%scan_samples)
    call read_int('max_iterations',opt%max_iterations); call read_int('max_points',opt%max_points)
    call read_real('growth_tol',opt%f_tol); call read_real('parameter_tol',opt%x_tol)
    call read_real('eigen_error_tol',eigen_tol); call read_real('max_eigen_jump',jump_tol)
    call read_real('continuation_step',continuation_step)
    adaptive_steps=.false.; minimum_step=1.e-5_real64
    vector_tracking=.false.; minimum_overlap=.95_real64
    sample_only=.false.
    reverse_outer=.false.;reverse_inner=.false.;anchor_row_seed=.false.;retry_eigenpairs=0
    if (cfg%has_key('neutral','adaptive_steps')) call cfg%get('neutral','adaptive_steps',adaptive_steps)
    if (cfg%has_key('neutral','vector_tracking')) call cfg%get('neutral','vector_tracking',vector_tracking)
    if (cfg%has_key('neutral','sample_only')) call cfg%get('neutral','sample_only',sample_only)
    if (cfg%has_key('neutral','reverse_outer')) call cfg%get('neutral','reverse_outer',reverse_outer)
    if (cfg%has_key('neutral','reverse_inner')) call cfg%get('neutral','reverse_inner',reverse_inner)
    if (cfg%has_key('neutral','anchor_row_seed')) call cfg%get('neutral','anchor_row_seed',anchor_row_seed)
    call read_int('retry_eigenpairs',retry_eigenpairs)
    call read_real('minimum_step',minimum_step)
    call read_real('minimum_overlap',minimum_overlap)
    if (sample_only .and. task/=33) then
      print *, 'sample_only is available only for task=33 field surveys.'
      return
    end if
    if ((reverse_inner .or. anchor_row_seed) .and. .not.sample_only) then
      print *, 'reverse_inner/anchor_row_seed require a sample_only field survey.'
      return
    end if
    if (retry_eigenpairs<0 .or. (retry_eigenpairs>1 .and. &
        (.not.vector_tracking .or. lst_mode/=LST_SPATIAL_MODE))) then
      print *, 'Candidate retries require spatial PEP and vector_tracking.'
      return
    end if
    call read_real('ds',opt%ds); call read_real('ds_min',opt%ds_min); call read_real('ds_max',opt%ds_max)
    call read_real('finite_difference',opt%finite_difference)
    call read_string('output_prefix',prefix)
    opt%bounds(:,1)=[outer_lo,inner_lo]; opt%bounds(:,2)=[outer_hi,inner_hi]
    opt%use_bisection=(task==32)
    if (.not.all(ieee_is_finite(opt%bounds)) .or. outer_lo>outer_hi .or. inner_lo>=inner_hi .or. &
        nrows<1 .or. opt%scan_samples<2 .or. opt%max_iterations<1 .or. &
        .not.ieee_is_finite(eigen_tol) .or. eigen_tol<=0._real64 .or. &
        .not.ieee_is_finite(jump_tol) .or. jump_tol<=0._real64 .or. &
        .not.ieee_is_finite(continuation_step) .or. continuation_step<=0._real64 .or. &
        .not.ieee_is_finite(minimum_step) .or. minimum_step<=0._real64 .or. &
        .not.ieee_is_finite(minimum_overlap) .or. minimum_overlap<=0._real64 .or. minimum_overlap>1._real64) then
      print *, 'Invalid neutral bounds, counts or tolerances.'
      return
    end if
    if (task==36 .and. outer_lo>=outer_hi) then
      print *, 'Arclength continuation requires a nonzero outer interval.'
      return
    end if
    if (outer_lo==outer_hi .and. nrows/=1) then
      print *, 'Equal outer bounds require outer_points=1.'
      return
    end if
    if (outer_id==5) then
      if (outer_lo<1 .or. outer_hi>in .or. abs(outer_lo-real(nint(outer_lo),real64))>1.e-9_real64 .or. &
          abs(outer_hi-real(nint(outer_hi),real64))>1.e-9_real64) then
        print *, 'iloc bounds must be integer station indices inside [1,in].'
        return
      end if
      if (nrows>1) then
        row_value=(outer_hi-outer_lo)/real(nrows-1,real64)
        if (abs(row_value-real(nint(row_value),real64))>1.e-9_real64 .or. row_value<1._real64) then
          print *, 'iloc spacing must be a positive integer; adjust outer_points.'
          return
        end if
      end if
    end if
    if ((outer_id==4 .and. outer_lo<=0) .or. (inner_id==4 .and. inner_lo<=0)) then
      print *, 'Reynolds numbers must be positive.'
      return
    end if
    ! A single scan row still needs a nonempty geometric bounding box in the core.
    if (outer_lo==outer_hi) then
      opt%bounds(1,1)=outer_lo-1._real64; opt%bounds(1,2)=outer_hi+1._real64
    end if
    call execute_command_line('mkdir -p out')
    have_outputs=.false.; roots_unit=0; field_unit=0; log_unit=0; ridge_unit=0
    open(newunit=log_unit,file=trim(prefix)//'_run.log',status='replace',iostat=ios)
    if (ios/=0) then
      print *, 'Cannot open neutral output prefix: ',trim(prefix)
      return
    end if
    write(log_unit,*) 'task=',task,' mode=',lst_mode,' axes=',trim(outer_name),' / ',trim(inner_name)
    write(log_unit,*) 'bounds=',outer_lo,outer_hi,inner_lo,inner_hi
    write(log_unit,*) 'growth_tol=',opt%f_tol,' eigen_error_tol=',eigen_tol,' max_eigen_jump=',jump_tol
    write(log_unit,*) 'adaptive_steps=',adaptive_steps,' minimum_step=',minimum_step
    write(log_unit,*) 'vector_tracking=',vector_tracking,' minimum_euclidean_overlap=',minimum_overlap
    write(log_unit,*) 'sample_only=',sample_only
    write(log_unit,*) 'reverse_outer=',reverse_outer,' reverse_inner=',reverse_inner,' anchor_row_seed=',anchor_row_seed
    write(log_unit,*) 'retry_eigenpairs=',retry_eigenpairs
    call open_table(roots_unit,trim(prefix)//'_roots.plt',ios)
    if (ios==0) call open_table(field_unit,trim(prefix)//'_field.plt',ios)
    if (ios==0) call open_table(ridge_unit,trim(prefix)//'_ridge.plt',ios)
    if (ios/=0) then
      write(log_unit,*) 'Cannot open result tables.'
      if (roots_unit/=0) close(roots_unit)
      if (field_unit/=0) close(field_unit)
      close(log_unit)
      return
    end if
    have_outputs=.true.; eval_count=0; failed_count=0; solve_count=0; overall=NC_NO_ROOT
    saved_tracking=lst_is_tracking()
    hard_failure=NC_SUCCESS
    allocate(history(10000)); nhistory=1
    allocate(history_vectors(10000),history_has_vector(10000)); history_has_vector=.false.
    history(1)%point=center; history(1)%eigen=seed; history(1)%status=NC_SUCCESS
    if (vector_tracking .and. ieee_is_finite(abs(seed)) .and. lst_eigen_error()<=eigen_tol) &
      call update_reference(history_vectors(1),history_has_vector(1),ok,.false.,seed)
    call lst_track_target(.true.)
    if (task==36) then
      p=center
      call read_real('start_outer',p(1)); call read_real('start_inner',p(2))
      branch=1; n=1
      call read_int('start_branch',branch); call read_int('direction',n)
      best_valid=.false.
      write(field_unit,'(a)') 'ZONE T="row_1"'
      call scan_roots(evaluate,p,2,inner_lo,inner_hi,seed,opt,roots,status)
      if (.not.any(status==[NC_SUCCESS,NC_NO_ROOT])) hard_failure=status
      write(log_unit,*) 'Initial scan status=',neutral_status_name(status)
      if (branch<1 .or. branch>size(roots)) then
        overall=status
        if (size(roots)>0) overall=NC_INVALID_INPUT
        write(log_unit,*) 'Initial neutral branch unavailable; roots=',size(roots),' status=',neutral_status_name(overall)
      else
        has_tangent=cfg%has_key('neutral','tangent_angle')
        if (has_tangent) then
          theta=0._real64; call read_real('tangent_angle',theta)
          theta=theta*acos(-1._real64)/180._real64
          initial_tangent=[cos(theta),sin(theta)]
          call trace_curve(evaluate,roots(branch)%point,roots(branch)%eigen,n,opt,points,status,initial_tangent)
        else
          call trace_curve(evaluate,roots(branch)%point,roots(branch)%eigen,n,opt,points,status)
        end if
        write(roots_unit,'(a)') 'ZONE T="arclength"'
        do n=1,size(points)
          call write_sample(roots_unit,points(n),branch)
        end do
        overall=status
        write(log_unit,*) 'Trace points=',size(points),' termination=',neutral_status_name(status)
      end if
    else
      next_seed=seed
      do row=1,nrows
        row_value=outer_lo
        if (nrows>1) row_value=outer_lo+(outer_hi-outer_lo)*real(row-1,real64)/real(nrows-1,real64)
        if (reverse_outer) row_value=outer_hi-(row_value-outer_lo)
        p=[row_value,center(2)]
        best_valid=.false.;have_row_anchor=.false.
        write(roots_unit,'(a,i0,a)') 'ZONE T="row_',row,'"'
        write(field_unit,'(a,i0,a)') 'ZONE T="row_',row,'"'
        write(ridge_unit,'(a,i0,a)') 'ZONE T="row_',row,'"'
        if (sample_only) then
          if (allocated(roots)) deallocate(roots)
          allocate(roots(0)); status=NC_NO_ROOT; row_seed=next_seed
          do n=1,opt%scan_samples
            p(2)=inner_lo+(inner_hi-inner_lo)*real(n-1,real64)/real(opt%scan_samples-1,real64)
            if (reverse_inner) p(2)=inner_hi-(p(2)-inner_lo)
            call evaluate(p,row_seed,sample_eigen,sample_growth,sample_valid)
            if (sample_valid) then
              row_seed=sample_eigen
              if (n==1) then
                row_anchor_eigen=sample_eigen;have_row_anchor=.true.
              end if
            else
              status=NC_EVALUATION_FAILED
            end if
          end do
        else
          call scan_roots(evaluate,p,2,inner_lo,inner_hi,next_seed,opt,roots,status)
        end if
        if (.not.any(status==[NC_SUCCESS,NC_NO_ROOT])) hard_failure=status
        do branch=1,size(roots)
          call write_sample(roots_unit,roots(branch),branch)
        end do
        if (size(roots)>0) overall=NC_SUCCESS
        if (best_valid) then
          if (any(task==[31,32])) call refine_peak()
          call write_sample(ridge_unit,best,0)
          if (.not.anchor_row_seed) next_seed=best%eigen
        end if
        if (anchor_row_seed .and. have_row_anchor) next_seed=row_anchor_eigen
        write(log_unit,*) 'row=',row,' outer=',row_value,' roots=',size(roots),' status=',neutral_status_name(status)
        flush(roots_unit); flush(field_unit); flush(log_unit)
      end do
    end if
    have_outputs=.false.
    do n=1,size(history_has_vector)
      if (history_has_vector(n)) then
        call VecDestroy(history_vectors(n),ierr)
        call check_vector_error(ierr)
      end if
    end do
    deallocate(history_vectors,history_has_vector)
    ! Restore the interactive local state, including the current eigenvector.
    call lst_track_target(saved_tracking)
    call install_state(saved_iloc,saved_re)
    ome=saved_ome; alf=saved_alf; bat=saved_bat
    call lst%solve(ome,alf,bat)
    ok=ieee_is_finite(real(ome,real64)) .and. ieee_is_finite(aimag(ome)) .and. &
      ieee_is_finite(real(alf,real64)) .and. ieee_is_finite(aimag(alf)) .and. lst_eigen_error()<=eigen_tol
    if (.not.ok) then
      hard_failure=NC_EVALUATION_FAILED
      write(log_unit,*) 'Restoring the original local eigenstate failed.'
    end if
    if (hard_failure/=NC_SUCCESS) overall=hard_failure
    write(log_unit,*) 'evaluations=',eval_count,' solves=',solve_count,' rejected=',failed_count, &
      ' result=',neutral_status_name(overall)
    close(roots_unit); close(field_unit); close(ridge_unit); close(log_unit)
    if (grid%has_assemble) call grid%finalize()
    if (flow%has_assemble) call flow%finalize()
    if (present(result_status)) result_status=overall
    print *, 'Neutral task=',task,' result=',trim(neutral_status_name(overall)),' prefix=',trim(prefix)
  contains
    subroutine read_real(key,value)
      character(*),intent(in) :: key
      real(real64),intent(inout) :: value
      if (cfg%has_key('neutral',key)) call cfg%get('neutral',key,value)
    end subroutine
    subroutine read_int(key,value)
      character(*),intent(in) :: key
      integer,intent(inout) :: value
      if (cfg%has_key('neutral',key)) call cfg%get('neutral',key,value)
    end subroutine
    subroutine read_string(key,value)
      character(*),intent(in) :: key
      character(*),intent(inout) :: value
      if (cfg%has_key('neutral',key)) call cfg%get('neutral',key,value)
    end subroutine
    integer function parameter_id(name) result(id)
      character(*),intent(in) :: name
      character(len=len(name)) :: canonical
      integer :: k,code
      canonical=adjustl(name)
      do k=1,len_trim(canonical)
        code=iachar(canonical(k:k))
        if (code>=65 .and. code<=90) canonical(k:k)=achar(code+32)
      end do
      select case(trim(canonical))
      case('alpha'); id=1
      case('beta'); id=2
      case('omega'); id=3
      case('re'); id=4
      case('iloc'); id=5
      case default; id=0
      end select
    end function
    real(real64) function parameter_value(id) result(value)
      integer,intent(in) :: id
      select case(id)
      case(1); value=real(saved_alf,real64)
      case(2); value=real(saved_bat,real64)
      case(3); value=real(saved_ome,real64)
      case(4); value=saved_re
      case(5); value=real(saved_iloc,real64)
      case default; value=0._real64
      end select
    end function
    subroutine make_state(point,eigen,o,a,b,r,i)
      real(real64),intent(in) :: point(2)
      complex(real64),intent(in) :: eigen
      complex(real64),intent(out) :: o,a,b
      real(real64),intent(out) :: r
      integer,intent(out) :: i
      integer :: axis,id
      o=saved_ome; a=saved_alf; b=saved_bat; r=saved_re; i=saved_iloc
      do axis=1,2
        id=inner_id
        if (axis==1) id=outer_id
        select case(id)
        case(1); a=cmplx(point(axis),0._real64,real64)
        case(2); b=cmplx(point(axis),0._real64,real64)
        case(3); o=cmplx(point(axis),0._real64,real64)
        case(4); r=point(axis)
        case(5); i=nint(point(axis))
        end select
      end do
      if (lst_mode==LST_TEMPORAL_MODE) then
        o=eigen
      else
        a=eigen
      end if
    end subroutine
    subroutine install_state(i,r)
      integer,intent(in) :: i
      real(real64),intent(in) :: r
      logical :: changed
      changed=(iloc/=i .or. Re/=r)
      if (iloc/=i) then
        iloc=i
        grid=assemble_lst_grid_profile(iloc); flow=assemble_lst_baseflow_profile(iloc)
        call lst%set_grid(grid); call lst%set_baseflow(flow)
      end if
      if (Re/=r) then
        Re=r; call lst%set_parameter('Re',Re)
      end if
      if (changed) call lst%setup()
    end subroutine
    subroutine evaluate(point,guess,eigen,g,valid)
      real(real64),intent(in) :: point(2)
      complex(real64),intent(in) :: guess
      complex(real64),intent(out) :: eigen
      real(real64),intent(out) :: g
      logical,intent(out) :: valid
      complex(real64) :: path_seed,previous_eigen,predictor
      real(real64) :: origin(2),trial(2),unit_step(2),distance,score,best_score
      real(real64) :: fraction,first_fraction,last_fraction,previous_fraction,increment,max_increment,trial_increment
      integer :: k,axis,id,nsteps,origin_index,attempts
      type(neutral_sample) :: sample
      Vec :: path_vector
      logical :: has_path_vector
      ! Recover a known point belonging to the requested eigenmode seed, then
      ! approach the requested point in small parameter steps. In particular the
      ! first scan endpoint must not jump directly from the configured seed.
      best_score=huge(1._real64); origin_index=1
      do k=1,nhistory
        score=abs(history(k)%eigen-guess)
        distance=sum(abs(history(k)%point-point)/max(1._real64,abs(center)))
        score=score+1.e-12_real64*distance
        if (score<best_score) then
          best_score=score; origin_index=k
        end if
      end do
      origin=history(origin_index)%point; path_seed=history(origin_index)%eigen
      has_path_vector=.false.
      if (vector_tracking .and. history_has_vector(origin_index)) then
        call VecDuplicate(history_vectors(origin_index),path_vector,ierr)
        call check_vector_error(ierr)
        call VecCopy(history_vectors(origin_index),path_vector,ierr)
        call check_vector_error(ierr)
        has_path_vector=.true.
      end if
      valid=.not.vector_tracking .or. has_path_vector
      eigen=cmplx(ieee_value(0._real64,ieee_quiet_nan),0._real64,real64)
      g=ieee_value(0._real64,ieee_quiet_nan)
      ! Change a discrete station at fixed inner parameter before differentiating
      ! the continuous eigenvalue branch in the inner direction.
      if (outer_id==5 .and. valid) then
        do while (nint(point(1))/=nint(origin(1)))
          trial=[origin(1)+sign(1._real64,point(1)-origin(1)),origin(2)]
          call advance(trial,path_seed,path_seed,eigen,g,valid,path_vector,has_path_vector)
          if (.not.valid) exit
          origin=trial; path_seed=eigen
        end do
      end if
      unit_step=continuation_step
      do axis=1,2
        id=inner_id
        if (axis==1) id=outer_id
        if (id==4) unit_step(axis)=.05_real64*saved_re
        if (id==5) unit_step(axis)=1._real64
      end do
      distance=maxval(abs(point-origin)/unit_step)
      valid=valid .and. ieee_is_finite(distance) .and. distance<=1000._real64
      nsteps=1
      if (valid) nsteps=max(1,ceiling(distance))
      eigen=cmplx(ieee_value(0._real64,ieee_quiet_nan),0._real64,real64)
      g=ieee_value(0._real64,ieee_quiet_nan)
      first_fraction=1._real64
      if (distance>0._real64) first_fraction=min(1._real64,.05_real64/distance)
      last_fraction=0._real64; previous_fraction=0._real64; previous_eigen=path_seed
      if (adaptive_steps .and. valid) then
        ! Retry a failed continuous step from the last accepted point. Neither
        ! residual nor eigenvalue-jump acceptance is relaxed during recovery.
        max_increment=1._real64/real(nsteps,real64)
        increment=first_fraction; attempts=0
        do while (last_fraction<1._real64)
          attempts=attempts+1
          if (attempts>4000) then
            valid=.false.; exit
          end if
          fraction=min(1._real64,last_fraction+increment)
          trial_increment=fraction-last_fraction
          trial=origin+(point-origin)*fraction
          predictor=path_seed
          if (last_fraction>previous_fraction) &
            predictor=path_seed+(path_seed-previous_eigen)*trial_increment/(last_fraction-previous_fraction)
          if (outer_id==5) trial(1)=real(nint(trial(1)),real64)
          call advance(trial,predictor,path_seed,eigen,g,valid,path_vector,has_path_vector)
          if (valid) then
            previous_eigen=path_seed; previous_fraction=last_fraction
            path_seed=eigen; last_fraction=fraction
            increment=min(max_increment,2._real64*trial_increment)
          else
            increment=.5_real64*trial_increment
            if (increment*maxval(abs(point-origin))<max(1.e-12_real64,minimum_step)) exit
            valid=.true.
          end if
        end do
      else
      do k=1,nsteps+1
        if (.not.valid) exit
        fraction=first_fraction
        if (k>1) fraction=first_fraction+(1._real64-first_fraction)*real(k-1,real64)/real(nsteps,real64)
        trial=origin+(point-origin)*fraction
        predictor=path_seed
        if (k>1 .and. last_fraction>previous_fraction) &
          predictor=path_seed+(path_seed-previous_eigen)*(fraction-last_fraction)/(last_fraction-previous_fraction)
        if (outer_id==5) trial(1)=real(nint(trial(1)),real64)
        call advance(trial,predictor,path_seed,eigen,g,valid,path_vector,has_path_vector)
        if (valid) then
          previous_eigen=path_seed; previous_fraction=last_fraction
          path_seed=eigen; last_fraction=fraction
        end if
        if (fraction>=1._real64) exit
      end do
      end if
      eval_count=eval_count+1
      sample%point=point; sample%eigen=eigen; sample%growth=g; sample%residual=abs(g)
      sample%status=NC_EVALUATION_FAILED
      if (valid) then
        sample%status=NC_SUCCESS
        if (nhistory<size(history)) then
          nhistory=nhistory+1; history(nhistory)=sample
        else
          history(nhistory)=sample
        end if
        if (vector_tracking .and. has_path_vector) then
          if (history_has_vector(nhistory)) then
            call VecDestroy(history_vectors(nhistory),ierr)
            call check_vector_error(ierr)
          end if
          call VecDuplicate(path_vector,history_vectors(nhistory),ierr)
          call check_vector_error(ierr)
          call VecCopy(path_vector,history_vectors(nhistory),ierr)
          call check_vector_error(ierr)
          history_has_vector(nhistory)=.true.
        end if
        if (.not.best_valid) then
          best=sample; best_valid=.true.
        else if (g>best%growth) then
          best=sample
        end if
      else
        failed_count=failed_count+1
      end if
      if (have_outputs) call write_sample(field_unit,sample,0)
      if (has_path_vector) then
        call VecDestroy(path_vector,ierr)
        call check_vector_error(ierr)
      end if
    end subroutine
    subroutine check_vector_error(error)
      PetscErrorCode, intent(in) :: error
      PetscErrorCode :: report_error
      if (error/=0) then
        report_error=error
        call PetscErrorF(report_error,__LINE__,"module_neutral_curve.f90")
        call MPIU_Abort(PETSC_COMM_WORLD,error)
      end if
    end subroutine
    subroutine update_reference(reference,has_reference,valid,enforce,expected,reference_eigen)
      Vec, intent(inout) :: reference
      logical, intent(inout) :: has_reference
      logical, intent(out) :: valid
      logical, intent(in) :: enforce
      complex(real64), intent(inout) :: expected
      complex(real64), intent(in), optional :: reference_eigen
      Vec :: current
      PetscScalar :: product,selected_eigen
      real(real64) :: selected_error
      PetscReal :: norm_reference,norm_current,overlap
      if (enforce .and. has_reference .and. present(reference_eigen)) then
        call lst_select_eigenvector_by_overlap(reference,reference_eigen,jump_tol,eigen_tol,minimum_overlap, &
          selected_eigen,current,selected_error,overlap,valid)
        if (valid) expected=selected_eigen
      else
        call lst_copy_selected_eigenvector(selected_eigen,current,selected_error,valid)
      end if
      if (.not.valid) return
      valid=selected_error<=eigen_tol .and. abs(selected_eigen-expected)<=1.e-10_real64*max(1._real64,abs(expected))
      if (.not.valid) then
        call VecDestroy(current,ierr)
        call check_vector_error(ierr)
        return
      end if
      call VecNorm(current,NORM_2,norm_current,ierr)
      call check_vector_error(ierr)
      valid=ieee_is_finite(norm_current) .and. norm_current>0._real64
      overlap=1._real64
      if (valid .and. has_reference .and. enforce) then
        call VecNorm(reference,NORM_2,norm_reference,ierr)
        call check_vector_error(ierr)
        call VecDot(reference,current,product,ierr)
        call check_vector_error(ierr)
        valid=ieee_is_finite(norm_reference) .and. norm_reference>0._real64
        if (valid) then
          overlap=abs(product/(norm_reference*norm_current))**2
          valid=ieee_is_finite(overlap) .and. overlap>=minimum_overlap
        end if
        write(log_unit,*) 'vector_overlap=',overlap,' accepted=',valid
      end if
      if (valid) then
        if (.not.has_reference) then
          call VecDuplicate(current,reference,ierr)
          call check_vector_error(ierr)
        end if
        call VecCopy(current,reference,ierr)
        call check_vector_error(ierr)
        has_reference=.true.
      end if
      call VecDestroy(current,ierr)
      call check_vector_error(ierr)
    end subroutine
    subroutine advance(point,guess,reference_eigen,eigen,g,valid,reference,has_reference)
      real(real64), intent(in) :: point(2)
      complex(real64), intent(in) :: guess,reference_eigen
      complex(real64), intent(out) :: eigen
      real(real64), intent(out) :: g
      logical, intent(out) :: valid
      Vec, intent(inout) :: reference
      logical, intent(inout) :: has_reference
      logical :: performed
      call solve_trial(point,guess,eigen,g,valid,performed)
      if (vector_tracking .and. performed) then
        ! The nearest candidate may fail while another converged pair matches
        ! the reference. Inspect the full returned candidate set before rejection.
        call update_reference(reference,has_reference,valid,.true.,eigen,reference_eigen)
        if (.not.valid .and. has_reference .and. retry_eigenpairs>1) then
          write(log_unit,*) 'candidate_retry=',retry_eigenpairs,' point=',point
          call lst_prepare_candidate_retry(retry_eigenpairs,reference)
          call solve_trial(point,guess,eigen,g,valid,performed)
          if (performed) call update_reference(reference,has_reference,valid,.true.,eigen,reference_eigen)
        end if
        if (valid) then
          if (lst_mode==LST_TEMPORAL_MODE) then
            g=aimag(eigen)
          else
            g=-aimag(eigen)
          end if
        end if
      else if (valid) then
        valid=abs(eigen-reference_eigen)<=jump_tol
      end if
    end subroutine
    subroutine solve_trial(point,guess,eigen,g,valid,performed)
      real(real64), intent(in) :: point(2)
      complex(real64), intent(in) :: guess
      complex(real64), intent(out) :: eigen
      real(real64), intent(out) :: g
      logical, intent(out) :: valid
      logical, intent(out) :: performed
      complex(real64) :: o,a,b
      real(real64) :: r
      integer :: i
      performed=.false.
      eigen=cmplx(ieee_value(0._real64,ieee_quiet_nan),0._real64,real64)
      g=ieee_value(0._real64,ieee_quiet_nan)
      call make_state(point,guess,o,a,b,r,i)
      valid=i>=1 .and. i<=in .and. r>0._real64
      if (.not.valid) return
      call install_state(i,r)
      call lst%solve(o,a,b)
      performed=.true.
      solve_count=solve_count+1
      if (lst_mode==LST_TEMPORAL_MODE) then
        eigen=o; g=aimag(o)
      else
        eigen=a; g=-aimag(a)
      end if
      valid=ieee_is_finite(real(eigen,real64)) .and. ieee_is_finite(aimag(eigen)) .and. &
        ieee_is_finite(g) .and. lst_eigen_error()<=eigen_tol
    end subroutine
    subroutine open_table(unit,filename,code)
      integer,intent(out) :: unit,code
      character(*),intent(in) :: filename
      unit=0
      open(newunit=unit,file=filename,status='replace',iostat=code)
      if (code==0) write(unit,'(a)') 'VARIABLES="outer","inner","iloc","x","Re",'// &
        '"alpha_r","alpha_i","beta_r","beta_i","omega_r","omega_i","growth","branch","status"'
    end subroutine
    subroutine write_sample(unit,sample,branch_index)
      integer,intent(in) :: unit,branch_index
      type(neutral_sample),intent(in) :: sample
      complex(real64) :: o,a,b
      real(real64) :: r
      integer :: i
      call make_state(sample%point,sample%eigen,o,a,b,r,i)
      write(unit,'(2(es24.15,1x),i0,1x,9(es24.15,1x),2(i0,1x))') sample%point,i, &
        xx(i,yb+1),r,a,b,o,sample%growth,branch_index,sample%status
    end subroutine
    subroutine refine_peak()
      real(real64) :: left,right,c,d,gc,gd,step
      complex(real64) :: ec,ed,peak_seed
      logical :: vc,vd
      integer :: iteration
      step=(inner_hi-inner_lo)/real(opt%scan_samples-1,real64)
      left=max(inner_lo,best%point(2)-step); right=min(inner_hi,best%point(2)+step)
      peak_seed=best%eigen
      do iteration=1,min(24,opt%max_iterations)
        if (right-left<=opt%x_tol) exit
        c=right-(right-left)*.6180339887498949_real64
        d=left+(right-left)*.6180339887498949_real64
        call evaluate([row_value,c],peak_seed,ec,gc,vc)
        call evaluate([row_value,d],peak_seed,ed,gd,vd)
        if (.not.vc .or. .not.vd) then
          hard_failure=NC_EVALUATION_FAILED
          write(log_unit,*) 'Peak refinement failed; retaining best valid sampled point.'
          exit
        end if
        if (gc>gd) then
          right=d; peak_seed=ec
        else
          left=c; peak_seed=ed
        end if
      end do
    end subroutine
  end subroutine run_neutral_curve
end module mod_neutral_curve
