! Standalone bounded neutral-root and pseudo-arclength algorithms.
! point=(outer coordinate,inner coordinate), bounds(:,1)=lower bounds,
! bounds(:,2)=upper bounds. evaluate must follow one continuous eigenmode;
! a successful callback must return finite eigen and growth values.
module mod_neutral_numerics
  use, intrinsic :: iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private
  integer, parameter, public :: dp=real64
  integer, parameter, public :: NC_SUCCESS=0, NC_NO_ROOT=1, NC_EVALUATION_FAILED=2, &
    NC_MAX_ITERATIONS=3, NC_BOUNDARY=4, NC_CLOSED=5, NC_STEP_TOO_SMALL=6, &
    NC_INVALID_INPUT=7, NC_DEGENERATE=8

  type, public :: neutral_options
    integer :: scan_samples=81, max_iterations=60, max_points=400, min_closure_points=12
    real(dp) :: f_tol=1.e-8_dp, x_tol=1.e-9_dp
    real(dp) :: ds=.01_dp, ds_min=1.e-5_dp, ds_max=.1_dp
    real(dp) :: finite_difference=1.e-5_dp
    logical :: use_bisection=.false.
    real(dp) :: bounds(2,2)=reshape([-huge(1._dp),-huge(1._dp),huge(1._dp),huge(1._dp)],[2,2])
  end type

  type, public :: neutral_sample
    real(dp) :: point(2)=0._dp
    complex(dp) :: eigen=(0._dp,0._dp)
    real(dp) :: growth=huge(1._dp), residual=huge(1._dp)
    integer :: status=NC_INVALID_INPUT, iterations=0
  end type

  abstract interface
    subroutine neutral_evaluator(point,seed,eigen,growth,ok)
      import dp
      real(dp), intent(in) :: point(2)
      complex(dp), intent(in) :: seed
      complex(dp), intent(out) :: eigen
      real(dp), intent(out) :: growth
      logical, intent(out) :: ok
    end subroutine
  end interface
  public :: neutral_evaluator, scan_roots, trace_curve, neutral_status_name

contains

  function neutral_status_name(status) result(name)
    integer, intent(in) :: status
    character(len=24) :: name
    select case(status)
    case(NC_SUCCESS); name='converged'
    case(NC_NO_ROOT); name='no bracketed root'
    case(NC_EVALUATION_FAILED); name='evaluation failed'
    case(NC_MAX_ITERATIONS); name='iteration/point limit'
    case(NC_BOUNDARY); name='parameter boundary'
    case(NC_CLOSED); name='closed curve'
    case(NC_STEP_TOO_SMALL); name='step too small'
    case(NC_INVALID_INPUT); name='invalid input'
    case(NC_DEGENERATE); name='zero/singular gradient'
    case default; name='unknown'
    end select
  end function

  logical function valid_options(options)
    type(neutral_options), intent(in) :: options
    valid_options=all(ieee_is_finite(options%bounds)) .and. &
      all(options%bounds(:,1)<options%bounds(:,2)) .and. &
      ieee_is_finite(options%f_tol) .and. options%f_tol>0._dp .and. &
      ieee_is_finite(options%x_tol) .and. options%x_tol>0._dp .and. &
      ieee_is_finite(options%finite_difference) .and. options%finite_difference>0._dp .and. &
      options%max_iterations>0
  end function

  logical function inside(point,options)
    real(dp), intent(in) :: point(2)
    type(neutral_options), intent(in) :: options
    inside=all(ieee_is_finite(point)) .and. all(point>=options%bounds(:,1)) .and. &
      all(point<=options%bounds(:,2))
  end function

  subroutine evaluate_checked(evaluate,point,seed,options,sample)
    procedure(neutral_evaluator) :: evaluate
    real(dp), intent(in) :: point(2)
    complex(dp), intent(in) :: seed
    type(neutral_options), intent(in) :: options
    type(neutral_sample), intent(out) :: sample
    logical :: ok
    sample=neutral_sample()
    sample%point=point
    if (.not.inside(point,options)) then
      sample%status=NC_BOUNDARY
      return
    end if
    if (.not.ieee_is_finite(real(seed,dp)) .or. .not.ieee_is_finite(aimag(seed))) then
      sample%status=NC_INVALID_INPUT
      return
    end if
    call evaluate(point,seed,sample%eigen,sample%growth,ok)
    if (.not.ok) then
      sample%status=NC_EVALUATION_FAILED
      return
    end if
    if (.not.ieee_is_finite(sample%growth) .or. &
        .not.ieee_is_finite(real(sample%eigen,dp)) .or. .not.ieee_is_finite(aimag(sample%eigen))) then
      sample%status=NC_EVALUATION_FAILED
      return
    end if
    sample%residual=abs(sample%growth)
    sample%status=NC_SUCCESS
  end subroutine

  ! Finds every sampled sign-change root and every sampled point whose residual
  ! meets f_tol. A finite scan cannot guarantee detection of tangential roots or
  ! multiple roots entirely between adjacent scan nodes. Failed evaluations cut
  ! the scan into separate intervals: no root bracket crosses such a point.
  ! status reports any failure even when roots(:) contains valid partial results.
  subroutine scan_roots(evaluate,fixed_point,axis,lower,upper,seed,options,roots,status)
    procedure(neutral_evaluator) :: evaluate
    real(dp), intent(in) :: fixed_point(2),lower,upper
    integer, intent(in) :: axis
    complex(dp), intent(in) :: seed
    type(neutral_options), intent(in) :: options
    type(neutral_sample), allocatable, intent(out) :: roots(:)
    integer, intent(out) :: status
    type(neutral_sample), allocatable :: work(:)
    type(neutral_sample) :: previous,current,root
    real(dp) :: point(2)
    complex(dp) :: current_seed
    integer :: i,count_roots,failure
    logical :: have_previous
    allocate(roots(0))
    status=NC_INVALID_INPUT
    if (.not.valid_options(options) .or. options%scan_samples<2) return
    if (axis<1 .or. axis>2) return
    if (.not.ieee_is_finite(lower) .or. .not.ieee_is_finite(upper)) return
    if (lower>=upper .or. lower<options%bounds(axis,1) .or. upper>options%bounds(axis,2)) return
    point=fixed_point
    point(axis)=lower
    if (.not.inside(point,options)) return
    allocate(work(2*options%scan_samples))
    count_roots=0
    failure=NC_SUCCESS
    have_previous=.false.
    current_seed=seed
    do i=1,options%scan_samples
      point=fixed_point
      point(axis)=lower+(upper-lower)*real(i-1,dp)/real(options%scan_samples-1,dp)
      call evaluate_checked(evaluate,point,current_seed,options,current)
      if (current%status/=NC_SUCCESS) then
        failure=current%status
        have_previous=.false.
        cycle
      end if
      current_seed=current%eigen
      if (current%residual<=options%f_tol) then
        call append_unique(current)
      else if (have_previous) then
        if (previous%residual>options%f_tol .and. &
            ((previous%growth<0._dp) .neqv. (current%growth<0._dp))) then
          call correct_bracket(evaluate,previous,current,axis,options,root)
          if (root%status==NC_SUCCESS) then
            call append_unique(root)
          else
            failure=root%status
          end if
        end if
      end if
      previous=current
      have_previous=.true.
    end do
    deallocate(roots)
    allocate(roots(count_roots))
    roots=work(:count_roots)
    status=NC_SUCCESS
    if (count_roots==0) status=NC_NO_ROOT
    if (failure/=NC_SUCCESS) status=failure
  contains
    subroutine append_unique(value)
      type(neutral_sample), intent(in) :: value
      if (count_roots>0) then
        if (abs(value%point(axis)-work(count_roots)%point(axis))<=options%x_tol) return
      end if
      count_roots=count_roots+1
      work(count_roots)=value
    end subroutine
  end subroutine

  subroutine correct_bracket(evaluate,left_input,right_input,axis,options,root)
    procedure(neutral_evaluator) :: evaluate
    type(neutral_sample), intent(in) :: left_input,right_input
    integer, intent(in) :: axis
    type(neutral_options), intent(in) :: options
    type(neutral_sample), intent(out) :: root
    type(neutral_sample) :: left,right
    real(dp) :: x,lo,hi,width,ratio,point(2)
    complex(dp) :: seed
    integer :: iteration
    left=left_input
    right=right_input
    root=left
    do iteration=1,options%max_iterations
      lo=left%point(axis)
      hi=right%point(axis)
      width=hi-lo
      ! Opposite signs make this ratio well-conditioned without multiplying f.
      ratio=abs(left%growth)/max(abs(left%growth),abs(right%growth))
      ratio=ratio/(ratio+abs(right%growth)/max(abs(left%growth),abs(right%growth)))
      x=lo+width*ratio
      if (x<=lo+.1_dp*width .or. x>=hi-.1_dp*width) x=lo+.5_dp*width
      if (options%use_bisection) x=lo+.5_dp*width
      if (x<=lo .or. x>=hi) exit
      point=left%point
      point(axis)=x
      seed=left%eigen
      if (x-lo>hi-x) seed=right%eigen
      call evaluate_checked(evaluate,point,seed,options,root)
      root%iterations=iteration
      if (root%status/=NC_SUCCESS) return
      if (root%residual<=options%f_tol) return
      if ((root%growth<0._dp) .eqv. (left%growth<0._dp)) then
        left=root
      else
        right=root
      end if
      ! Interval width alone never certifies a neutral root.
    end do
    root%status=NC_MAX_ITERATIONS
  end subroutine

  subroutine estimate_gradient(evaluate,center,options,gradient,status)
    procedure(neutral_evaluator) :: evaluate
    type(neutral_sample), intent(in) :: center
    type(neutral_options), intent(in) :: options
    real(dp), intent(out) :: gradient(2)
    integer, intent(out) :: status
    type(neutral_sample) :: plus,minus
    real(dp) :: point(2),h
    integer :: axis
    gradient=0._dp
    status=NC_SUCCESS
    do axis=1,2
      h=options%finite_difference*max(1._dp,abs(center%point(axis)))
      point=center%point
      point(axis)=center%point(axis)+h
      call evaluate_checked(evaluate,point,center%eigen,options,plus)
      point(axis)=center%point(axis)-h
      call evaluate_checked(evaluate,point,center%eigen,options,minus)
      if (plus%status==NC_SUCCESS .and. minus%status==NC_SUCCESS) then
        gradient(axis)=(plus%growth-minus%growth)/(2._dp*h)
      else if (plus%status==NC_SUCCESS .and. minus%status==NC_BOUNDARY) then
        gradient(axis)=(plus%growth-center%growth)/h
      else if (minus%status==NC_SUCCESS .and. plus%status==NC_BOUNDARY) then
        gradient(axis)=(center%growth-minus%growth)/h
      else
        status=NC_EVALUATION_FAILED
        if (plus%status==NC_BOUNDARY .and. minus%status==NC_BOUNDARY) status=NC_BOUNDARY
        return
      end if
    end do
    if (.not.all(ieee_is_finite(gradient))) status=NC_EVALUATION_FAILED
    if (norm2(gradient)<=tiny(1._dp)) status=NC_DEGENERATE
  end subroutine

  ! Initial Newton projection may adjust both coordinates. Later correction is
  ! confined to the normal through the predictor, i.e. the pseudo-arclength
  ! constraint dot(point-predictor,tangent)=0 is preserved exactly.
  subroutine project_neutral(evaluate,initial,seed,options,result,normal,max_distance)
    procedure(neutral_evaluator) :: evaluate
    real(dp), intent(in) :: initial(2)
    complex(dp), intent(in) :: seed
    type(neutral_options), intent(in) :: options
    type(neutral_sample), intent(out) :: result
    real(dp), optional, intent(in) :: normal(2),max_distance
    type(neutral_sample) :: trial
    real(dp) :: gradient(2),delta(2),direction(2),derivative,distance,fraction
    integer :: iteration,status,backtrack,last_failure
    call evaluate_checked(evaluate,initial,seed,options,result)
    if (result%status/=NC_SUCCESS) return
    do iteration=1,options%max_iterations
      result%iterations=iteration-1
      if (result%residual<=options%f_tol) return
      call estimate_gradient(evaluate,result,options,gradient,status)
      if (status/=NC_SUCCESS) then
        result%status=status
        return
      end if
      direction=gradient/norm2(gradient)
      if (present(normal)) direction=normal
      derivative=dot_product(gradient,direction)
      if (abs(derivative)<=sqrt(epsilon(1._dp))*norm2(gradient)) then
        result%status=NC_DEGENERATE
        return
      end if
      delta=(-result%growth/derivative)*direction
      distance=norm2(delta)
      if (distance>options%ds_max) delta=delta*(options%ds_max/distance)
      fraction=1._dp
      last_failure=NC_DEGENERATE
      do backtrack=1,12
        if (present(max_distance)) then
          if (norm2(result%point+fraction*delta-initial)>max_distance) then
            fraction=.5_dp*fraction
            cycle
          end if
        end if
        call evaluate_checked(evaluate,result%point+fraction*delta,result%eigen,options,trial)
        if (trial%status==NC_SUCCESS) then
          if (trial%residual<result%residual .or. trial%residual<=options%f_tol) exit
        else
          last_failure=trial%status
        end if
        fraction=.5_dp*fraction
      end do
      if (backtrack>12) then
        result%status=last_failure
        return
      end if
      result=trial
      result%iterations=iteration
      if (result%residual<=options%f_tol) return
      ! Small steps with nonzero residual are failures, never false convergence.
      if (norm2(fraction*delta)<=options%x_tol .and. result%residual>options%f_tol) then
        result%status=NC_MAX_ITERATIONS
        return
      end if
    end do
    result%status=NC_MAX_ITERATIONS
    result%iterations=options%max_iterations
  end subroutine

  subroutine trace_curve(evaluate,start_point,seed,direction,options,points,status,initial_tangent)
    procedure(neutral_evaluator) :: evaluate
    real(dp), intent(in) :: start_point(2)
    complex(dp), intent(in) :: seed
    integer, intent(in) :: direction
    type(neutral_options), intent(in) :: options
    type(neutral_sample), allocatable, intent(out) :: points(:)
    integer, intent(out) :: status
    real(dp), optional, intent(in) :: initial_tangent(2)
    type(neutral_sample), allocatable :: work(:)
    type(neutral_sample) :: current,candidate
    real(dp) :: gradient(2),normal(2),tangent(2),next_tangent(2),first_tangent(2),predictor(2)
    real(dp) :: step,travel,distance
    integer :: count_points,gradient_status,last_failure
    allocate(points(0))
    status=NC_INVALID_INPUT
    if (.not.valid_options(options)) return
    if (abs(direction)/=1 .or. options%max_points<2 .or. options%min_closure_points<4) return
    if (.not.ieee_is_finite(options%ds) .or. .not.ieee_is_finite(options%ds_min) .or. &
        .not.ieee_is_finite(options%ds_max)) return
    if (options%ds_min<=0._dp .or. options%ds<options%ds_min .or. options%ds>options%ds_max) return
    allocate(work(options%max_points))
    call project_neutral(evaluate,start_point,seed,options,current)
    status=current%status
    if (status/=NC_SUCCESS) return
    count_points=1
    work(1)=current
    call estimate_gradient(evaluate,current,options,gradient,gradient_status)
    if (gradient_status/=NC_SUCCESS) then
      status=gradient_status
      call finish()
      return
    end if
    normal=gradient/norm2(gradient)
    tangent=[-normal(2),normal(1)]
    if (present(initial_tangent)) then
      if (.not.all(ieee_is_finite(initial_tangent)) .or. norm2(initial_tangent)<=tiny(1._dp)) then
        status=NC_INVALID_INPUT
        call finish()
        return
      end if
      if (abs(dot_product(tangent,initial_tangent))<=sqrt(epsilon(1._dp))*norm2(initial_tangent)) then
        status=NC_DEGENERATE
        call finish()
        return
      end if
      if (dot_product(tangent,initial_tangent)<0._dp) tangent=-tangent
    end if
    tangent=real(direction,dp)*tangent
    first_tangent=tangent
    travel=0._dp
    step=options%ds
    status=NC_MAX_ITERATIONS
    do while(count_points<options%max_points)
      last_failure=NC_STEP_TOO_SMALL
      do
        predictor=current%point+step*tangent
        if (.not.inside(predictor,options)) then
          last_failure=NC_BOUNDARY
        else
          call project_neutral(evaluate,predictor,current%eigen,options,candidate,normal,2._dp*step)
          last_failure=candidate%status
          if (candidate%status==NC_SUCCESS) exit
        end if
        if (step*.5_dp<options%ds_min) exit
        step=step*.5_dp
      end do
      if (last_failure/=NC_SUCCESS) then
        status=last_failure
        exit
      end if
      distance=norm2(candidate%point-current%point)
      if (distance<=options%x_tol) then
        status=NC_STEP_TOO_SMALL
        exit
      end if
      call estimate_gradient(evaluate,candidate,options,gradient,gradient_status)
      if (gradient_status/=NC_SUCCESS) then
        status=gradient_status
        exit
      end if
      normal=gradient/norm2(gradient)
      next_tangent=[-normal(2),normal(1)]
      if (dot_product(next_tangent,tangent)<0._dp) next_tangent=-next_tangent
      travel=travel+distance
      if (count_points>=options%min_closure_points .and. travel>6._dp*options%ds) then
        if (norm2(candidate%point-work(1)%point)<=.8_dp*step .and. &
            dot_product(next_tangent,first_tangent)>.5_dp) then
          count_points=count_points+1
          work(count_points)=work(1)
          status=NC_CLOSED
          exit
        end if
      end if
      count_points=count_points+1
      work(count_points)=candidate
      current=candidate
      tangent=next_tangent
      step=min(options%ds,1.25_dp*step)
    end do
    call finish()
  contains
    subroutine finish()
      deallocate(points)
      allocate(points(count_points))
      points=work(:count_points)
    end subroutine
  end subroutine
end module mod_neutral_numerics
