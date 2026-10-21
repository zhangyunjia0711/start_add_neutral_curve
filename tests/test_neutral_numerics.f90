program test_neutral_numerics
  use mod_neutral_numerics
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  implicit none
  type(neutral_options) :: options
  type(neutral_sample), allocatable :: roots(:),points(:),reverse_points(:)
  integer :: mode,status,checks
  real(dp) :: expected,orientation
  complex(dp), parameter :: seed=(3._dp,0._dp)
  checks=0
  options%bounds(:,1)=[-3._dp,-2._dp]
  options%bounds(:,2)=[3._dp,2._dp]
  options%scan_samples=53
  options%f_tol=1.e-10_dp
  options%x_tol=1.e-12_dp
  mode=1
  call scan_roots(analytical,[.5_dp,0._dp],2,-1.5_dp,1.5_dp,seed,options,roots,status)
  call check(status==NC_SUCCESS,'ellipse scan status')
  call check(size(roots)==2,'ellipse has two neutral branches')
  expected=sqrt(1._dp-.5_dp**2/4._dp)
  call check(abs(roots(1)%point(2)+expected)<1.e-8_dp,'lower ellipse root')
  call check(abs(roots(2)%point(2)-expected)<1.e-8_dp,'upper ellipse root')
  call check(all(roots%residual<=options%f_tol),'root residuals')
  options%use_bisection=.true.
  call scan_roots(analytical,[.5_dp,0._dp],2,-1.5_dp,1.5_dp,seed,options,roots,status)
  call check(status==NC_SUCCESS .and. size(roots)==2,'pure bisection finds both branches')
  call check(all(roots%residual<=options%f_tol),'bisection residuals')
  options%use_bisection=.false.

  options%ds=.08_dp
  options%ds_max=.2_dp
  options%max_points=400
  call trace_curve(analytical,[2.1_dp,.05_dp],seed,1,options,points,status)
  call check(status==NC_CLOSED,'ellipse closes from initially nonneutral point')
  call check(size(points)>50,'ellipse has meaningful number of continuation points')
  call check(all(points%residual<=options%f_tol),'all trace residuals certified')
  call check(norm2(points(1)%point-points(size(points))%point)<1.e-14_dp,'closed endpoints coincide')
  call check(min_coordinate(points,1)<-1.95_dp .and. max_coordinate(points,1)>1.95_dp,'both x turning points')
  call check(min_coordinate(points,2)<-.98_dp .and. max_coordinate(points,2)>.98_dp,'both y turning points')
  orientation=points(2)%point(2)-points(1)%point(2)
  call trace_curve(analytical,[2.1_dp,.05_dp],seed,-1,options,reverse_points,status)
  call check(status==NC_CLOSED,'reverse trace closes')
  call check(orientation*(reverse_points(2)%point(2)-reverse_points(1)%point(2))<0._dp,'direction reverses')
  call trace_curve(analytical,[2.1_dp,.05_dp],seed,1,options,reverse_points,status,[0._dp,1._dp])
  call check(status==NC_CLOSED,'initial tangent supported')
  call check(reverse_points(2)%point(2)>reverse_points(1)%point(2),'initial tangent orientation')

  mode=2
  call scan_roots(analytical,[0._dp,0._dp],2,-1._dp,1._dp,seed,options,roots,status)
  call check(status==NC_EVALUATION_FAILED .and. size(roots)==0,'NaN gap is not bridged')
  mode=3
  call scan_roots(analytical,[0._dp,0._dp],2,-1._dp,1._dp,seed,options,roots,status)
  call check(status==NC_EVALUATION_FAILED .and. size(roots)==0,'failed gap is not bridged')
  mode=4
  call scan_roots(analytical,[0._dp,0._dp],2,-1._dp,1._dp,seed,options,roots,status)
  call check(status==NC_NO_ROOT .and. size(roots)==0,'unsampled tangency is not falsely claimed')
  mode=5
  call scan_roots(analytical,[0._dp,0._dp],2,-1.5_dp,1.5_dp,seed,options,roots,status)
  call check(status==NC_SUCCESS .and. size(roots)==4,'all four disjoint roots returned')
  mode=6
  call scan_roots(analytical,[0._dp,0._dp],2,-1._dp,1._dp,seed,options,roots,status)
  call check(status==NC_MAX_ITERATIONS .and. size(roots)==0,'discontinuity never passes on bracket width')
  mode=7
  call scan_roots(analytical,[0._dp,0._dp],2,-1._dp,1._dp,seed,options,roots,status)
  call check(status==NC_NO_ROOT .and. size(roots)==0,'constant growth has no roots')
  call trace_curve(analytical,[0._dp,0._dp],seed,1,options,points,status)
  call check(status==NC_DEGENERATE .and. size(points)==0,'flat initial gradient fails explicitly')
  mode=9
  call scan_roots(analytical,[0._dp,0._dp],2,-1._dp,1._dp,seed,options,roots,status)
  call check(status==NC_EVALUATION_FAILED .and. size(roots)==0,'nonfinite eigenvalue rejected')

  mode=1
  options%max_points=5
  call trace_curve(analytical,[2._dp,0._dp],seed,1,options,points,status)
  call check(status==NC_MAX_ITERATIONS .and. size(points)==5,'trace point limit')
  options%max_points=400
  options%bounds(:,1)=[1.5_dp,-.2_dp]
  options%bounds(:,2)=[2.5_dp,.2_dp]
  call trace_curve(analytical,[2._dp,0._dp],seed,1,options,points,status)
  call check(status==NC_BOUNDARY,'trace parameter boundary')
  call check(all(points%residual<=options%f_tol),'bounded partial trace retains valid points')
  options%bounds(:,1)=[-3._dp,-2._dp]
  options%bounds(:,2)=[3._dp,2._dp]
  mode=8
  call trace_curve(analytical,[2._dp,0._dp],seed,-1,options,points,status)
  call check(status==NC_EVALUATION_FAILED,'trace callback failure is explicit')
  call check(all(points%residual<=options%f_tol),'failed partial trace retains valid points')
  options%scan_samples=1
  call scan_roots(analytical,[0._dp,0._dp],2,-1._dp,1._dp,seed,options,roots,status)
  call check(status==NC_INVALID_INPUT,'invalid scan options')
  print '(a,i0,a)','PASS: ',checks,' analytical neutral-numerics checks'
contains
  subroutine check(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not.condition) then
      print '(a,a)','FAIL: ',label
      error stop 1
    end if
    checks=checks+1
  end subroutine

  real(dp) function min_coordinate(values,axis)
    type(neutral_sample), intent(in) :: values(:)
    integer, intent(in) :: axis
    integer :: i
    min_coordinate=huge(1._dp)
    do i=1,size(values)
      min_coordinate=min(min_coordinate,values(i)%point(axis))
    end do
  end function

  real(dp) function max_coordinate(values,axis)
    type(neutral_sample), intent(in) :: values(:)
    integer, intent(in) :: axis
    integer :: i
    max_coordinate=-huge(1._dp)
    do i=1,size(values)
      max_coordinate=max(max_coordinate,values(i)%point(axis))
    end do
  end function

  subroutine analytical(point,initial_seed,eigen,growth,ok)
    real(dp), intent(in) :: point(2)
    complex(dp), intent(in) :: initial_seed
    complex(dp), intent(out) :: eigen
    real(dp), intent(out) :: growth
    logical, intent(out) :: ok
    ok=.true.
    growth=1._dp-point(1)**2/4._dp-point(2)**2
    select case(mode)
    case(2,3)
      growth=point(2)
      if (abs(point(2))<.2_dp) then
        if (mode==2) growth=ieee_value(1._dp,ieee_quiet_nan)
        if (mode==3) ok=.false.
      end if
    case(4)
      growth=(point(2)-.331_dp)**2
    case(5)
      growth=(point(2)**2-.25_dp)*(point(2)**2-1._dp)
    case(6)
      growth=1._dp
      if (point(2)<.231_dp) growth=-1._dp
    case(7)
      growth=1._dp
    case(8)
      if (point(2)>.25_dp) ok=.false.
    end select
    eigen=cmplx(3._dp+.1_dp*point(1),growth,dp)
    if (mode==9) eigen=cmplx(ieee_value(1._dp,ieee_quiet_nan),growth,dp)
    if (abs(initial_seed)>1.e6_dp) ok=.false.
  end subroutine
end program
