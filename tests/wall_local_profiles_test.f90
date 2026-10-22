program test_wall_local_profiles
   use iso_fortran_env, only: real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use mod_wall_local_coordinates
   use global_parameters, only: xx,yy,qq,drx,drz,in,jn,yb,ye,rank,world_size=>size, &
                                eta_y,eta_yy,dq0dy,dq0dyy
   implicit none
   real(real64) :: wx(7),wy(7),s(7),t(2),normal(2),hint(2),theta,nvalue
   real(real64) :: n(16),q(5,16),qn(5,16),qnn(5,16),en(16),enn(16),expected(5,16)
   real(real64) :: reference_en(16),reference_enn(16),reference_qn(5,16),reference_qnn(5,16)
   integer :: i,j,unit
   character(32) :: argument

   call get_command_argument(1,argument)
   if (trim(argument)=='reject_adjoint') then
      call configure_wall_local_coordinates('wall_local',.true.,.true.,.true.,.true.)
      stop 0 ! WILL_FAIL test must fail if the rejection above did not occur.
   endif
   if (trim(argument)=='reject_nonlocal') then
      call configure_wall_local_coordinates('wall_local',.false.,.true.,.false.,.true.)
      stop 0
   endif

   wx=[0.0_real64,0.1_real64,0.3_real64,0.6_real64,1.0_real64,1.5_real64,2.1_real64]
   wy=1.0_real64; hint=[0.0_real64,0.1_real64]
   call frozen_wall_basis(wx,wy,2,hint,t,normal,s)
   if (any(t/=[1.0_real64,0.0_real64]).or.any(normal/=[0.0_real64,1.0_real64])) &
      error stop 'Exact planar basis was not preserved.'
   do i=1,7
      theta=0.01_real64*(i-1)
      wx(i)=-cos(theta); wy(i)=sin(theta)
   enddo
   call frozen_wall_basis(wx,wy,1,[-0.1_real64,0.0_real64],t,normal,s)
   if (maxval(abs(t-[0.0_real64,1.0_real64]))>1e-7_real64) error stop 'Nose tangent inaccurate.'
   if (maxval(abs(normal-[-1.0_real64,0.0_real64]))>1e-7_real64) error stop 'Nose normal inaccurate.'
   do j=1,16
      n(j)=0.1_real64*(j-1)+0.002_real64*(j-1)**2
   enddo
   q(1,:)=n+2.0_real64; q(2,:)=n*n; q(3,:)=0.0_real64; q(4,:)=3.0_real64; q(5,:)=1.0_real64
   call projected_normal_derivatives(n,q,en,enn,qn,qnn)
   if (maxval(abs(qn(1,:)-1.0_real64))>1e-11_real64) error stop 'Linear profile first derivative failed.'
   if (maxval(abs(qnn(1,:)))>1e-10_real64) error stop 'Linear profile second derivative failed.'
   if (maxval(abs(qn(2,3:14)-2*n(3:14)))>1e-10_real64) error stop 'Quadratic normal derivative failed.'
   if (maxval(abs(qnn(2,3:14)-2.0_real64))>1e-9_real64) error stop 'Quadratic second derivative failed.'

   in=7; jn=16; yb=0; ye=15; rank=0; world_size=1
   allocate(xx(in,jn),yy(in,jn),qq(5,in,jn),drx(in),drz(in))
   allocate(eta_y(in,jn),eta_yy(in,jn),dq0dy(5,in,jn),dq0dyy(5,in,jn))
   drx=0.0_real64; drz=0.0_real64
   do i=1,in
      theta=0.01_real64*(i-1)
      do j=1,jn
         xx(i,j)=wx(i)*(1+n(j)); yy(i,j)=wy(i)*(1+n(j))
         qq(1,i,j)=2.0_real64
         qq(2,i,j)=sin(theta)*n(j)-cos(theta)*0.1_real64*n(j)
         qq(3,i,j)=cos(theta)*n(j)+sin(theta)*0.1_real64*n(j)
         qq(4,i,j)=3.0_real64; qq(5,i,j)=1.0_real64
      enddo
   enddo
   eta_y=1.0_real64; eta_yy=0.0_real64; dq0dy=0.0_real64; dq0dyy=0.0_real64
   call configure_wall_local_coordinates('wall_local',.true.,.true.,.false.,.true.)
   q=qq(:,3,:); qn=0.0_real64; qnn=0.0_real64; en=0.0_real64; enn=0.0_real64
   call wall_local_grid(3,n,en,enn)
   call wall_local_flow(3,q,qn,qnn)
   if (.not.all(ieee_is_finite(qn)).or..not.all(ieee_is_finite(qnn))) &
      error stop 'Curved profile derivatives are not finite.'
   if (maxval(abs(q(2,:)-n))>1e-7_real64.or.maxval(abs(q(3,:)-0.1_real64*n))>1e-7_real64) &
      error stop 'Whole velocity rotation failed.'

   do i=1,in
      xx(i,:)=real(i,real64)
      do j=1,jn
         nvalue=0.1_real64*(j-1)+0.002_real64*(j-1)**2
         yy(i,j)=1.0_real64+nvalue
         qq(:,i,j)=[2.0_real64,nvalue,0.1_real64*nvalue,3.0_real64,1.0_real64]
      enddo
   enddo
   eta_y=7.0_real64; eta_yy=8.0_real64; dq0dy=9.0_real64; dq0dyy=10.0_real64
   call configure_wall_local_coordinates('wall_local',.true.,.true.,.false.,.true.)
   q=qq(:,3,:); expected=q; qn=dq0dy(:,3,:); qnn=dq0dyy(:,3,:)
   en=eta_y(3,:); enn=eta_yy(3,:)
   reference_en=en; reference_enn=enn; reference_qn=qn; reference_qnn=qnn
   call wall_local_grid(3,n,en,enn)
   call wall_local_flow(3,q,qn,qnn)
   if (any(q/=expected).or.any(en/=reference_en).or.any(enn/=reference_enn).or. &
       any(qn/=reference_qn).or.any(qnn/=reference_qnn)) error stop 'Identity flat profile was changed.'
   if (abs(n(1))>0.0_real64) error stop 'Local wall origin was not zero.'
   open(newunit=unit,file='wall_local_profile_metadata.csv',status='old')
   close(unit,status='delete')
   print *, 'wall_local geometry, normal derivatives, rotation, and exact flat compatibility passed.'
end program test_wall_local_profiles
