! Explicit frozen-wall-basis, zero-curvature profile coordinates for direct LST.
! This is a local tangent-plane approximation, not a nonparallel nose model.
module mod_wall_local_coordinates
   use iso_fortran_env, only: real64
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use global_parameters, only: cfg, xx, yy, qq, drx, drz, in, jn, yb, ye, &
                                rank, world_size => size, global_eta_y => eta_y, global_eta_yy => eta_yy, &
                                global_qy => dq0dy, global_qyy => dq0dyy
   use mod_dx, only: dx_fdm4nd_center, dxx_fdm4nd_center
   implicit none
   private
   public :: configure_wall_local_coordinates, wall_local_enabled
   public :: wall_local_grid, wall_local_flow
   public :: frozen_wall_basis, projected_normal_derivatives

   logical, save :: configured=.false., enabled=.false., consistent_normal_derivatives=.false.
   logical, save :: cached_identity=.false., metadata_started=.false.
   integer, save :: cached_station=0
   real(real64), allocatable, save :: wall_s(:), local_n(:), local_eta_n(:), local_eta_nn(:)
   real(real64), allocatable, save :: local_q(:,:), local_qn(:,:), local_qnn(:,:)
   logical, allocatable, save :: reported(:)
   real(real64), save :: tangent(2), normal(2), maximum_tangential_drift

contains

   subroutine configure_wall_local_coordinates(name,direct_local,spatial,adjoint,comoving)
      character(*), intent(in) :: name
      logical, intent(in) :: direct_local,spatial,adjoint,comoving
      character(len=len(name)) :: selection
      integer :: k, code
      selection=trim(adjustl(name))
      do k=1,len_trim(selection)
         code=iachar(selection(k:k))
         if (code>=iachar('A').and.code<=iachar('Z')) selection(k:k)=achar(code+32)
      enddo
      select case(trim(selection))
      case('global')
         enabled=.false.
      case('wall_local','wall_local_consistent')
         if (.not.direct_local.or..not.spatial.or.adjoint.or..not.comoving) &
            error stop 'wall_local requires direct local spatial LST with COMOVING_CS; no adjoint/LPSE support.'
         if (world_size/=1) error stop 'wall_local currently requires one MPI rank.'
         enabled=.true.
      case default
         error stop 'lst.profile_coordinates must be global, wall_local or wall_local_consistent.'
      end select
      consistent_normal_derivatives=trim(selection)=='wall_local_consistent'
      configured=.true.
      cached_station=0
      metadata_started=.false.
      if (allocated(wall_s)) deallocate(wall_s,local_n,local_eta_n,local_eta_nn,local_q,local_qn,local_qnn,reported)
   end subroutine configure_wall_local_coordinates

   logical function wall_local_enabled() result(active)
      character(256) :: selection
      if (.not.configured) then
         ! Other executables share these profile modules, but their old adapters
         ! cannot silently enable a transform lacking a validated direct-LST context.
         selection='global'
         if (cfg%has_key('lst','profile_coordinates')) call cfg%get('lst','profile_coordinates',selection)
         if (trim(selection)/='global') &
            error stop 'profile_coordinates requested outside the supported direct spatial START_LST adapter.'
      endif
      active=enabled
   end function wall_local_enabled

   subroutine frozen_wall_basis(wall_x,wall_y,ix,outward_hint,t,n,s)
      real(real64), intent(in) :: wall_x(:),wall_y(:),outward_hint(2)
      integer, intent(in) :: ix
      real(real64), intent(out) :: t(2),n(2),s(:)
      real(real64) :: a(2),b(2),cross,scale,z(5),weights(5),denom,product,value
      integer :: count,left,right,k,m,l,total
      logical :: straight
      total=size(wall_x)
      if (total<3.or.size(wall_y)/=total.or.size(s)/=total) &
         error stop 'wall_local requires at least three wall geometry points.'
      if (ix<1.or.ix>total) error stop 'wall_local station is outside wall geometry.'
      s(1)=0.0_real64
      do k=2,total
         scale=sqrt((wall_x(k)-wall_x(k-1))**2+(wall_y(k)-wall_y(k-1))**2)
         if (.not.ieee_is_finite(scale).or.scale<=tiny(scale)) &
            error stop 'wall_local wall geometry has coincident or nonfinite points.'
         s(k)=s(k-1)+scale
      enddo
      left=max(1,min(ix-1,total-2))
      a=[wall_x(left+1)-wall_x(left),wall_y(left+1)-wall_y(left)]
      b=[wall_x(left+2)-wall_x(left+1),wall_y(left+2)-wall_y(left+1)]
      cross=a(1)*b(2)-a(2)*b(1)
      straight=abs(cross)<=64*epsilon(1.0_real64)*sqrt(sum(a*a)*sum(b*b))
      if (straight) then
         ! Preserve exact flat columns, including the first point after a junction.
         t=a+b
      else
         count=min(5,total)
         left=max(1,min(ix-2,total-count+1))
         right=left+count-1
         scale=maxval(abs(s(left:right)-s(ix)))
         z(1:count)=(s(left:right)-s(ix))/scale
         weights=0.0_real64
         ! Derivative of the local Lagrange interpolant at s(ix); <=5 points.
         do k=1,count
            denom=1.0_real64
            do l=1,count
               if (l/=k) denom=denom*(z(k)-z(l))
            enddo
            value=0.0_real64
            do m=1,count
               if (m==k) cycle
               product=1.0_real64
               do l=1,count
                  if (l/=k.and.l/=m) product=product*(-z(l))
               enddo
               value=value+product
            enddo
            weights(k)=value/denom/scale
         enddo
         t=[sum(weights(1:count)*wall_x(left:right)),sum(weights(1:count)*wall_y(left:right))]
      endif
      scale=sqrt(sum(t*t))
      if (.not.ieee_is_finite(scale).or.scale<=tiny(scale)) error stop 'wall_local wall tangent is invalid.'
      t=t/scale
      where(abs(t)<32*epsilon(1.0_real64)) t=0.0_real64
      t=t/sqrt(sum(t*t))
      n=[-t(2),t(1)]
      if (sum(n*outward_hint)<0.0_real64) n=-n
      if (abs(sum(n*outward_hint))<=tiny(scale)) error stop 'wall_local outward grid direction is tangential.'
   end subroutine frozen_wall_basis

   subroutine projected_normal_derivatives(n,q,eta_n,eta_nn,qn,qnn)
      real(real64), intent(in) :: n(:),q(:,:)
      real(real64), intent(out) :: eta_n(:),eta_nn(:),qn(:,:),qnn(:,:)
      real(real64) :: nj(size(n)),njj(size(n)),qj(size(n)),qjj(size(n))
      integer :: count,k
      count=size(n)
      if (count<5.or.size(q,2)/=count) error stop 'wall_local requires at least five normal points.'
      if (.not.all(ieee_is_finite(n)).or..not.all(ieee_is_finite(q))) &
         error stop 'wall_local normal coordinates or baseflow contain nonfinite values.'
      if (any(n(2:count)<=n(1:count-1))) error stop 'wall_local projected normal coordinate is not increasing.'
      call dx_fdm4nd_center(n,1,count,nj,1,count)
      call dxx_fdm4nd_center(n,1,count,njj,1,count)
      if (any(nj<=0.0_real64).or..not.all(ieee_is_finite(nj))) &
         error stop 'wall_local normal coordinate derivative is not positive.'
      eta_n=1.0_real64/nj
      eta_nn=-njj/nj**3
      do k=1,size(q,1)
         call dx_fdm4nd_center(q(k,:),1,count,qj,1,count)
         call dxx_fdm4nd_center(q(k,:),1,count,qjj,1,count)
         qn(k,:)=qj*eta_n
         qnn(k,:)=qjj*eta_n**2+qj*eta_nn
      enddo
      if (.not.all(ieee_is_finite(eta_n)).or..not.all(ieee_is_finite(eta_nn)).or. &
          .not.all(ieee_is_finite(qn)).or..not.all(ieee_is_finite(qnn))) &
         error stop 'wall_local recomputed normal derivatives are not finite.'
   end subroutine projected_normal_derivatives

   subroutine prepare_wall_local(ix)
      integer, intent(in) :: ix
      integer :: k,unit,status
      real(real64) :: hint(2),dx(size(yy,2)),dy(size(yy,2)),flat_scale
      if (.not.wall_local_enabled()) error stop 'wall_local preparation called while disabled.'
      if (cached_station==ix) return
      if (world_size/=1.or.yb/=0.or.ye/=jn-1) &
         error stop 'wall_local requires one rank and the full wall-normal domain.'
      if (any(abs(drx)>1e-14_real64).or.any(abs(drz)>1e-14_real64)) &
         error stop 'wall_local tangent-plane model requires zero explicit curvature.'
      if (.not.allocated(wall_s)) then
         allocate(wall_s(in),local_n(jn),local_eta_n(jn),local_eta_nn(jn))
         allocate(local_q(5,jn),local_qn(5,jn),local_qnn(5,jn),reported(in))
         reported=.false.
      endif
      hint=[xx(ix,2)-xx(ix,1),yy(ix,2)-yy(ix,1)]
      call frozen_wall_basis(xx(:,1),yy(:,1),ix,hint,tangent,normal,wall_s)
      dx=xx(ix,:)-xx(ix,1)
      dy=yy(ix,:)-yy(ix,1)
      local_n=normal(1)*dx+normal(2)*dy
      local_n(1)=0.0_real64
      maximum_tangential_drift=maxval(abs(tangent(1)*dx+tangent(2)*dy))
      local_q=qq(:,ix,:)
      local_q(2,:)=tangent(1)*qq(2,ix,:)+tangent(2)*qq(3,ix,:)
      local_q(3,:)=normal(1)*qq(2,ix,:)+normal(2)*qq(3,ix,:)
      call projected_normal_derivatives(local_n,local_q,local_eta_n,local_eta_nn,local_qn,local_qnn)
      flat_scale=max(1.0_real64,maxval(abs(xx(ix,:))))
      cached_identity=tangent(1)==1.0_real64.and.tangent(2)==0.0_real64.and. &
                      normal(1)==0.0_real64.and.normal(2)==1.0_real64.and. &
                      maxval(abs(dx))<=16*epsilon(flat_scale)*flat_scale.and..not.consistent_normal_derivatives
      if (cached_identity) then
         local_eta_n=global_eta_y(ix,:)
         local_eta_nn=global_eta_yy(ix,:)
         local_qn=global_qy(:,ix,:)
         local_qnn=global_qyy(:,ix,:)
         if (.not.all(ieee_is_finite(local_eta_n)).or..not.all(ieee_is_finite(local_eta_nn)).or. &
             .not.all(ieee_is_finite(local_qn)).or..not.all(ieee_is_finite(local_qnn))) &
            error stop 'wall_local identity flat profile has nonfinite original normal derivatives.'
      endif
      cached_station=ix
      if (.not.reported(ix).and.rank==0) then
         if (.not.metadata_started) then
            open(newunit=unit,file='wall_local_profile_metadata.csv',status='replace',action='write',iostat=status)
            if (status/=0) error stop 'Cannot create wall_local profile metadata.'
            write(unit,'(a)') '# wall_local: frozen wall basis; zero explicit curvature; direct spatial local LST only.'
            write(unit,'(a)') '# x_global is unchanged in standard solver tables; s_polyline is independent metadata.'
            write(unit,'(a)') '# Tangential baseflow derivatives are unavailable and set to zero for the frozen profile.'
            write(unit,'(a)') 'iloc,x_global_N,y_global_N,s_polyline_N,t_x,t_y,n_x,n_y,n_min,n_max,'// &
               'max_tangential_drift,max_abs_Vn,max_abs_Qn,max_abs_Qnn,min_delta_n,identity_flat,normal_derivatives_finite'
            close(unit)
            metadata_started=.true.
         endif
         open(newunit=unit,file='wall_local_profile_metadata.csv',status='old',position='append',action='write',iostat=status)
         if (status/=0) error stop 'Cannot append wall_local profile metadata.'
         write(unit,'(i0,14(",",es25.16e3),",",l1,",T")') ix,xx(ix,1),yy(ix,1),wall_s(ix), &
            tangent,normal,minval(local_n),maxval(local_n),maximum_tangential_drift, &
            maxval(abs(local_q(3,:))),maxval(abs(local_qn)),maxval(abs(local_qnn)), &
            minval(local_n(2:jn)-local_n(1:jn-1)),cached_identity
         close(unit)
         reported(ix)=.true.
      endif
   end subroutine prepare_wall_local

   subroutine wall_local_grid(ix,n,eta_n,eta_nn)
      integer, intent(in) :: ix
      real(real64), intent(out) :: n(:)
      real(real64), intent(inout) :: eta_n(:),eta_nn(:)
      call prepare_wall_local(ix)
      n=local_n
      ! The original Cartesian metrics are already the exact desired discrete
      ! normal metrics on identity flat columns. Keep them for regression.
      if (.not.cached_identity) then
         eta_n=local_eta_n
         eta_nn=local_eta_nn
      endif
   end subroutine wall_local_grid

   subroutine wall_local_flow(ix,q,qn,qnn)
      integer, intent(in) :: ix
      real(real64), intent(inout) :: q(:,:),qn(:,:),qnn(:,:)
      call prepare_wall_local(ix)
      if (.not.cached_identity) then
         q=local_q
         qn=local_qn
         qnn=local_qnn
      endif
   end subroutine wall_local_flow
end module mod_wall_local_coordinates
