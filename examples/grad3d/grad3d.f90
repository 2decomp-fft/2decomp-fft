!! SPDX-License-Identifier: BSD-3-Clause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the gradient of a periodic field using global
! transposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program grad3d

   use mpi

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use MPI
#if defined(_GPU)
   use cudafor
   use openacc
#endif

   implicit none

   integer, parameter :: nx_base = 64, ny_base = 64, nz_base = 64
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer :: nargin, arg, FNLength, status, DecInd
   character(len=80) :: InputFN

   real(mytype), parameter :: lx = 1.0_mytype
   real(mytype), parameter :: ly = 1.0_mytype
   real(mytype), parameter :: lz = 1.0_mytype

   real(mytype) :: dx, dy, dz

   real(mytype), allocatable, dimension(:, :, :) :: phi1, phi2, phi3
   real(mytype), allocatable, dimension(:, :, :) :: dphiX, dphiY, dphiz
   real(mytype), allocatable, dimension(:, :, :) :: wk2, wk3

   integer :: ierror

   integer :: xlast, ylast, zlast


   logical :: passing, all_pass

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   ! Now we can check if user put some inputs
   ! Handle input file like a boss -- GD
   nargin = command_argument_count()
   if ((nargin == 0) .or. (nargin == 2) .or. (nargin == 5)) then
      do arg = 1, nargin
         call get_command_argument(arg, InputFN, FNLength, status)
         read (InputFN, *, iostat=status) DecInd
         if (arg == 1) then
            p_row = DecInd
         elseif (arg == 2) then
            p_col = DecInd
         elseif (arg == 3) then
            nx = DecInd
         elseif (arg == 4) then
            ny = DecInd
         elseif (arg == 5) then
            nz = DecInd
         end if
      end do
   else
      ! nrank not yet computed we need to avoid write
      ! for every rank
      call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)
      if (nrank == 0) then
         print *, "This Test takes no inputs or 2 inputs as"
         print *, "  1) p_row (default=0)"
         print *, "  2) p_col (default=0)"
         print *, "or 5 inputs as"
         print *, "  1) p_row (default=0)"
         print *, "  2) p_col (default=0)"
         print *, "  3) nx "
         print *, "  4) ny "
         print *, "  5) nz "
         print *, "Number of inputs is not correct and the defult settings"
         print *, "will be used"
      end if
   end if
   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   dx = lx/real(nx,mytype)
   dy = ly/real(ny,mytype)
   dz = lz/real(nz,mytype)

   call allocate_var()
   !$acc data create(phi2,phi3,wk2,wk3) copy(phi1,dphiX, dphiY, dphiZ)
   call init_phi()
   call compute_grad()
   call test_derX(dphiX)
   call test_derY()
   call test_derZ()
   !$acc end data

   if (nrank == 0) then
      write (*, *) '-----------------------------------------------'
      write (*, *) "End GRAD calculation check all pass: ", all_pass
      write (*, *) '==============================================='
   end if

   call finalize()

   call decomp_2d_finalize

   if (.not. all_pass) call decomp_2d_abort(1, "Error in grad3d")

   call MPI_FINALIZE(ierror)

contains

   !=====================================================================
   ! Initialize
   !=====================================================================
   subroutine allocate_var()

      use decomp_2d

      implicit none

      logical, parameter :: glob = .false.

      ! initialise u,v,w with random numbers in X-pencil
      call alloc_x(phi1, glob)
      call alloc_x(dphiX, glob)
      call alloc_x(dphiY, glob)
      call alloc_x(dphiZ, glob)

      ! Working array used more than once
      call alloc_y(phi2, glob)
      call alloc_y(wk2, glob)

      call alloc_z(phi3, glob)
      call alloc_z(wk3, glob)

      all_pass = .true.

   end subroutine allocate_var

   !=====================================================================
   ! Initialize the scalar field
   !=====================================================================
   subroutine init_phi()

      implicit none

      integer :: i, j, k
      real(mytype), parameter :: twopi=2._mytype*acos(-1._mytype)
      real(mytype) :: x, y, z
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Initial velocity
      !$acc kernels default(present)
      do concurrent (k=1:xe3, j=1:xe2, i=1:xe1)
         z = (k + xs3 - 2) * dz
         y = (j + xs2 - 2) * dy
         x = (i + xs1 - 2) * dx
         phi1(i, j, k) = -2._mytype * cos(twopi * (x / lx)) * cos(twopi * (y / ly)) * sin(twopi * (z / lz))
      enddo
      !$acc end kernels

   end subroutine init_phi
   !=====================================================================
   ! Finalize with deallocation of arrays
   !=====================================================================
   subroutine finalize()

      implicit none

      deallocate (phi1, phi2, phi3)
      deallocate (dphiX, dphiY, dphiZ)
      deallocate (wk2, wk3)

   end subroutine finalize
   !=====================================================================
   ! Calculate gradient using global transposition
   !=====================================================================
   subroutine compute_grad()

      implicit none

      ! Compute X derivative
      call derx(dphiX,phi1,dx,xsize(1),xsize(2),xsize(3))

      ! Compute Y derivative
      call transpose_x_to_y(phi1,phi2)
      call dery(wk2,phi2,dy,ysize(1),ysize(2),ysize(3))
      call transpose_y_to_x(wk2,dphiY)

      ! Compute Z derivative
      call transpose_y_to_z(phi2,phi3)
      call derx(wk3,phi3,dz,zsize(1),zsize(2),zsize(3))
      call transpose_z_to_y(wk3,wk2)
      call transpose_y_to_x(wk2,dphiZ)

   end subroutine compute_grad
   !=====================================================================
   ! Calculate gradient in X (data in X-pencil)
   !=====================================================================
   subroutine derx(df,ff,delta,nx,ny,nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx,ny,nz) :: df
      real(mytype), intent(in), dimension(nx,ny,nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff/delta

      !$acc kernels default(present)
      do concurrent (k=1:nz, j=1:ny)
         df(1,j,k) = coeff*(ff(2,j,k)-ff(nx,j,k)) 
         do concurrent (i=2:nx-1)
            df(i,j,k) = coeff*(ff(i+1,j,k)-ff(i-1,j,k)) 
         enddo
         df(nx,j,k) = coeff*(ff(1,j,k)-ff(nx-1,j,k))
      enddo
      !$acc end kernels

   end subroutine derx 
   !=====================================================================
   ! Calculate gradient in Y (data in Y-pencil)
   !=====================================================================
   subroutine dery(df,ff,delta,nx,ny,nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx,ny,nz) :: df
      real(mytype), intent(in), dimension(nx,ny,nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff/delta

      !$acc kernels default(present)
      do concurrent (k=1:nz)
        do concurrent (i=1:nx)
           df(i,1,k) = coeff*(ff(i,2,k)-ff(i,ny,k))
        enddo
        do concurrent (j=2:ny-1, i=1:nx)
           df(i,j,k) = coeff*(ff(i,j+1,k)-ff(i,j-1,k))
        enddo
        do concurrent (i=1:nx)
           df(i,ny,k) = coeff*(ff(i,1,k)-ff(i,ny-1,k)) 
        enddo
      enddo
      !$acc end kernels

   end subroutine dery 
   !=====================================================================
   ! Calculate gradient in Z (data in Z-pencil)
   !=====================================================================
   subroutine derz(df,ff,delta,nx,ny,nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx,ny,nz) :: df
      real(mytype), intent(in), dimension(nx,ny,nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype
      
      coeff = coeff/delta

      !$acc kernels default(present)
      do concurrent (j=1:ny, i=1:nx)
         df(i,j,1) = coeff*(ff(i,j,2)-ff(i,j,nz  ))
      enddo
      do concurrent (k=2:nz-1, j=1:ny, i=1:nx)
         df(i,j,k) = coeff*(ff(i,j,k+1)-ff(i,j,k-1))
      enddo
      do concurrent (j=1:ny, i=1:nx)
         df(i,j,nz) = coeff*(ff(i,j,1)-ff(i,j,nz-1))
      enddo
      !$acc end kernels

   end subroutine derz 
   !=====================================================================
   ! Calculate divergence using halo-cell exchange (data in X-pencil)
   !=====================================================================
   subroutine test_derX(df)

      implicit none
      ! Arguments
      real(mytype), intent(out), dimension(xsize(1),xsize(2),xsize(3)) :: df

      integer :: i, j, k
      real(mytype), parameter :: twopi=2._mytype*acos(-1._mytype)
      real(mytype) :: x, y, z
      real(mytype) :: dphi, dphi_num
      real(mytype) :: error = 0._mytype
      real(mytype) :: err_all = 0._mytype
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Initial velocity
      !$acc parallel loop default(present) reduction(+:error)
      do k=1,xe3
        do j=1,xe2
          do i=1,xe1
             z = (k + xs3 - 2) * dz
             y = (j + xs2 - 2) * dy
             x = (i + xs1 - 2) * dx
             dphi = -2._mytype * twopi * sin (twopi * (x / lx)) * cos(twopi * (y / ly)) * sin(twopi * (z / lz))
             dphi_num = df(i,j,k)
             error = error + (dphi-dphi_num)*(dphi-dphi_num)
          enddo
        enddo
      enddo
      !$acc end loop
      call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      err_all = sqrt(err_all) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

      if (nrank == 0) then
         write (*, *) 'DX error / mesh point: ', err_all
      endif

   end subroutine test_derX

   !=====================================================================
   ! Calculate divergence using halo-cell exchange (data in Y-pencil)
   !=====================================================================
   subroutine test_derY()

      implicit none


   end subroutine test_derY

   !=====================================================================
   ! Calculate divergence using halo-cell exchange (data in Z-pencil)
   !=====================================================================
   subroutine test_derZ()

      implicit none

   end subroutine test_derZ

end program grad3d
