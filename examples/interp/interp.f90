!! SPDX-License-Identifier: BSD-3-Clause
program interp

   use mpi
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_interp
   use decomp_2d_mpi
   use decomp_2d_testing
#if defined(_GPU)
   use cudafor
   use openacc
#endif

   implicit none

   integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot

   ! 3D arrays on the main grid
   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3

   ! decomp_info object for the coarse grid
   type(decomp_info) :: grid
   ! 3D arrays on the coarse grid
   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b

   real(mytype), parameter :: eps = 1.0E-7_mytype

   integer :: i, j, k, ierror
   integer :: st1, st2, st3
   integer :: en1, en2, en3

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   ! Now we can check if user put some inputs
   call decomp_2d_testing_init(p_row, p_col, nx, ny, nz)

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   ! Default grid. Allocate arrays
   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)

   ! Init x
   st1 = xstart(1); en1 = xend(1)
   st2 = xstart(2); en2 = xend(2)
   st3 = xstart(3); en3 = xend(3)
   ! original x-pencil based data
   !$acc data copy(u1)
   !$acc parallel loop default(present)
   do k = st3, en3
      do j = st2, en2
         do i = st1, en1
            u1(i, j, k) = real(i, mytype) + real(j-1, mytype) + real(k-1, mytype)
         end do
      end do
   end do
   !$acc end loop
   !$acc update self(u1)
   !$acc end data

   ! Init y
   st1 = ystart(1); en1 = yend(1)
   st2 = ystart(2); en2 = yend(2)
   st3 = ystart(3); en3 = yend(3)
   ! original x-pencil based data
   !$acc data copy(u2)
   !$acc parallel loop default(present)
   do k = st3, en3
      do j = st2, en2
         do i = st1, en1
            u2(i, j, k) = real(i-1, mytype) + real(j, mytype) + real(k-1, mytype)
         end do
      end do
   end do
   !$acc end loop
   !$acc update self(u2)
   !$acc end data

   ! Init z
   st1 = zstart(1); en1 = zend(1)
   st2 = zstart(2); en2 = zend(2)
   st3 = zstart(3); en3 = zend(3)
   ! original x-pencil based data
   !$acc data copy(u3)
   !$acc parallel loop default(present)
   do k = st3, en3
      do j = st2, en2
         do i = st1, en1
            u3(i, j, k) = real(i-1, mytype) + real(j-1, mytype) + real(k, mytype)
         end do
      end do
   end do
   !$acc end loop
   !$acc update self(u3)
   !$acc end data

   ! Print 1D grids
   call show(decomp_main, u1, u2, u3)

   ! Interpolate on a slightly smaller grid
   call decomp_info_init(nx-1, ny-1, nz-1, grid)
   ! Allocate memory
   call alloc_x(u1b, grid)
   call alloc_y(u2b, grid)
   call alloc_z(u3b, grid)
   ! Interpolation from default "decomp_main" to "grid"
   !$acc data copyin(u1, u2, u3), copy(u1b, u2b, u3b)
   call decomp_2d_interp_var3d(1, u1, u1b, grid)
   call decomp_2d_interp_var3d(2, u2, u2b, grid)
   call decomp_2d_interp_var3d(3, u3, u3b, grid)
   !$acc update self(u1b, u2b, u3b)
   !$acc end data
   ! Print 1D grids
   call show(grid, u1b, u2b, u3b)
   ! Free memory
   deallocate (u1b, u2b, u3b)
   call decomp_info_finalize(grid)

   ! Interpolate on a coarse grid
   call decomp_info_init(max(nx/10, p_row), &
                         max(ny/10, max(p_row, p_col)), &
                         max(nz/10, p_col), &
                         grid)
   ! Allocate memory
   call alloc_x(u1b, grid)
   call alloc_y(u2b, grid)
   call alloc_z(u3b, grid)
   ! Interpolation from default "decomp_main" to "grid"
   !$acc data copyin(u1, u2, u3), copy(u1b, u2b, u3b)
   call decomp_2d_interp_var3d(1, u1, u1b, grid)
   call decomp_2d_interp_var3d(2, u2, u2b, grid)
   call decomp_2d_interp_var3d(3, u3, u3b, grid)
   !$acc update self(u1b, u2b, u3b)
   !$acc end data
   ! Print 1D grids
   call show(grid, u1b, u2b, u3b)
   ! Free memory
   deallocate (u1b, u2b, u3b)
   call decomp_info_finalize(grid)

   ! Free memory
   deallocate (u1, u2, u3)

   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

contains

   subroutine show(info, x, y, z)

      implicit none

      type(decomp_info), intent(in) :: info
      real(mytype), dimension(:,:,:), intent(in) :: x, y, z

      integer :: ierr

      ! Only the initial array will have the x grid located at j=1 and k=1
      if (info%xst(2) == 1 .and. info%xst(3) == 1) then
         write(*,*) "x data at j=1 and k=1"
         write(*,*) "  Size ", info%xsz(1)
         write(*,*) real(x(:,1,1), 4)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! Only the initial array will have the y grid located at i=1 and k=1
      if (info%yst(1) == 1 .and. info%yst(3) == 1) then
         write(*,*) "y data at i=1 and k=1"
         write(*,*) "  Size ", info%ysz(2)
         write(*,*) real(y(1,:,1), 4)
      end if                 
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      ! Only the initial array will have the z grid located at i=1 and j=1
      if (info%zst(1) == 1 .and. info%zst(2) == 1) then
         write(*,*) "z data at i=1 and j=1"
         write(*,*) "  Size ", info%zsz(3)
         write(*,*) real(z(1,1,:), 4)
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

   end subroutine show

end program interp
