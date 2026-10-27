!! SPDX-License-Identifier: BSD-3-Clause
program io_subdomain

   use mpi
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io
   use decomp_2d_io_object_mpi
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

#ifdef COMPLEX_TEST
   complex(mytype), allocatable, dimension(:, :, :) :: data1
   complex(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
#else
   real(mytype), allocatable, dimension(:, :, :) :: data1
   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   integer :: i, j, k, m, ierror
   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3

   logical :: dir_exists

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

   call decomp_2d_io_init()

   ! ***** global data *****
   allocate (data1(nx, ny, nz))
   m = 1
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
#ifdef COMPLEX_TEST
            data1(i, j, k) = cmplx(real(m, mytype), real(nx * ny * nz - m, mytype))
#else
            data1(i, j, k) = real(m, mytype)
#endif
            m = m + 1
         end do
      end do
   end do

   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)

   xst1 = xstart(1); xen1 = xend(1)
   xst2 = xstart(2); xen2 = xend(2)
   xst3 = xstart(3); xen3 = xend(3)
   ! original x-pencil based data
   !$acc data copyin(data1) copy(u1,u2,u3)
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1(i, j, k) = data1(i, j, k)
         end do
      end do
   end do
   !$acc end loop

   ! transpose
   call transpose_x_to_y(u1, u2)
   call transpose_y_to_z(u2, u3)
   !$acc update self(u1)
   !$acc update self(u2)
   !$acc update self(u3)
   !$acc end data

   ! write to disk
   if (nrank == 0) then
      inquire (file="out", exist=dir_exists)
      if (.not. dir_exists) then
         call execute_command_line("mkdir out 2> /dev/null")
      end if
   end if

   ! Write total array
   call decomp_2d_write_one(1, u1, 'ref_one_1.dat', opt_dirname='out')
   call decomp_2d_write_one(2, u2, 'ref_one_2.dat', opt_dirname='out')
   call decomp_2d_write_one(3, u3, 'ref_one_3.dat', opt_dirname='out')

   ! Write total array using subdomain
   call decomp_2d_write_subdomain(1, u1, 'one_1.dat', (/1, 1, 1/), (/nx, ny, nz/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'one_2.dat', (/1, 1, 1/), (/nx, ny, nz/), opt_dirname='out')
   call decomp_2d_write_subdomain(3, u3, 'one_3.dat', (/1, 1, 1/), (/nx, ny, nz/), opt_dirname='out')

   ! Check
   if (nrank == 0) then
      call execute_command_line("diff out/ref_one_1.dat out/one_1.dat", exitstat=i)
      call execute_command_line("diff out/ref_one_2.dat out/one_2.dat", exitstat=j)
      call execute_command_line("diff out/ref_one_3.dat out/one_3.dat", exitstat=k)
      if (i /= 0 .or. j /=0 .or. k /= 0) then
         print *, "Subdomain failed to reproduce write_one"
         stop 1
      end if
   end if

   ! Write planes
   call decomp_2d_write_plane(1, u1, 'ref_plane_1.dat', opt_iplane=nx / 2, opt_dirname='out', opt_reduce_prec=.false.)
   call decomp_2d_write_plane(2, u2, 'ref_plane_2.dat', opt_iplane=ny / 2, opt_dirname='out', opt_reduce_prec=.false.)
   call decomp_2d_write_plane(3, u3, 'ref_plane_3.dat', opt_iplane=nz / 2, opt_dirname='out', opt_reduce_prec=.false.)

   ! Write planes using subdomain
   call decomp_2d_write_subdomain(1, u1, 'plane_1.dat', (/nx/2, 1, 1/), (/nx/2, ny, nz/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'plane_2.dat', (/1, ny/2, 1/), (/nx, ny/2, nz/), opt_dirname='out')
   call decomp_2d_write_subdomain(3, u3, 'plane_3.dat', (/1, 1, nz/2/), (/nx, ny, nz/2/), opt_dirname='out')

   ! Check
   if (nrank == 0) then
      call execute_command_line("diff out/ref_plane_1.dat out/plane_1.dat", exitstat=i)
      call execute_command_line("diff out/ref_plane_2.dat out/plane_2.dat", exitstat=j)
      call execute_command_line("diff out/ref_plane_3.dat out/plane_3.dat", exitstat=k)
      if (i /= 0 .or. j /=0 .or. k /= 0) then
         print *, "Subdomain failed to reproduce write_plane"
         stop 2
      end if
   end if

   ! Split the main domain into blocks
   call decomp_2d_write_subdomain(2, u2, 'x1_y1_z1.dat', (/1, 1, 1/), (/nx/2, ny/2, nz/2/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'x1_y1_z2.dat', (/1, 1, nz/2+1/), (/nx/2, ny/2, nz/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'x1_y2_z1.dat', (/1, ny/2+1, 1/), (/nx/2, ny, nz/2/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'x1_y2_z2.dat', (/1, ny/2+1, nz/2+1/), (/nx/2, ny, nz/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'x2_y1_z1.dat', (/nx/2+1, 1, 1/), (/nx, ny/2, nz/2/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'x2_y1_z2.dat', (/nx/2+1, 1, nz/2+1/), (/nx, ny/2, nz/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'x2_y2_z1.dat', (/nx/2+1, ny/2+1, 1/), (/nx, ny, nz/2/), opt_dirname='out')
   call decomp_2d_write_subdomain(2, u2, 'x2_y2_z2.dat', (/nx/2+1, ny/2+1, nz/2+1/), (/nx, ny, nz/), opt_dirname='out')

   ! Check ?
   ! Reading subdomain is not possible currently

   deallocate (u1, u2, u3)
   deallocate (data1)

   call decomp_2d_io_fin
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_subdomain
