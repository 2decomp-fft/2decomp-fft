!! SPDX-License-Identifier: BSD-3-Clause
!!
!! FIXME The issue below is specific to GPU and should be discussed in a dedicated github issue
!!
!! NB in case of GPU only the writing in the aligned pencil (i.e. X for a 1 array) is performed.
!! IO subrotines needs update for non managed GPU case
!!
program io_plane_test

   use mpi
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io
   use decomp_2d_io_adios
   use decomp_2d_io_object_adios
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

   real(mytype), allocatable, dimension(:, :, :) :: data1
   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3

   integer :: i, j, k, m, ierror
   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3

   type(d2d_io_adios) :: io

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

   call decomp_2d_register_plane("x_pencil-x_plane.dat", 1, real_type, opt_reduce_prec=.false.)
   call decomp_2d_register_plane("y_pencil-y_plane.dat", 2, real_type, opt_reduce_prec=.false.)
   call decomp_2d_register_plane("z_pencil-z_plane.dat", 3, real_type, opt_reduce_prec=.false.)

   ! ***** global data *****
   allocate (data1(nx, ny, nz))
   m = 1
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            data1(i, j, k) = real(m, mytype)
            m = m + 1
         end do
      end do
   end do

   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)

   ! For GPU we port the global data create the different pencil arrays
   ! Move back to host the arrays for writing on disk
   xst1 = xstart(1); xen1 = xend(1)
   xst2 = xstart(2); xen2 = xend(2)
   xst3 = xstart(3); xen3 = xend(3)

   !$acc data copyin(data1) copy(u1,u2,u3)
   ! original X-pensil based data
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1(i, j, k) = data1(i, j, k)
         end do
      end do
   end do
   !$acc end loop
   call transpose_x_to_y(u1, u2)
   call transpose_y_to_z(u2, u3)
   !$acc update self(u1)
   !$acc update self(u2)
   !$acc update self(u3)
   !$acc end data

   ! X-pencil data
   call io%open_start(decomp_2d_write_mode)
   call decomp_2d_adios_write_plane(io, u1, 'x_pencil-x_plane.dat', &
                                    opt_iplane=nx / 2, &
                                    opt_ipencil=1, &
                                    opt_reduce_prec=.false.)
   ! Y-pencil data
   call decomp_2d_adios_write_plane(io, u2, 'y_pencil-y_plane.dat', &
                                    opt_iplane=ny / 2, &
                                    opt_ipencil=2, &
                                    opt_reduce_prec=.false.)
   ! Z-pencil data
   call decomp_2d_adios_write_plane(io, u3, 'z_pencil-z_plane.dat', &
                                    opt_iplane=nz / 2, &
                                    opt_ipencil=3, &
                                    opt_reduce_prec=.false.)
   call io%end_close

   deallocate (u1, u2, u3)
   deallocate (data1)
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_plane_test
