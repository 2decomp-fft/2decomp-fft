!! SPDX-License-Identifier: BSD-3-Clause
program io_test

   use mpi
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io
   use decomp_2d_io_object_mpi
   use decomp_2d_mpi
   use decomp_2d_profiler
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

   integer :: ierr

#ifdef COMPLEX_TEST
   complex(mytype), allocatable, dimension(:, :, :) :: data1

   complex(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#else
   real(mytype), allocatable, dimension(:, :, :) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   type(d2d_io_mpi), save :: io1, io2, io3
   integer, save :: req1, req2, req3

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

   ! Make sure the initial value of the MPI request is valid
   req1 = MPI_REQUEST_NULL
   req2 = MPI_REQUEST_NULL
   req3 = MPI_REQUEST_NULL

   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)

   call alloc_x(u1b, .true.)
   call alloc_y(u2b, .true.)
   call alloc_z(u3b, .true.)

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

   !
   call decomp_2d_write_one(1, u1, 'u1.dat', opt_reduce_prec=.false., opt_dirname='out', opt_nb_req=req1, opt_io=io1)
   call decomp_2d_write_one(2, u2, 'u2.dat', opt_reduce_prec=.false., opt_dirname='out', opt_nb_req=req2, opt_io=io2)
   call decomp_2d_write_one(3, u3, 'u3.dat', opt_reduce_prec=.false., opt_dirname='out', opt_nb_req=req3, opt_io=io3)

   !
   ! Some computation can happen here as long as u1, u2 and u3 are not modified
   !

   call MPI_Barrier(MPI_COMM_WORLD, ierr)
   if (decomp_profiler_io) call decomp_profiler_start("io_wait")
   call decomp_2d_mpi_wait(req1)
   if (decomp_profiler_io) call decomp_profiler_end("io_wait")
   call io1%close
   if (decomp_profiler_io) call decomp_profiler_start("io_wait")
   call decomp_2d_mpi_wait(req2)
   if (decomp_profiler_io) call decomp_profiler_end("io_wait")
   call io2%close
   if (decomp_profiler_io) call decomp_profiler_start("io_wait")
   call decomp_2d_mpi_wait(req3)
   if (decomp_profiler_io) call decomp_profiler_end("io_wait")
   call io3%close

   ! read back to different arrays
   u1b = 0; u2b = 0; u3b = 0
   call decomp_2d_read_one(1, u1b, 'u1.dat', opt_reduce_prec=.false., opt_dirname='out')
   call decomp_2d_read_one(2, u2b, 'u2.dat', opt_reduce_prec=.false., opt_dirname='out')
   call decomp_2d_read_one(3, u3b, 'u3.dat', opt_reduce_prec=.false., opt_dirname='out')

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   ! compare
   call check("one file, multiple fields")

   deallocate (u1, u2, u3)
   deallocate (u1b, u2b, u3b)
   deallocate (data1)

   call decomp_2d_io_fin
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

contains

   subroutine check(stage)

      character(len=*), intent(in) :: stage

      integer :: ierr

      if (nrank == 0) then
         print *, "Checking "//stage
      end if
      call MPI_Barrier(MPI_COMM_WORLD, ierr)

      do k = xstart(3), xend(3)
         do j = xstart(2), xend(2)
            do i = xstart(1), xend(1)
               if (abs((u1(i, j, k) - u1b(i, j, k))) > eps) then
                  print *, i, j, k, u1(i, j, k), u1b(i, j, k)
                  call decomp_2d_abort(1, "x_plane")
               end if
            end do
         end do
      end do

      do k = ystart(3), yend(3)
         do j = ystart(2), yend(2)
            do i = ystart(1), yend(1)
               if (abs((u2(i, j, k) - u2b(i, j, k))) > eps) then
                  print *, i, j, k, u2(i, j, k), u2b(i, j, k)
                  call decomp_2d_abort(2, "y_plane")
               end if
            end do
         end do
      end do

      do k = zstart(3), zend(3)
         do j = zstart(2), zend(2)
            do i = zstart(1), zend(1)
               if (abs((u3(i, j, k) - u3b(i, j, k))) > eps) then
                  print *, i, j, k, u3(i, j, k), u3b(i, j, k)
                  call decomp_2d_abort(3, "z_plane")
               end if
            end do
         end do
      end do

      ! Also check against the global data array
      do k = xstart(3), xend(3)
         do j = xstart(2), xend(2)
            do i = xstart(1), xend(1)
               if (abs(data1(i, j, k) - u1b(i, j, k)) > eps) then
                  print *, i, j, k, data1(i, j, k), u1b(i, j, k)
                  call decomp_2d_abort(4, "x_plane")
               end if
            end do
         end do
      end do

      do k = ystart(3), yend(3)
         do j = ystart(2), yend(2)
            do i = ystart(1), yend(1)
               if (abs(data1(i, j, k) - u2b(i, j, k)) > eps) then
                  print *, i, j, k, data1(i, j, k), u2b(i, j, k)
                  call decomp_2d_abort(5, "y_plane")
               end if
            end do
         end do
      end do

      do k = zstart(3), zend(3)
         do j = zstart(2), zend(2)
            do i = zstart(1), zend(1)
               if (abs(data1(i, j, k) - u3b(i, j, k)) > eps) then
                  print *, i, j, k, data1(i, j, k), u3b(i, j, k)
                  call decomp_2d_abort(6, "z_plane")
               end if
            end do
         end do
      end do

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (nrank == 0) then
         print *, "Checking "//stage//" pass!"
      end if

   end subroutine check

end program io_test
