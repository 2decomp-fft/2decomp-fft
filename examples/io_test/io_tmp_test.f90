!!! io_tmp_test.f90
!!
!! Tests writing from temporary arrays using the flush functionality, only relevant for ADIOS2.
!!
!! SPDX-License-Identifier: BSD-3-Clause
program io_test

   use mpi
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d
   use decomp_2d_io
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
   integer :: nargin, arg, FNLength, status, DecInd
   character(len=80) :: InputFN

   integer :: ierr

#ifdef COMPLEX_TEST
   complex(mytype), allocatable, dimension(:, :, :) :: data1

   complex(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
   complex(mytype), allocatable, dimension(:, :, :) :: v1, v2, v3
   complex(mytype), allocatable, dimension(:, :, :) :: v1b, v2b, v3b

   complex(mytype), allocatable, dimension(:, :, :) :: tmp1, tmp2, tmp3
#else
   real(mytype), allocatable, dimension(:, :, :) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
   real(mytype), allocatable, dimension(:, :, :) :: v1, v2, v3
   real(mytype), allocatable, dimension(:, :, :) :: v1b, v2b, v3b

   real(mytype), allocatable, dimension(:, :, :) :: tmp1, tmp2, tmp3
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   character(len=*), parameter :: io_name = "test-io"
   character(len=*), parameter :: io_restart = "restart-io"

   integer :: i, j, k, m, ierror
   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3

#ifndef ADIOS2
   logical :: dir_exists
#endif

   integer, parameter :: output2D = 0 ! Which plane to write in 2D (0 for 3D)

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

   call decomp_2d_io_init()
   call decomp_2d_init_io(io_name)
   call decomp_2d_register_variable(io_name, "u1.dat", 1, 0, output2D, mytype)
   call decomp_2d_register_variable(io_name, "u2.dat", 2, 0, output2D, mytype)
   call decomp_2d_register_variable(io_name, "u3.dat", 3, 0, output2D, mytype)
   call decomp_2d_register_variable(io_name, "v1.dat", 1, 0, output2D, mytype)
   call decomp_2d_register_variable(io_name, "v2.dat", 2, 0, output2D, mytype)
   call decomp_2d_register_variable(io_name, "v3.dat", 3, 0, output2D, mytype)
   call decomp_2d_init_io(io_restart)
   call decomp_2d_register_variable(io_restart, "u1.dat", 1, 0, output2D, mytype)
   call decomp_2d_register_variable(io_restart, "u2.dat", 2, 0, output2D, mytype)
   call decomp_2d_register_variable(io_restart, "u3.dat", 3, 0, output2D, mytype)
   call decomp_2d_register_variable(io_restart, "v1.dat", 1, 0, output2D, mytype)
   call decomp_2d_register_variable(io_restart, "v2.dat", 2, 0, output2D, mytype)
   call decomp_2d_register_variable(io_restart, "v3.dat", 3, 0, output2D, mytype)

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

   call alloc_x(u1b, .true.)
   call alloc_y(u2b, .true.)
   call alloc_z(u3b, .true.)

   call alloc_x(v1, .true.)
   call alloc_y(v2, .true.)
   call alloc_z(v3, .true.)

   call alloc_x(v1b, .true.)
   call alloc_y(v2b, .true.)
   call alloc_z(v3b, .true.)

   call alloc_x(tmp1, .true.)
   call alloc_y(tmp2, .true.)
   call alloc_z(tmp3, .true.)

   xst1 = xstart(1); xen1 = xend(1)
   xst2 = xstart(2); xen2 = xend(2)
   xst3 = xstart(3); xen3 = xend(3)
   ! original x-pencil based data
   !$acc data copyin(data1) copy(u1,u2,u3,v1,v2,v3)
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1(i, j, k) = data1(i, j, k)
            v1(i, j, k) = 10 * data1(i, j, k)
         end do
      end do
   end do
   !$acc end loop

   ! transpose
   call transpose_x_to_y(u1, u2)
   call transpose_y_to_z(u2, u3)
   call transpose_x_to_y(v1, v2)
   call transpose_y_to_z(v2, v3)
   !$acc update self(u1)
   !$acc update self(u2)
   !$acc update self(u3)
   !$acc update self(v1)
   !$acc update self(v2)
   !$acc update self(v3)
   !$acc end data

   ! write to disk
#ifndef ADIOS2
   if (nrank == 0) then
      inquire (file="out", exist=dir_exists)
      if (.not. dir_exists) then
         call execute_command_line("mkdir out 2> /dev/null")
      end if
   end if
#endif

   ! Standard I/O pattern - file per field
#ifdef ADIOS2
   call decomp_2d_open_io(io_name, "out", decomp_2d_write_mode)
   call decomp_2d_start_io(io_name, "out")
#endif
   ! Copy data to temporary memory and write from the temporary memory

   tmp1(:, :, :) = u1(:, :, :)
   call decomp_2d_write_one(1, tmp1, 'out', 'u1.dat', 0, io_name, opt_deferred_writes=.false.)
   tmp1(:, :, :) = v1(:, :, :)
   call decomp_2d_write_one(1, tmp1, 'out', 'v1.dat', 0, io_name, opt_deferred_writes=.false.)

   tmp2(:, :, :) = u2(:, :, :)
   call decomp_2d_write_one(2, tmp2, 'out', 'u2.dat', 0, io_name, opt_deferred_writes=.false.)
   tmp2(:, :, :) = v2(:, :, :)
   call decomp_2d_write_one(2, tmp2, 'out', 'v2.dat', 0, io_name, opt_deferred_writes=.false.)

   tmp3(:, :, :) = u3(:, :, :)
   call decomp_2d_write_one(3, tmp3, 'out', 'u3.dat', 0, io_name, opt_deferred_writes=.false.)
   tmp3(:, :, :) = v3(:, :, :)
   call decomp_2d_write_one(3, tmp3, 'out', 'v3.dat', 0, io_name, opt_deferred_writes=.false.)
#ifdef ADIOS2
   call decomp_2d_end_io(io_name, "out")
   call decomp_2d_close_io(io_name, "out")
#endif

   ! read back to different arrays
#ifdef ADIOS2
   call decomp_2d_open_io(io_name, "out", decomp_2d_read_mode)
   call decomp_2d_start_io(io_name, "out")
#endif
   call decomp_2d_read_one(1, u1b, 'out', 'u1.dat', io_name, reduce_prec=.false.)
   call decomp_2d_read_one(2, u2b, 'out', 'u2.dat', io_name, reduce_prec=.false.)
   call decomp_2d_read_one(3, u3b, 'out', 'u3.dat', io_name, reduce_prec=.false.)

   call decomp_2d_read_one(1, v1b, 'out', 'v1.dat', io_name, reduce_prec=.false.)
   call decomp_2d_read_one(2, v2b, 'out', 'v2.dat', io_name, reduce_prec=.false.)
   call decomp_2d_read_one(3, v3b, 'out', 'v3.dat', io_name, reduce_prec=.false.)
#ifdef ADIOS2
   call decomp_2d_end_io(io_name, "out")
   call decomp_2d_close_io(io_name, "out")
#endif

   ! compare
   call check("file per field")

   ! Checkpoint I/O pattern - multiple fields per file
   call decomp_2d_open_io(io_name, "checkpoint", decomp_2d_write_mode)
   call decomp_2d_start_io(io_name, "checkpoint")

   tmp1(:, :, :) = u1(:, :, :)
   call decomp_2d_write_one(1, tmp1, 'checkpoint', 'u1.dat', 0, io_name, opt_deferred_writes=.false.)
   tmp1(:, :, :) = v1(:, :, :)
   call decomp_2d_write_one(1, tmp1, 'checkpoint', 'v1.dat', 0, io_name, opt_deferred_writes=.false.)

   tmp2(:, :, :) = u2(:, :, :)
   call decomp_2d_write_one(2, tmp2, 'checkpoint', 'u2.dat', 0, io_name, opt_deferred_writes=.false.)
   tmp2(:, :, :) = v2(:, :, :)
   call decomp_2d_write_one(2, tmp2, 'checkpoint', 'v2.dat', 0, io_name, opt_deferred_writes=.false.)

   tmp3(:, :, :) = u3(:, :, :)
   call decomp_2d_write_one(3, tmp3, 'checkpoint', 'u3.dat', 0, io_name, opt_deferred_writes=.false.)
   tmp3(:, :, :) = v3(:, :, :)
   call decomp_2d_write_one(3, tmp3, 'checkpoint', 'v3.dat', 0, io_name, opt_deferred_writes=.false.)

   call decomp_2d_end_io(io_name, "checkpoint")
   call decomp_2d_close_io(io_name, "checkpoint")

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   ! read back to different arrays
   ! XXX: For the MPI-IO backend the order of reading must match the order of writing!
   u1b = 0; u2b = 0; u3b = 0
   v1b = 0; v2b = 0; v3b = 0
   call decomp_2d_open_io(io_restart, "checkpoint", decomp_2d_read_mode)
   call decomp_2d_start_io(io_restart, "checkpoint")

   call decomp_2d_read_one(1, u1b, 'checkpoint', 'u1.dat', io_restart, reduce_prec=.false.)
   call decomp_2d_read_one(1, v1b, 'checkpoint', 'v1.dat', io_restart, reduce_prec=.false.)

   call decomp_2d_read_one(2, u2b, 'checkpoint', 'u2.dat', io_restart, reduce_prec=.false.)
   call decomp_2d_read_one(2, v2b, 'checkpoint', 'v2.dat', io_restart, reduce_prec=.false.)

   call decomp_2d_read_one(3, u3b, 'checkpoint', 'u3.dat', io_restart, reduce_prec=.false.)
   call decomp_2d_read_one(3, v3b, 'checkpoint', 'v3.dat', io_restart, reduce_prec=.false.)

   call decomp_2d_end_io(io_restart, "checkpoint")
   call decomp_2d_close_io(io_restart, "checkpoint")

   call MPI_Barrier(MPI_COMM_WORLD, ierr)

   ! compare
   call check("one file, multiple fields")

   deallocate (u1, u2, u3)
   deallocate (u1b, u2b, u3b)
   deallocate (v1, v2, v3)
   deallocate (v1b, v2b, v3b)
   deallocate (data1)
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
                  print *, u1(i, j, k), u1b(i, j, k)
                  stop 1
               end if
               if (abs((v1(i, j, k) - v1b(i, j, k))) > eps) then
                  print *, v1(i, j, k), v1b(i, j, k)
                  stop 11
               end if
            end do
         end do
      end do

      do k = ystart(3), yend(3)
         do j = ystart(2), yend(2)
            do i = ystart(1), yend(1)
               if (abs((u2(i, j, k) - u2b(i, j, k))) > eps) then
                  print *, u2(i, j, k), u2b(i, j, k)
                  stop 2
               end if
               if (abs((v2(i, j, k) - v2b(i, j, k))) > eps) stop 22
            end do
         end do
      end do

      do k = zstart(3), zend(3)
         do j = zstart(2), zend(2)
            do i = zstart(1), zend(1)
               if (abs((u3(i, j, k) - u3b(i, j, k))) > eps) stop 3
               if (abs((v3(i, j, k) - v3b(i, j, k))) > eps) stop 33
            end do
         end do
      end do

      ! Also check against the global data array
      do k = xstart(3), xend(3)
         do j = xstart(2), xend(2)
            do i = xstart(1), xend(1)
               if (abs(data1(i, j, k) - u1b(i, j, k)) > eps) stop 4
               if (abs(10 * data1(i, j, k) - v1b(i, j, k)) > eps) stop 44
            end do
         end do
      end do

      do k = ystart(3), yend(3)
         do j = ystart(2), yend(2)
            do i = ystart(1), yend(1)
               if (abs((data1(i, j, k) - u2b(i, j, k))) > eps) stop 5
               if (abs((10 * data1(i, j, k) - v2b(i, j, k))) > eps) stop 55
            end do
         end do
      end do

      do k = zstart(3), zend(3)
         do j = zstart(2), zend(2)
            do i = zstart(1), zend(1)
               if (abs((data1(i, j, k) - u3b(i, j, k))) > eps) stop 6
               if (abs((10 * data1(i, j, k) - v3b(i, j, k))) > eps) stop 66
            end do
         end do
      end do

      call MPI_Barrier(MPI_COMM_WORLD, ierr)
      if (nrank == 0) then
         print *, "Checking "//stage//" pass!"
      end if

   end subroutine check

end program io_test
