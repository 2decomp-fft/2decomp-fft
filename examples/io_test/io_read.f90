program io_read

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

#ifdef COMPLEX_TEST
   complex(mytype), allocatable, dimension(:, :, :) :: data1

   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#else
   real(mytype), allocatable, dimension(:, :, :) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   character(len=*), parameter :: io_name = "test-io"
#ifndef ADIOS2
   logical ::file_exists1, file_exists2, file_exists3
#endif

   integer :: i, j, k, m, ierror

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

#ifndef ADIOS2
   if (nrank == 0) then
      inquire (file="out/u1.dat", exist=file_exists1)
      inquire (file="out/u2.dat", exist=file_exists2)
      inquire (file="out/u3.dat", exist=file_exists3)
      if (.not. (file_exists1 .and. file_exists2 .and. file_exists3)) then
         call decomp_2d_abort(1, "Error, data 'out/u<1,2,3>.dat' must exist before running io_read test case!")
      end if
   end if
#endif

   call decomp_2d_io_init()
   call decomp_2d_init_io(io_name)
   call decomp_2d_register_variable(io_name, "u1.dat", 1, 0, output2D, mytype)
   call decomp_2d_register_variable(io_name, "u2.dat", 2, 0, output2D, mytype)
   call decomp_2d_register_variable(io_name, "u3.dat", 3, 0, output2D, mytype)

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

   call alloc_x(u1b, .true.)
   call alloc_y(u2b, .true.)
   call alloc_z(u3b, .true.)

   ! read back to different arrays
#ifdef ADIOS2
   call decomp_2d_open_io(io_name, "out", decomp_2d_read_mode)
   call decomp_2d_start_io(io_name, "out")
#endif
   call decomp_2d_read_one(1, u1b, 'out', 'u1.dat', io_name, reduce_prec=.false.)
   call decomp_2d_read_one(2, u2b, 'out', 'u2.dat', io_name, reduce_prec=.false.)
   call decomp_2d_read_one(3, u3b, 'out', 'u3.dat', io_name, reduce_prec=.false.)
#ifdef ADIOS2
   call decomp_2d_end_io(io_name, "out")
   call decomp_2d_close_io(io_name, "out")
#endif

   ! Check against the global data array
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            if (abs((data1(i, j, k) - u1b(i, j, k))) > eps) stop 4
         end do
      end do
   end do

   do k = ystart(3), yend(3)
      do j = ystart(2), yend(2)
         do i = ystart(1), yend(1)
            if (abs((data1(i, j, k) - u2b(i, j, k))) > eps) stop 5
         end do
      end do
   end do

   do k = zstart(3), zend(3)
      do j = zstart(2), zend(2)
         do i = zstart(1), zend(1)
            if (abs((data1(i, j, k) - u3b(i, j, k))) > eps) stop 6
         end do
      end do
   end do

   deallocate (u1b, u2b, u3b)
   deallocate (data1)
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_read
