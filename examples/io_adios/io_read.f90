program io_read

   use mpi
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io
   use decomp_2d_io_family
   use decomp_2d_io_object
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

   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#else
   real(mytype), allocatable, dimension(:, :, :) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   character(len=*), parameter :: io_name = "test-io"
   type(d2d_io_family), save :: io_family
   type(d2d_io), save :: io
#ifndef ADIOS2
   logical ::file_exists1, file_exists2, file_exists3
#endif

   integer :: i, j, k, m, ierror

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
   call decomp_2d_io_register_var3d("u1.dat", 1, real_type)
   call io_family%init(io_name)
   call io_family%register_var3d("u2.dat", 2, real_type)
   call io_family%register_var3d("u3.dat", 3, real_type)

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
   call io%open_start("out", decomp_2d_read_mode)
#endif
   call decomp_2d_read_one(1, u1b, 'u1.dat', opt_dirname='out', opt_io=io)
#ifdef ADIOS2
   call io%end_close
#endif
#ifdef ADIOS2
   call io%open_start("out", decomp_2d_read_mode, opt_family=io_family)
#endif
   call decomp_2d_read_one(2, u2b, 'u2.dat', opt_dirname='out', opt_io=io)
   call decomp_2d_read_one(3, u3b, 'u3.dat', opt_dirname='out', opt_io=io)
#ifdef ADIOS2
   call io%end_close
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
   call io_family%fin
   call decomp_2d_io_fin
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_read
