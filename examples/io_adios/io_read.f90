program io_read

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

#ifdef COMPLEX_TEST
   complex(mytype), allocatable, dimension(:, :, :) :: data1

   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u2c, u3b
#else
   real(mytype), allocatable, dimension(:, :, :) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u2c, u3b
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   character(len=*), parameter :: io_name = "test-io"
   type(d2d_io_family), save :: io_family
   type(d2d_io_adios), save :: io

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

   call decomp_2d_io_init()
   call io_family%init(io_name)

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
   call alloc_y(u2c, .true.)
   call alloc_z(u3b, .true.)

   ! read back to different arrays
   call decomp_2d_adios_read_var(io, 1, u1b, 'u1.dat')
   call decomp_2d_adios_read_var(io, 2, u2b, 'u2.dat')
   call io%end_close
   call io%open_start(decomp_2d_read_mode, opt_family=io_family)
   call decomp_2d_adios_read_var(io, 2, u2c, 'u2.dat', opt_step_start=0)
   call decomp_2d_adios_read_var(io, 3, u3b, 'u3.dat', opt_step_start=0)
   call io%end_close

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

   do k = ystart(3), yend(3)
      do j = ystart(2), yend(2)
         do i = ystart(1), yend(1)
            if (abs((data1(i, j, k) - u2c(i, j, k))) > eps) stop 5
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

   deallocate (u1b, u2b, u2c, u3b)
   deallocate (data1)
   call io_family%fin
   call decomp_2d_io_fin
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_read
