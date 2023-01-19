program io_read

   use mpi

   use decomp_2d
   use decomp_2d_io

   implicit none

   integer, parameter :: nx = 17, ny = 13, nz = 11
   ! use different number of processes
   integer :: p_row = 0, p_col = 0

#ifdef COMPLEX_TEST
   complex(mytype), dimension(nx, ny, nz) :: data1

   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#else
   real(mytype), dimension(nx, ny, nz) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   integer :: i, j, k, m, ierror

   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)
   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   ! ***** global data *****
   m = 1
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
#ifdef COMPLEX_TEST
            data1(i, j, k) = cmplx(real(m, mytype), real(nx*ny*nz - m, mytype))
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
   call decomp_2d_read_one(1, u1b, '.', 'u1.dat', 'test', reduce_prec=.false.)
   call decomp_2d_read_one(2, u2b, '.', 'u2.dat', 'test', reduce_prec=.false.)
   call decomp_2d_read_one(3, u3b, '.', 'u3.dat', 'test', reduce_prec=.false.)

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
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_read
