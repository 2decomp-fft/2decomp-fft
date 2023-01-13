program io_test

   use mpi

   use decomp_2d
   use decomp_2d_io

   implicit none

   integer, parameter :: nx = 17, ny = 13, nz = 11
   integer :: p_row = 0, p_col = 0

#ifdef COMPLEX_TEST
   complex(mytype), dimension(nx, ny, nz) :: data1

   complex(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   complex(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#else
   real(mytype), dimension(nx, ny, nz) :: data1

   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   real(mytype), allocatable, dimension(:, :, :) :: u1b, u2b, u3b
#endif

   real(mytype), parameter :: eps = 1.0E-7_mytype

   integer :: i, j, k, m, ierror

   call MPI_INIT(ierror)
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

   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)

   call alloc_x(u1b, .true.)
   call alloc_y(u2b, .true.)
   call alloc_z(u3b, .true.)

   ! original x-pencil based data
   !$acc data copyin(data1,xstart,xend) copy(u1,u2,u3) 
   !$acc parallel loop default(present)
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
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
   call decomp_2d_write_one(1, u1, '.', 'u1.dat', 0, 'test')
   call decomp_2d_write_one(2, u2, '.', 'u2.dat', 0, 'test')
   call decomp_2d_write_one(3, u3, '.', 'u3.dat', 0, 'test')

   ! read back to different arrays
   call decomp_2d_read_one(1, u1b, '.', 'u1.dat', 'test', reduce_prec=.false.)
   call decomp_2d_read_one(2, u2b, '.', 'u2.dat', 'test', reduce_prec=.false.)
   call decomp_2d_read_one(3, u3b, '.', 'u3.dat', 'test', reduce_prec=.false.)

   ! compare
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            if (abs((u2(i, j, k) - u2b(i, j, k))) > eps) stop 1
         end do
      end do
   end do

   do k = ystart(3), yend(3)
      do j = ystart(2), yend(2)
         do i = ystart(1), yend(1)
            if (abs((u2(i, j, k) - u2b(i, j, k))) > eps) stop 2
         end do
      end do
   end do

   do k = zstart(3), zend(3)
      do j = zstart(2), zend(2)
         do i = zstart(1), zend(1)
            if (abs((u3(i, j, k) - u3b(i, j, k))) > eps) stop 3
         end do
      end do
   end do

   ! Also check against the global data array
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            if (abs(data1(i, j, k) - u1b(i, j, k)) > eps) stop 4
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

   deallocate (u1, u2, u3)
   deallocate (u1b, u2b, u3b)
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_test
