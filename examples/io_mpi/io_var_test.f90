!! SPDX-License-Identifier: BSD-3-Clause
! Sample application to test the read/write_var sets of routines
! in the IO library

program io_var_test

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

   real(mytype), parameter :: eps = 1.0E-7

   ! for global data
   real(mytype), allocatable, dimension(:, :, :) :: data1
   real(mytype), allocatable, dimension(:, :, :) :: data1_large
   complex(mytype), allocatable, dimension(:, :, :) :: cdata1

   ! for distributed data
   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
   real(mytype), allocatable, dimension(:, :, :) :: u1l, u2l, u3l
   complex(mytype), allocatable, dimension(:, :, :) :: cu1, cu2, cu3

   ! another copy
   real(mytype), allocatable, dimension(:, :, :) :: u1_b, u2_b, u3_b
   real(mytype), allocatable, dimension(:, :, :) :: u1l_b, u2l_b, u3l_b
   complex(mytype), allocatable, dimension(:, :, :) :: cu1_b, cu2_b, cu3_b

   real(mytype), allocatable, dimension(:) :: tmp
   complex(mytype), allocatable, dimension(:) :: ctmp
   integer, allocatable, dimension(:) :: itmp

   type(d2d_io_mpi) :: io

   integer :: xst1, xst2, xst3
   integer :: xen1, xen2, xen3

   TYPE(DECOMP_INFO) :: large

   integer :: i, j, k, m, ierror
   character(len=15) :: filename

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

   call decomp_2d_init(nx, ny, nz, p_row, p_col, complex_pool=.true.)

   call decomp_2d_testing_log()

   call decomp_2d_io_init()

   ! also create a data set over a large domain
   call decomp_info_init(nx * 2, ny * 2, nz * 2, large)

   ! initialise global data
   allocate (data1(nx, ny, nz))
   allocate (cdata1(nx, ny, nz))
   allocate (data1_large(nx * 2, ny * 2, nz * 2))
   m = 1
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            data1(i, j, k) = real(m, mytype)
            cdata1(i, j, k) = cmplx(real(m, mytype), real(m, mytype), kind=mytype)
            m = m + 1
         end do
      end do
   end do

   m = 1
   do k = 1, nz * 2
      do j = 1, ny * 2
         do i = 1, nx * 2
            data1_large(i, j, k) = real(m, mytype)
            m = m + 1
         end do
      end do
   end do

   ! allocate memory
   call alloc_x(u1, .true.)
   call alloc_y(u2, .true.)
   call alloc_z(u3, .true.)
   call alloc_x(u1_b, .true.)
   call alloc_y(u2_b, .true.)
   call alloc_z(u3_b, .true.)

   call alloc_x(cu1, .true.)
   call alloc_y(cu2, .true.)
   call alloc_z(cu3, .true.)
   call alloc_x(cu1_b, .true.)
   call alloc_y(cu2_b, .true.)
   call alloc_z(cu3_b, .true.)

   call alloc_x(u1l, large, .true.)
   call alloc_y(u2l, large, .true.)
   call alloc_z(u3l, large, .true.)
   call alloc_x(u1l_b, large, .true.)
   call alloc_y(u2l_b, large, .true.)
   call alloc_z(u3l_b, large, .true.)

   ! distribute the data
   !$acc data copyin(data1,cdata1,data1_large) copy(u1,u2,u3,u1l,u2l,u3l,cu1,cu2,cu3)
   xst1 = xstart(1); xen1 = xend(1)
   xst2 = xstart(2); xen2 = xend(2)
   xst3 = xstart(3); xen3 = xend(3)
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1(i, j, k) = data1(i, j, k)
            cu1(i, j, k) = cdata1(i, j, k)
         end do
      end do
   end do
   !$acc end loop
   xst1 = large%xst(1); xen1 = large%xen(1)
   xst2 = large%xst(2); xen2 = large%xen(2)
   xst3 = large%xst(3); xen3 = large%xen(3)
   !$acc parallel loop default(present)
   do k = xst3, xen3
      do j = xst2, xen2
         do i = xst1, xen1
            u1l(i, j, k) = data1_large(i, j, k)
         end do
      end do
   end do
   !$acc end loop

   ! transpose
   call transpose_x_to_y(u1, u2)
   call transpose_y_to_z(u2, u3)
   call transpose_x_to_y(u1l, u2l, large)
   call transpose_y_to_z(u2l, u3l, large)
   call transpose_x_to_y(cu1, cu2)
   call transpose_y_to_z(cu2, cu3)
   !$acc update self (u1)
   !$acc update self (u2)
   !$acc update self (u3)
   !$acc update self (u1l)
   !$acc update self (u2l)
   !$acc update self (u3l)
   !$acc update self (cu1)
   !$acc update self (cu2)
   !$acc update self (cu3)
   !$acc end data

   ! open file for IO
   write (filename, '(A,I3.3)') 'io_var_data.', nproc
   call io%open(filename, decomp_2d_write_mode)

   ! test writing scalar data
   allocate (tmp(2))
   tmp(1) = 1._mytype
   tmp(2) = 2._mytype
   allocate (ctmp(3))
   ctmp(1) = cmplx(1.0, 1.0, mytype)
   ctmp(2) = cmplx(2.0, 2.0, mytype)
   ctmp(3) = cmplx(3.0, 3.0, mytype)
   allocate (itmp(3))
   call decomp_2d_write_scalar(io, 2, tmp)
   call decomp_2d_write_scalar(io, 3, ctmp)
   call decomp_2d_write_scalar(io, 3, (/nx, ny, nz/))

   ! test the IO routines by writing all data to disk
   call decomp_2d_write_var(io, 1, u1)
   call decomp_2d_write_var(io, 2, u2)
   call decomp_2d_write_var(io, 3, u3)
   call decomp_2d_write_var(io, 1, u1l, opt_decomp=large)
   call decomp_2d_write_var(io, 2, u2l, opt_decomp=large)
   call decomp_2d_write_var(io, 3, u3l, opt_decomp=large)
   call decomp_2d_write_var(io, 1, cu1)
   call decomp_2d_write_var(io, 2, cu2)
   call decomp_2d_write_var(io, 3, cu3)

   if (nrank == 0) write (*, *) 'disp=', io%disp

   call io%close

   ! read data back in from file
   call io%open(filename, decomp_2d_read_mode)
   ! skip the scalars (2 real, 3 cmplx, 3 int)
#ifdef DOUBLE_PREC
   ! if double precision: 2*8+3*8*2+3*4
   io%disp = 76_MPI_OFFSET_KIND
#else
   ! if single precision: 2*4+3*4*2+3*4
   io%disp = 44_MPI_OFFSET_KIND
#endif

   call decomp_2d_read_var(io, 1, u1_b)
   call decomp_2d_read_var(io, 2, u2_b)
   call decomp_2d_read_var(io, 3, u3_b)
   call decomp_2d_read_var(io, 1, u1l_b, opt_decomp=large)
   call decomp_2d_read_var(io, 2, u2l_b, opt_decomp=large)
   call decomp_2d_read_var(io, 3, u3l_b, opt_decomp=large)
   call decomp_2d_read_var(io, 1, cu1_b)
   call decomp_2d_read_var(io, 2, cu2_b)
   call decomp_2d_read_var(io, 3, cu3_b)

   io%disp = 0_MPI_OFFSET_KIND
   call decomp_2d_read_scalar(io, 2, tmp)
   call decomp_2d_read_scalar(io, 3, ctmp)
   call decomp_2d_read_scalar(io, 3, itmp)
   if (nrank == 0) then
      write (*, '(2F8.3)') tmp
      write (*, 20) ctmp
20    format(3(:, '(', F5.2, ',', F5.2, ')'))
      write (*, '(A,3I5)') 'nx,ny,nz', itmp
   end if

   call io%close
   deallocate (tmp, ctmp, itmp)

   ! validate the data
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            if (abs(u1(i, j, k) - u1_b(i, j, k)) > eps) stop 1
            if (abs(cu1(i, j, k) - cu1_b(i, j, k)) > eps) stop 2
         end do
      end do
   end do

   do k = ystart(3), yend(3)
      do j = ystart(2), yend(2)
         do i = ystart(1), yend(1)
            if (abs(u2(i, j, k) - u2_b(i, j, k)) > eps) stop 3
            if (abs(cu2(i, j, k) - cu2_b(i, j, k)) > eps) stop 4
         end do
      end do
   end do

   do k = zstart(3), zend(3)
      do j = zstart(2), zend(2)
         do i = zstart(1), zend(1)
            if (abs(u3(i, j, k) - u3_b(i, j, k)) > eps) stop 5
            if (abs(cu3(i, j, k) - cu3_b(i, j, k)) > eps) stop 6
         end do
      end do
   end do

   do k = large%xst(3), large%xen(3)
      do j = large%xst(2), large%xen(2)
         do i = large%xst(1), large%xen(1)
            if (abs(u1l(i, j, k) - u1l_b(i, j, k)) > eps) stop 7
         end do
      end do
   end do

   do k = large%yst(3), large%yen(3)
      do j = large%yst(2), large%yen(2)
         do i = large%yst(1), large%yen(1)
            if (abs(u2l(i, j, k) - u2l_b(i, j, k)) > eps) stop 8
         end do
      end do
   end do

   do k = large%zst(3), large%zen(3)
      do j = large%zst(2), large%zen(2)
         do i = large%zst(1), large%zen(1)
            if (abs(u3l(i, j, k) - u3l_b(i, j, k)) > eps) stop 9
         end do
      end do
   end do

   if (nrank == 0) write (*, *) 'passed self test'

   ! clean up
   deallocate (u1, u2, u3)
   deallocate (u1l, u2l, u3l)
   deallocate (cu1, cu2, cu3)
   deallocate (u1_b, u2_b, u3_b)
   deallocate (u1l_b, u2l_b, u3l_b)
   deallocate (cu1_b, cu2_b, cu3_b)
   deallocate (data1, cdata1, data1_large)
   call decomp_info_finalize(large)
   call decomp_2d_io_fin
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program io_var_test
