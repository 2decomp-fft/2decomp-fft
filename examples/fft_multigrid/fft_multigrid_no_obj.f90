!! SPDX-License-Identifier: BSD-3-Clause

module reference_data

   public

contains

   ! This subroutines provides the reference real data
   pure function reference_rdata(i, j, k, nx, ny, nz)

      !$acc routine seq
      use decomp_2d_constants, only: mytype

      implicit none

      integer, intent(in) :: i, j, k, nx, ny, nz
      real(mytype) :: reference_rdata

      reference_rdata = real(i, mytype) / real(nx, mytype) &
                        * real(j, mytype) / real(ny, mytype) &
                        * real(k, mytype) / real(nz, mytype)

   end function reference_rdata

   ! This subroutine provides the reference complex data
   pure function reference_cdata(i, j, k, nx, ny, nz)

      !$acc routine seq
      use decomp_2d_constants, only: mytype

      implicit none

      integer, intent(in) :: i, j, k, nx, ny, nz
      complex(mytype) :: reference_cdata

      reference_cdata = cmplx(reference_rdata(i, j, k, nx, ny, nz), &
                              -3._mytype * reference_rdata(i, j, k, nx, ny, nz), &
                              mytype)

   end function reference_cdata

end module reference_data

! Run a test for c2c transform several times
subroutine ntest_c2c(nt)

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_fft
   use decomp_2d_mpi
   use MPI
   use reference_data
#if defined(_GPU)
   use cudafor
   use cufft
   use openacc
#endif

   implicit none

   integer, intent(in) :: nt

   ! Local variables
   integer :: nx, ny, nz, st1, en1, st2, en2, st3, en3
   integer :: i, j, k, itest, ierror, format
   type(decomp_info), pointer :: ph => null()
   complex(mytype), dimension(:, :, :), allocatable :: in, out
   double precision :: t1, t2, tt
   real(mytype) :: error

   ! Get the decomp_info object describing the arrays
   ph => decomp_2d_fft_get_ph()
   nx = ph%xsz(1); ny = ph%ysz(2); nz = ph%zsz(3)

   ! Get the format
   format = decomp_2d_fft_get_format()

   ! Allocate memory
   if (format == PHYSICAL_IN_X) then
      call alloc_x(in, ph, opt_global=.true.)
      call alloc_z(out, ph, opt_global=.true.)
   else
      call alloc_z(in, ph, opt_global=.true.)
      call alloc_x(out, ph, opt_global=.true.)
   end if

   ! Initialise input
   if (format == PHYSICAL_IN_X) then
      st1 = ph%xst(1); en1 = ph%xen(1)
      st2 = ph%xst(2); en2 = ph%xen(2)
      st3 = ph%xst(3); en3 = ph%xen(3)
   else
      st1 = ph%zst(1); en1 = ph%zen(1)
      st2 = ph%zst(2); en2 = ph%zen(2)
      st3 = ph%zst(3); en3 = ph%zen(3)
   end if
   do k = st3, en3
      do j = st2, en2
         do i = st1, en1
            in(i, j, k) = reference_cdata(i, j, k, nx, ny, nz)
         end do
      end do
   end do

   !$acc data copyin(in), copyout(out)

   ! First iteration is out of the loop with a dedicated timer
   t1 = MPI_WTIME()
   call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
   t1 = MPI_WTIME() - t1
   t2 = MPI_WTIME()
   call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
   t2 = MPI_WTIME() - t2
   ! Rescale input
   !$acc kernels
   in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
   !$acc end kernels
   ! MPI reduction of timers
   call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   if (nrank == 0) then
      write (*, *) '===== c2c interface ====='
      write (*, *) 'Grid ', nx, ny, nz
      if (format == PHYSICAL_IN_X) then
         write (*, *) 'Format PHYSICAL_IN_X'
      else
         write (*, *) 'Format PHYSICAL_IN_Z'
      end if
      write (*, *) ''
      write (*, *) 'First iteration with dedicated timer'
      write (*, *) '     time (sec): ', t1 / dble(nproc), t2 / dble(nproc)
      write (*, *) ''
   end if

   ! Init the timer and run the tests
   t1 = 0.d0
   t2 = 0.d0
   do itest = 1, nt

      ! Forward FFT
      tt = MPI_WTIME()
      call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
      tt = MPI_WTIME() - tt
      t1 = t1 + tt

      ! Inverse FFT
      tt = MPI_WTIME()
      call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
      tt = MPI_WTIME() - tt
      t2 = t2 + tt

      ! Rescale input
      !$acc kernels
      in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
      !$acc end kernels

   end do

#if defined(_GPU)
   ierror = cudaDeviceSynchronize()
#endif

   ! MPI reduction of timers
   call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   t1 = t1 / dble(nproc) / dble(nt)
   t2 = t2 / dble(nproc) / dble(nt)

   ! Compute the sqrt of the L2 error
   error = 0._mytype
   if (format == PHYSICAL_IN_X) then
      st1 = ph%xst(1); en1 = ph%xen(1)
      st2 = ph%xst(2); en2 = ph%xen(2)
      st3 = ph%xst(3); en3 = ph%xen(3)
   else
      st1 = ph%zst(1); en1 = ph%zen(1)
      st2 = ph%zst(2); en2 = ph%zen(2)
      st3 = ph%zst(3); en3 = ph%zen(3)
   end if
   !$acc parallel loop default(present) reduction(+:error)
   do k = st3, en3
      do j = st2, en2
         do i = st1, en1
            error = error + (real(in(i, j, k), mytype) - reference_rdata(i, j, k, nx, ny, nz))**2 &
                    + (aimag(in(i, j, k)) - reference_rdata(i, j, k, nx, ny, nz))**2
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
   error = sqrt(error) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

   ! Print some stuff
   if (nrank == 0) then
      write (*, *) 'error / mesh point: ', error
      write (*, *) 'Avg time (sec): ', t1, t2
      write (*, *) ''
   end if

   !$acc end data

   ! Free memory, nullify pointers
   deallocate (in, out)
   nullify (ph)

end subroutine ntest_c2c

! Run a test for r2c transform several times
subroutine ntest_r2c(nt)

   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_constants
   use decomp_2d_mpi
   use MPI
   use reference_data
#if defined(_GPU)
   use cudafor
   use cufft
   use openacc
#endif

   implicit none

   integer, intent(in) :: nt

   ! Local variables
   integer :: nx, ny, nz, st1, en1, st2, en2, st3, en3
   integer :: i, j, k, itest, ierror, format
   type(decomp_info), pointer :: ph => null(), sp => null()
   real(mytype), dimension(:, :, :), allocatable :: in
   complex(mytype), dimension(:, :, :), allocatable :: out
   double precision :: t1, t2, tt
   real(mytype) :: error

   ! Get the decomp_info object describing the arrays
   ph => decomp_2d_fft_get_ph()
   sp => decomp_2d_fft_get_sp()
   nx = ph%xsz(1); ny = ph%ysz(2); nz = ph%zsz(3)

   ! Get the format
   format = decomp_2d_fft_get_format()

   ! Allocate memory
   if (format == PHYSICAL_IN_X) then
      call alloc_x(in, ph, opt_global=.true.)
      call alloc_z(out, sp, opt_global=.true.)
   else
      call alloc_z(in, ph, opt_global=.true.)
      call alloc_x(out, sp, opt_global=.true.)
   end if

   ! Initialise input
   if (format == PHYSICAL_IN_X) then
      st1 = ph%xst(1); en1 = ph%xen(1)
      st2 = ph%xst(2); en2 = ph%xen(2)
      st3 = ph%xst(3); en3 = ph%xen(3)
   else
      st1 = ph%zst(1); en1 = ph%zen(1)
      st2 = ph%zst(2); en2 = ph%zen(2)
      st3 = ph%zst(3); en3 = ph%zen(3)
   end if
   do k = st3, en3
      do j = st2, en2
         do i = st1, en1
            in(i, j, k) = reference_rdata(i, j, k, nx, ny, nz)
         end do
      end do
   end do

   !$acc data copyin(in), copyout(out)

   ! First iteration is out of the loop with a dedicated timer
   t1 = MPI_WTIME()
   call decomp_2d_fft_3d(in, out)
   t1 = MPI_WTIME() - t1
   t2 = MPI_WTIME()
   call decomp_2d_fft_3d(out, in)
   t2 = MPI_WTIME() - t2
   ! Rescale input
   !$acc kernels
   in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
   !$acc end kernels
   ! MPI reduction of timers
   call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   if (nrank == 0) then
      write (*, *) '===== r2c interface ====='
      write (*, *) 'Grid ', nx, ny, nz
      if (format == PHYSICAL_IN_X) then
         write (*, *) 'Format PHYSICAL_IN_X'
      else
         write (*, *) 'Format PHYSICAL_IN_Z'
      end if
      write (*, *) ''
      write (*, *) 'First iteration with dedicated timer'
      write (*, *) '     time (sec): ', t1 / dble(nproc), t2 / dble(nproc)
      write (*, *) ''
   end if

   ! Init the timer and run the tests
   t1 = 0.d0
   t2 = 0.d0
   do itest = 1, nt

      ! Forward FFT
      tt = MPI_WTIME()
      call decomp_2d_fft_3d(in, out)
      tt = MPI_WTIME() - tt
      t1 = t1 + tt

      ! Inverse FFT
      tt = MPI_WTIME()
      call decomp_2d_fft_3d(out, in)
      tt = MPI_WTIME() - tt
      t2 = t2 + tt

      ! Rescale input
      !$acc kernels
      in = in / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
      !$acc end kernels

   end do

#if defined(_GPU)
   ierror = cudaDeviceSynchronize()
#endif

   ! MPI reduction of timers
   call MPI_ALLREDUCE(MPI_IN_PLACE, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   call MPI_ALLREDUCE(MPI_IN_PLACE, t2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
   t1 = t1 / dble(nproc) / dble(nt)
   t2 = t2 / dble(nproc) / dble(nt)

   ! Compute the sqrt of the L2 error
   error = 0._mytype
   if (format == PHYSICAL_IN_X) then
      st1 = ph%xst(1); en1 = ph%xen(1)
      st2 = ph%xst(2); en2 = ph%xen(2)
      st3 = ph%xst(3); en3 = ph%xen(3)
   else
      st1 = ph%zst(1); en1 = ph%zen(1)
      st2 = ph%zst(2); en2 = ph%zen(2)
      st3 = ph%zst(3); en3 = ph%zen(3)
   end if
   !$acc parallel loop default(present) reduction(+:error)
   do k = st3, en3
      do j = st2, en2
         do i = st1, en1
            error = error + (in(i, j, k) - reference_rdata(i, j, k, nx, ny, nz))**2
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
   error = sqrt(error) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

   ! Print some stuff
   if (nrank == 0) then
      write (*, *) 'error / mesh point: ', error
      write (*, *) 'Avg time (sec): ', t1, t2
      write (*, *) ''
   end if

   !$acc end data

   ! Free memory, nullify pointers
   deallocate (in, out)
   nullify (sp)
   nullify (ph)

end subroutine ntest_r2c

program fft_multigrid_no_obj

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_fft
   use decomp_2d_mpi
   use decomp_2d_testing
   use MPI
#if defined(_GPU)
   use cudafor
   use cufft
   use openacc
#endif

   implicit none

   interface
      subroutine ntest_c2c(ntest)
         implicit none
         integer, intent(in) :: ntest
      end subroutine ntest_c2c

      subroutine ntest_r2c(ntest)
         implicit none
         integer, intent(in) :: ntest
      end subroutine ntest_r2c
   end interface

   integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot

   integer :: ntest = 10  ! repeat test this times
   integer :: ngrid = 4
   integer :: ierror, igrid

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   ! Now we can check if user put some inputs
   call decomp_2d_testing_init(p_row, p_col, nx, ny, nz, ntest)

   call decomp_2d_init(nx + 1, ny + 1, nz + 1, p_row, p_col)

   call decomp_2d_testing_log()

   ! Set the number of FFT grids
   call decomp_2d_fft_set_ngrid(ngrid)

   ! Init each FFT grid
   call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, 1)
   if (ngrid > 1) call decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz, 2)
   if (ngrid > 2) call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny + 2, nz + 16, 3)
   if (ngrid > 3) call decomp_2d_fft_init(PHYSICAL_IN_Z, nx + 16, ny + 2, nz, 4)

   ! Test each grid
   do igrid = 1, ngrid
      ! Use a specific grid
      call decomp_2d_fft_use_grid(igrid)
      ! Test the grid
      call ntest_c2c(ntest)
      call ntest_r2c(ntest)
   end do

   call decomp_2d_fft_finalize
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program fft_multigrid_no_obj
