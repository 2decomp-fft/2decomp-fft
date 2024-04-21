!! SPDX-License-Identifier: BSD-3-Clause
program fft_c2c_z_skip

   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_testing
   use MPI
#if defined(_GPU)
   use cudafor
   use cufft
   use openacc
#endif

   implicit none

   integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot

   integer :: ntest = 10  ! repeat test this times

   type(decomp_info), pointer :: ph => null()
   complex(mytype), allocatable, dimension(:, :, :) :: in, out

   real(mytype) :: dr, di, error
   integer :: ierror, i, j, k, m
   integer :: zst1, zst2, zst3
   integer :: zen1, zen2, zen3
   double precision :: n1, flops, t1, t2, t3, t4
   logical, dimension(3) :: skip_c2c = [.true., .false., .true.]

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

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Test the c2c interface
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call decomp_2d_fft_init(PHYSICAL_IN_Z, opt_skip_XYZ_c2c=skip_c2c) ! non-default Z-pencil input

   ph => decomp_2d_fft_get_ph()
   !  input is Z-pencil data
   ! output is X-pencil data
   call alloc_z(in, ph, .true.)
   call alloc_x(out, ph, .true.)
   zst1 = zstart(1); zen1 = zend(1)
   zst2 = zstart(2); zen2 = zend(2)
   zst3 = zstart(3); zen3 = zend(3)

   ! initilise input
   do k = zst3, zen3
      do j = zst2, zen2
         do i = zst1, zen1
            dr = real(i, mytype) / real(nx, mytype) * real(j, mytype) &
                 / real(ny, mytype) * real(k, mytype) / real(nz, mytype)
            di = dr
            in(i, j, k) = cmplx(dr, di, mytype)
         end do
      end do
   end do

   !$acc data copyin(in) copy(out)
   ! First iterations out of the counting loop
   t1 = MPI_WTIME()
   call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
   t2 = MPI_WTIME() - t1
   t3 = MPI_WTIME()
   call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
   t4 = MPI_WTIME() - t3
   call MPI_ALLREDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t1 = t1 / dble(nproc)
   call MPI_ALLREDUCE(t4, t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t3 = t3 / dble(nproc)
   if (nrank == 0) then
      write (*, *) '===== c2c interface ====='
      write (*, *) 'First iteration with dedicated timer'
      write (*, *) '     time (sec): ', t1, t3
      write (*, *) ''
   end if
   ! Init the time
   t2 = 0.d0
   t4 = 0.d0
   !$acc kernels
   if (.not. skip_c2c(1)) in = in / real(nx, mytype)
   if (.not. skip_c2c(2)) in = in / real(ny, mytype)
   if (.not. skip_c2c(3)) in = in / real(nz, mytype)
   !$acc end kernels
   do m = 1, ntest

      ! forward FFT
      t1 = MPI_WTIME()
      call decomp_2d_fft_3d(in, out, DECOMP_2D_FFT_FORWARD)
      t2 = t2 + MPI_WTIME() - t1

      ! inverse FFT
      t3 = MPI_WTIME()
      call decomp_2d_fft_3d(out, in, DECOMP_2D_FFT_BACKWARD)
      t4 = t4 + MPI_WTIME() - t3

      ! normalisation - note 2DECOMP&FFT doesn't normalise
      !$acc kernels
      if (.not. skip_c2c(1)) in = in / real(nx, mytype)
      if (.not. skip_c2c(2)) in = in / real(ny, mytype)
      if (.not. skip_c2c(3)) in = in / real(nz, mytype)
      !$acc end kernels

   end do
#if defined(_GPU)
   ierror = cudaDeviceSynchronize()
#endif

   call MPI_ALLREDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t1 = t1 / dble(nproc) / dble(ntest)
   call MPI_ALLREDUCE(t4, t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t3 = t3 / dble(nproc) / dble(ntest)

   ! checking accuracy
   error = 0._mytype
   !$acc parallel loop default(present) reduction(+:error)
   do k = zst3, zen3
      do j = zst2, zen2
         do i = zst1, zen1
            dr = real(i, mytype) / real(nx, mytype) * real(j, mytype) &
                 / real(ny, mytype) * real(k, mytype) / real(nz, mytype)
            di = dr
            dr = dr - real(in(i, j, k), mytype)
            di = di - aimag(in(i, j, k))
            error = error + sqrt(dr * dr + di * di)
         end do
      end do
   end do
   !$acc end loop
   call MPI_ALLREDUCE(MPI_IN_PLACE, error, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
   error = error / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

   ! Abort if the error is too high
   ! A large enough value is needed for the generic backend
   if (error > epsilon(error) * 50 * ntest) then
      if (nrank == 0) write (*, *) 'error / mesh point: ', error
      call decomp_2d_abort(__FILE__, __LINE__, int(log10(error)), "c2c Z test")
   end if

   if (nrank == 0) then
      write (*, *) 'error / mesh point: ', error
      write (*, *) 'Avg time (sec): ', t1, t3
      n1 = real(nx) * real(ny) * real(nz)
      n1 = n1**(1.d0 / 3.d0)
      ! 5n*log(n) flops per 1D FFT of size n using Cooley-Tukey algorithm
      flops = 5.d0 * n1 * log(n1) / log(2.d0)
      ! 3 sets of 1D FFTs for 3 directions, each having n^2 1D FFTs
      flops = flops * 3.d0 * n1**2
      flops = 2.d0 * flops / (t1 + t3)
      write (*, *) 'GFLOPS : ', flops / 1000.d0**3
      write (*, *) '   '
      write (*, *) 'fft_c2c_z completed '
   end if
   !$acc end data

   deallocate (in, out)
   nullify (ph)
   call decomp_2d_fft_finalize
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program fft_c2c_z_skip
