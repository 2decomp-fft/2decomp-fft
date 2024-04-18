!! SPDX-License-Identifier: BSD-3-Clause
program fft_r2c_z

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

   type(decomp_info), pointer :: ph => null(), sp => null()
   complex(mytype), allocatable, dimension(:, :, :) :: out
   real(mytype), allocatable, dimension(:, :, :) :: in_r

   real(mytype) :: dr, error
   integer :: ierror, i, j, k, m
   integer :: zst1, zst2, zst3
   integer :: zen1, zen2, zen3
   double precision :: t1, t2, t3, t4
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
   ! Test the r2c/c2r interface
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call decomp_2d_fft_init(PHYSICAL_IN_Z, opt_skip_XYZ_c2c=skip_c2c) ! non-default Z-pencil input

   ph => decomp_2d_fft_get_ph()
   sp => decomp_2d_fft_get_sp()
   !  input is Z-pencil data
   ! output is X-pencil data
   call alloc_z(in_r, ph, .true.)
   call alloc_x(out, sp, .true.)
   zst1 = zstart(1); zen1 = zend(1)
   zst2 = zstart(2); zen2 = zend(2)
   zst3 = zstart(3); zen3 = zend(3)

   ! initilise input
   do k = zst3, zen3
      do j = zst2, zen2
         do i = zst1, zen1
            in_r(i, j, k) = real(i, mytype) / real(nx, mytype) * real(j, mytype) &
                            / real(ny, mytype) * real(k, mytype) / real(nz, mytype)
         end do
      end do
   end do

   !$acc data copyin(in_r) copy(out)
   ! First iterations out of the counting loop
   t1 = MPI_WTIME()
   call decomp_2d_fft_3d(in_r, out)
   t2 = MPI_WTIME() - t1
   t3 = MPI_WTIME()
   call decomp_2d_fft_3d(out, in_r)
   t4 = MPI_WTIME() - t3
   call MPI_ALLREDUCE(t2, t1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t1 = t1 / dble(nproc)
   call MPI_ALLREDUCE(t4, t3, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t3 = t3 / dble(nproc)
   if (nrank == 0) then
      write (*, *) '===== r2c/c2r interface ====='
      write (*, *) 'First iteration with dedicated timer'
      write (*, *) '     time (sec): ', t1, t3
      write (*, *) ''
   end if
   ! Init the time
   t2 = 0.d0
   t4 = 0.d0
   !$acc kernels
   if (.not. skip_c2c(1)) in_r = in_r / real(nx, mytype)
   if (.not. skip_c2c(2)) in_r = in_r / real(ny, mytype)
   in_r = in_r / real(nz, mytype) ! r2c / c2r + physical_in_z
   !$acc end kernels

   do m = 1, ntest

      ! 3D r2c FFT
      t1 = MPI_WTIME()
      call decomp_2d_fft_3d(in_r, out)
      t2 = t2 + MPI_WTIME() - t1

      ! 3D inverse FFT
      t3 = MPI_WTIME()
      call decomp_2d_fft_3d(out, in_r)
      t4 = t4 + MPI_WTIME() - t3

      ! normalisation - note 2DECOMP&FFT doesn't normalise
      !$acc kernels
      if (.not. skip_c2c(1)) in_r = in_r / real(nx, mytype)
      if (.not. skip_c2c(2)) in_r = in_r / real(ny, mytype)
      in_r = in_r / real(nz, mytype) ! r2c / c2r + physical_in_z
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
            error = error + abs(in_r(i, j, k) - dr)
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
      if (nrank == 0) write (*, *) 'Incorrect normalization due to ignored request to skip the r2c/c2r transform'
      call decomp_2d_abort(__FILE__, __LINE__, int(log10(error)), "r2c X test")
   end if

   if (nrank == 0) then
      write (*, *) 'error / mesh point: ', error
      write (*, *) 'Avg time (sec): ', t1, t3
      write (*, *) '   '
      write (*, *) 'fft_r2c_z completed '
      write (*, *) '   '
   end if
   !$acc end data

   deallocate (in_r)
   deallocate (out)
   nullify (ph)
   nullify (sp)
   call decomp_2d_fft_finalize
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program fft_r2c_z
