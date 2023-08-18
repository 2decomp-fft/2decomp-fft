!! SPDX-License-Identifier: BSD-3-Clause
program fft_c2c_z

   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_constants
   use decomp_2d_mpi
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
   integer :: nargin, arg, FNLength, status, DecInd
   character(len=80) :: InputFN

   integer :: ntest = 10  ! repeat test this times

   type(decomp_info), pointer :: ph => null()
   complex(mytype), allocatable, dimension(:, :, :) :: in, out

   real(mytype) :: dr, di, error, err_all
   integer :: ierror, i, j, k, m
   integer :: zst1, zst2, zst3
   integer :: zen1, zen2, zen3
   double precision :: n1, flops, t1, t2, t3, t4

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
   if ((nargin == 0) .or. (nargin == 2) .or. (nargin == 5) .or. (nargin == 6)) then
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
         elseif (arg == 6) then
            ntest = DecInd
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
         print *, "or 5-6 inputs as"
         print *, "  1) p_row (default=0)"
         print *, "  2) p_col (default=0)"
         print *, "  3) nx "
         print *, "  4) ny "
         print *, "  5) nz "
         print *, "  6) n iterations (optional)"
         print *, "Number of inputs is not correct and the defult settings"
         print *, "will be used"
      end if
   end if

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Test the c2c interface
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call decomp_2d_fft_init(PHYSICAL_IN_Z) ! non-default Z-pencil input

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
      write (*, *) 'It 0 time (sec): ', t1, t3
   endif
   ! Init the time
   t2 = 0.d0
   t4 = 0.d0
   !$acc kernels
   in = in / real(nx, mytype) / real(ny, mytype) / real(nz, mytype)
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
      in = in / real(nx, mytype) / real(ny, mytype) / real(nz, mytype)
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
   call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
   err_all = err_all / real(nx, mytype) / real(ny, mytype) / real(nz, mytype)

   if (nrank == 0) then
      write (*, *) 'error / mesh point: ', err_all
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

end program fft_c2c_z
