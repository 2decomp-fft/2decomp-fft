program fft_c2c_z

   use decomp_2d
   use decomp_2d_fft
   use decomp_2d_constants
   use MPI
#if defined(_GPU)
   use cudafor
   use cufft
   use openacc
#endif

   implicit none

   !integer, parameter :: nx_base=4, ny_base=2, nz_base=3
   integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot

   integer, parameter :: ntest = 10  ! repeat test this times

   type(decomp_info), pointer :: ph => null()
   complex(mytype), allocatable, dimension(:, :, :) :: in, out

   real(mytype) :: dr, di, error, err_all, n1, flops
   integer :: ierror, i, j, k, m
   real(mytype) :: t1, t2, t3, t4

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot/4) + 1
   nx = nx_base*resize_domain
   ny = ny_base*resize_domain
   nz = nz_base*resize_domain
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

   ! initilise input
   do k = zstart(3), zend(3)
      do j = zstart(2), zend(2)
         do i = zstart(1), zend(1)
            dr = real(i, mytype)/real(nx, mytype)*real(j, mytype) &
                 /real(ny, mytype)*real(k, mytype)/real(nz, mytype)
            di = dr
            in(i, j, k) = cmplx(dr, di, mytype)
         end do
      end do
   end do

   t2 = 0._mytype
   t4 = 0._mytype
   !$acc data copyin(in,zstart,zend) copy(out)
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
      in = in/real(nx, mytype)/real(ny, mytype)/real(nz, mytype)
      !$acc end kernels

   end do
#if defined(_GPU)
   ierror = cudaDeviceSynchronize()
#endif

   call MPI_ALLREDUCE(t2, t1, 1, real_type, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t1 = t1/real(nproc, mytype)
   call MPI_ALLREDUCE(t4, t3, 1, real_type, MPI_SUM, &
                      MPI_COMM_WORLD, ierror)
   t3 = t3/real(nproc, mytype)

   ! checking accuracy
   error = 0._mytype
   !$acc parallel loop default(present) reduction(+:error)
   do k = zstart(3), zend(3)
      do j = zstart(2), zend(2)
         do i = zstart(1), zend(1)
            dr = real(i, mytype)/real(nx, mytype)*real(j, mytype) &
                 /real(ny, mytype)*real(k, mytype)/real(nz, mytype)
            di = dr
            dr = dr - real(in(i, j, k), mytype)
            di = di - aimag(in(i, j, k))
            error = error + sqrt(dr*dr + di*di)
         end do
      end do
   end do
  !$acc end loop
   call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
   err_all = err_all/real(nx, mytype)/real(ny, mytype)/real(nz, mytype)

   if (nrank == 0) then
      write (*, *) '===== c2c interface ====='
      write (*, *) 'error / mesh point: ', err_all
      write (*, *) 'time (sec): ', t1, t3
      n1 = real(nx, mytype)*real(ny, mytype)*real(nz, mytype)
      n1 = n1**(1._mytype/3._mytype)
      ! 5n*log(n) flops per 1D FFT of size n using Cooley-Tukey algorithm
      flops = 5._mytype*n1*log(n1)/log(2.0_mytype)
      ! 3 sets of 1D FFTs for 3 directions, each having n^2 1D FFTs
      flops = flops*3._mytype*n1**2
      flops = 2._mytype*flops/((t1 + t3)/real(NTEST, mytype))
      write (*, *) 'GFLOPS : ', flops/1000._mytype**3
   end if
   !$acc end data

   deallocate (in, out)
   nullify (ph)
   call decomp_2d_fft_finalize
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program fft_c2c_z
