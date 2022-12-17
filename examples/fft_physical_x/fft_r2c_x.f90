program fft_r2c_x

   use decomp_2d
   use decomp_2d_fft
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

   integer, parameter :: ntest = 10  ! repeat test this times

   type(decomp_info), pointer :: ph => null(), sp => null()
   complex(mytype), allocatable, dimension(:, :, :) :: out
   real(mytype), allocatable, dimension(:, :, :) :: in_r

   real(mytype) :: dr, error, err_all
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
   ! Test the r2c/c2r interface
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   call decomp_2d_fft_init(PHYSICAL_IN_X)

   ph => decomp_2d_fft_get_ph()
   sp => decomp_2d_fft_get_sp()
   !  input is X-pencil data
   ! output is Z-pencil data
   call alloc_x(in_r, ph, .true.)
   call alloc_z(out, sp, .true.)

   ! initilise input
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            in_r(i, j, k) = real(i, mytype)/real(nx, mytype)*real(j, mytype) &
                            /real(ny, mytype)*real(k, mytype)/real(nz, mytype)
         end do
      end do
   end do

   t2 = 0._mytype
   t4 = 0._mytype
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
      in_r = in_r/(real(nx, mytype)*real(ny, mytype)*real(nz, mytype))
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
   do k = xstart(3), xend(3)
      do j = xstart(2), xend(2)
         do i = xstart(1), xend(1)
            dr = real(i, mytype)/real(nx, mytype)*real(j, mytype) &
                 /real(ny, mytype)*real(k, mytype)/real(nz, mytype)
            error = error + abs(in_r(i, j, k) - dr)
         end do
      end do
   end do
   call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
   err_all = err_all/real(nx, mytype)/real(ny, mytype)/real(nz, mytype)

   if (nrank == 0) then
      write (*, *) '===== r2c/c2r interface ====='
      write (*, *) 'error / mesh point: ', err_all
      write (*, *) 'time (sec): ', t1, t3
   end if

   deallocate (in_r, out)
   nullify (ph)
   nullify (sp)
   call decomp_2d_fft_finalize
   call decomp_2d_finalize
   call MPI_FINALIZE(ierror)

end program fft_r2c_x

