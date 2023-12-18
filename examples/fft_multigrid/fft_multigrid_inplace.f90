!! SPDX-License-Identifier: BSD-3-Clause

program fft_multigrid_inplace

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_fft
   use decomp_2d_mpi
   use decomp_2d_testing
   use MPI
   use mod_ntest
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
   integer :: ngrid = 10
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
   if (ngrid > 0) call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, 1, &
                                          opt_inplace=.false., &
                                          opt_inplace_r2c=.false., &
                                          opt_inplace_c2r=.false.)
   if (ngrid > 1) call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, 2, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.false., &
                                          opt_inplace_c2r=.false.)
   if (ngrid > 2) call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, 3, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.true., &
                                          opt_inplace_c2r=.false.)
   if (ngrid > 3) call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, 4, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.false., &
                                          opt_inplace_c2r=.true.)
   if (ngrid > 4) call decomp_2d_fft_init(PHYSICAL_IN_X, nx, ny, nz, 5, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.true., &
                                          opt_inplace_c2r=.true.)
   if (ngrid > 5) call decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz, 6, &
                                          opt_inplace=.false., &
                                          opt_inplace_r2c=.false., &
                                          opt_inplace_c2r=.false.)
   if (ngrid > 6) call decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz, 7, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.false., &
                                          opt_inplace_c2r=.false.)
   if (ngrid > 7) call decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz, 8, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.true., &
                                          opt_inplace_c2r=.false.)
   if (ngrid > 8) call decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz, 9, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.false., &
                                          opt_inplace_c2r=.true.)
   if (ngrid > 9) call decomp_2d_fft_init(PHYSICAL_IN_Z, nx, ny, nz, 10, &
                                          opt_inplace=.true., &
                                          opt_inplace_r2c=.true., &
                                          opt_inplace_c2r=.true.)

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

end program fft_multigrid_inplace
