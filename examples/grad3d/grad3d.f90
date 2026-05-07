!! SPDX-License-Identifier: BSD-3-Clause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the gradient of a periodic field using global
! transposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program grad3d

   use mpi

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_testing
   use MPI
#if defined(_OPENMP_GPU)
   use omp_lib
#elif defined(_GPU)
   use cudafor
   use openacc
#endif

   implicit none

   integer, parameter :: nx_base = 64, ny_base = 64, nz_base = 64
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer :: ierror
   logical :: all_pass
   logical, parameter :: do_io = .false.
   double precision :: t0
   double precision :: t_total_start, t_total_end
   double precision :: t_compute_start, t_compute_end
   double precision :: t_io_start, t_io_end
   double precision :: t_init_phi
   double precision :: t_derx, t_dery, t_derz
   double precision :: t_x2y, t_y2x_y, t_y2z, t_z2y, t_y2x_z
   double precision :: t_test_x, t_test_y, t_test_z

   real(mytype), parameter :: lx = 1.0_mytype
   real(mytype), parameter :: ly = 1.0_mytype
   real(mytype), parameter :: lz = 1.0_mytype

   real(mytype) :: dx, dy, dz
   real(mytype) :: error_ref

   real(mytype), allocatable, dimension(:, :, :) :: phi1, phi2, phi3
   real(mytype), allocatable, dimension(:, :, :) :: dphiX, dphiY, dphiz
   real(mytype), allocatable, dimension(:, :, :) :: wk2, wk3

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

   ! Avoid writing decomp_2d_setup.log from this timing-oriented example.
   decomp_log = D2D_LOG_STDOUT

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   dx = lx / real(nx, mytype)
   dy = ly / real(ny, mytype)
   dz = lz / real(nz, mytype)

   if (nrank == 0) then
      write (*, *) '-----------------------------------------------'
      write (*, *) "Mesh Resolution ", nx, ny, nz
      if (.not. do_io) write (*, *) "I/O output disabled in grad3d"
      write (*, *) '-----------------------------------------------'
   end if

   call allocate_var()
   call reset_timers()
   call MPI_BARRIER(MPI_COMM_WORLD, ierror)
   t_total_start = MPI_WTIME()
#if defined(_OPENMP_GPU)
   !$omp target data map(alloc: phi1,dphiX,dphiY,dphiZ,phi2,phi3,wk2,wk3)
#elif defined(_GPU)
   !$acc data create(phi2,phi3,wk2,wk3) copy(phi1,dphiX, dphiY, dphiZ)
#endif
   call MPI_BARRIER(MPI_COMM_WORLD, ierror)
   t_compute_start = MPI_WTIME()
   t0 = MPI_WTIME()
   call init_phi()
   t_init_phi = t_init_phi + (MPI_WTIME() - t0)
   call compute_grad()
   t0 = MPI_WTIME()
   call test_derX(dphiX)
   t_test_x = t_test_x + (MPI_WTIME() - t0)
   t0 = MPI_WTIME()
   call test_derY(dphiY)
   t_test_y = t_test_y + (MPI_WTIME() - t0)
   t0 = MPI_WTIME()
   call test_derZ(dphiZ)
   t_test_z = t_test_z + (MPI_WTIME() - t0)
   call MPI_BARRIER(MPI_COMM_WORLD, ierror)
   t_compute_end = MPI_WTIME()

   if (do_io) then
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      t_io_start = MPI_WTIME()
      call write_data()
      call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      t_io_end = MPI_WTIME()
   else
      t_io_start = t_compute_end
      t_io_end = t_compute_end
   end if
#if defined(_OPENMP_GPU)
   !$omp end target data
#elif defined(_GPU)
   !$acc end data
#endif
   call MPI_BARRIER(MPI_COMM_WORLD, ierror)
   t_total_end = MPI_WTIME()

   if (nrank == 0) then
      write (*, *) '-----------------------------------------------'
      write (*, *) "End GRAD calculation check all pass: ", all_pass
   end if
   call print_timing("init_phi", t_init_phi)
   call print_timing("derx", t_derx)
   call print_timing("transpose_x_to_y", t_x2y)
   call print_timing("dery", t_dery)
   call print_timing("transpose_y_to_x for dphiY", t_y2x_y)
   call print_timing("transpose_y_to_z", t_y2z)
   call print_timing("derz", t_derz)
   call print_timing("transpose_z_to_y", t_z2y)
   call print_timing("transpose_y_to_x for dphiZ", t_y2x_z)
   call print_timing("y->z pack", d2d_t_y2z_pack)
   call print_timing("y->z mpi", d2d_t_y2z_mpi)
   call print_timing("y->z unpack", d2d_t_y2z_unpack)
   call print_timing("z->y pack", d2d_t_z2y_pack)
   call print_timing("z->y mpi", d2d_t_z2y_mpi)
   call print_timing("z->y unpack", d2d_t_z2y_unpack)
   call print_counter_pair("y->z use_device_ptr hits / MPI calls", &
                           d2d_n_y2z_use_device_ptr, d2d_n_y2z_mpi_calls)
   call print_counter_pair("z->y use_device_ptr hits / MPI calls", &
                           d2d_n_z2y_use_device_ptr, d2d_n_z2y_mpi_calls)
   call print_timing("test_derX", t_test_x)
   call print_timing("test_derY", t_test_y)
   call print_timing("test_derZ", t_test_z)
   call print_timing("compute+validation total", t_compute_end - t_compute_start)
   call print_timing("I/O total", t_io_end - t_io_start)
   call print_timing("grad3d wall time", t_total_end - t_total_start)
   if (nrank == 0) write (*, *) '==============================================='

   call finalize()

   call decomp_2d_finalize

   if (.not. all_pass) call decomp_2d_abort(1, "Error in grad3d")

   call MPI_FINALIZE(ierror)

contains

   subroutine reset_timers()

      implicit none

      t_init_phi = 0.0d0
      t_derx = 0.0d0
      t_dery = 0.0d0
      t_derz = 0.0d0
      t_x2y = 0.0d0
      t_y2x_y = 0.0d0
      t_y2z = 0.0d0
      t_z2y = 0.0d0
      t_y2x_z = 0.0d0
      t_test_x = 0.0d0
      t_test_y = 0.0d0
      t_test_z = 0.0d0

   end subroutine reset_timers

   subroutine print_timing(label, value)

      implicit none

      character(len=*), intent(in) :: label
      double precision, intent(in) :: value

      double precision :: value_max

      call MPI_REDUCE(value, value_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
      if (nrank == 0) write (*, '(A,1X,A,1X,ES12.5)') "Timing max over ranks [s]:", trim(label), value_max

   end subroutine print_timing

   subroutine print_counter_pair(label, value1, value2)

      implicit none

      character(len=*), intent(in) :: label
      integer, intent(in) :: value1, value2

      integer :: value1_sum, value2_sum

      call MPI_REDUCE(value1, value1_sum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      call MPI_REDUCE(value2, value2_sum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror)
      if (nrank == 0) write (*, '(A,1X,A,1X,I0,1X,I0)') "Timing counter sum over ranks:", trim(label), value1_sum, value2_sum

   end subroutine print_counter_pair

   !=====================================================================
   ! Initialize
   !=====================================================================
   subroutine allocate_var()

      use decomp_2d

      implicit none

      logical, parameter :: glob = .false.

      ! Allocate main variables in X-pencil
      call alloc_x(phi1, glob)
      call alloc_x(dphiX, glob)
      call alloc_x(dphiY, glob)
      call alloc_x(dphiZ, glob)

      ! Working array used more than once
      call alloc_y(phi2, glob)
      call alloc_y(wk2, glob)

      call alloc_z(phi3, glob)
      call alloc_z(wk3, glob)

      all_pass = .true.

   end subroutine allocate_var

   !=====================================================================
   ! Initialize the scalar field
   !=====================================================================
   subroutine init_phi()

      implicit none

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Scalar field
#if defined(_OPENMP_GPU)
      !$omp target teams distribute parallel do collapse(3) private(x,y,z) map(present,alloc:phi1)
#elif defined(_GPU)
      !$acc kernels default(present)
#endif
      do k = 1, xe3
         do j = 1, xe2
            do i = 1, xe1
               z = (k + xs3 - 2) * dz
               y = (j + xs2 - 2) * dy
               x = (i + xs1 - 2) * dx
               phi1(i, j, k) = -2._mytype * cos(twopi * (x / lx)) * cos(twopi * (y / ly)) * sin(twopi * (z / lz))
            end do
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
#elif defined(_GPU)
      !$acc end kernels
#endif

   end subroutine init_phi
   !=====================================================================
   ! Finalize with deallocation of arrays
   !=====================================================================
   subroutine finalize()

      implicit none

      deallocate (phi1, phi2, phi3)
      deallocate (dphiX, dphiY, dphiZ)
      deallocate (wk2, wk3)

   end subroutine finalize
   !=====================================================================
   ! Calculate gradient using global transposition
   !=====================================================================
   subroutine compute_grad()

      implicit none
      double precision :: t

      ! Compute X derivative
      t = MPI_WTIME()
      call derx(dphiX, phi1, dx, xsize(1), xsize(2), xsize(3))
      t_derx = t_derx + (MPI_WTIME() - t)

      ! Compute Y derivative
      t = MPI_WTIME()
      call transpose_x_to_y(phi1, phi2)
      t_x2y = t_x2y + (MPI_WTIME() - t)
      t = MPI_WTIME()
      call dery(wk2, phi2, dy, ysize(1), ysize(2), ysize(3))
      t_dery = t_dery + (MPI_WTIME() - t)
      t = MPI_WTIME()
      call transpose_y_to_x(wk2, dphiY)
      t_y2x_y = t_y2x_y + (MPI_WTIME() - t)

      ! Compute Z derivative
      t = MPI_WTIME()
      call transpose_y_to_z(phi2, phi3)
      t_y2z = t_y2z + (MPI_WTIME() - t)
      t = MPI_WTIME()
      call derz(wk3, phi3, dz, zsize(1), zsize(2), zsize(3))
      t_derz = t_derz + (MPI_WTIME() - t)
      t = MPI_WTIME()
      call transpose_z_to_y(wk3, wk2)
      t_z2y = t_z2y + (MPI_WTIME() - t)
      t = MPI_WTIME()
      call transpose_y_to_x(wk2, dphiZ)
      t_y2x_z = t_y2x_z + (MPI_WTIME() - t)

   end subroutine compute_grad
   !=====================================================================
   ! Calculate gradient in X (data in X-pencil)
   !=====================================================================
   subroutine derx(df, ff, delta, nx, ny, nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx, ny, nz) :: df
      real(mytype), intent(in), dimension(nx, ny, nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff / delta

#if defined(_OPENMP_GPU)
      !$omp target teams distribute parallel do collapse(2) map(present,alloc:df,ff)
#elif defined(_GPU)
      !$acc kernels default(present)
#endif
      do k = 1, nz
         do j = 1, ny
            df(1, j, k) = coeff * (ff(2, j, k) - ff(nx, j, k))
            do i = 2, nx - 1
               df(i, j, k) = coeff * (ff(i + 1, j, k) - ff(i - 1, j, k))
            end do
            df(nx, j, k) = coeff * (ff(1, j, k) - ff(nx - 1, j, k))
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
#elif defined(_GPU)
      !$acc end kernels
#endif

   end subroutine derx
   !=====================================================================
   ! Calculate gradient in Y (data in Y-pencil)
   !=====================================================================
   subroutine dery(df, ff, delta, nx, ny, nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx, ny, nz) :: df
      real(mytype), intent(in), dimension(nx, ny, nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff / delta

#if defined(_OPENMP_GPU)
      !$omp target teams distribute parallel do collapse(3) map(present,alloc:df,ff)
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               if (j == 1) then
                  df(i, j, k) = coeff * (ff(i, 2, k) - ff(i, ny, k))
               else if (j == ny) then
                  df(i, j, k) = coeff * (ff(i, 1, k) - ff(i, ny - 1, k))
               else
                  df(i, j, k) = coeff * (ff(i, j + 1, k) - ff(i, j - 1, k))
               end if
            end do
         end do
      end do
      !$omp end target teams distribute parallel do
#elif defined(_GPU)
      !$acc kernels default(present)
      do k = 1, nz
         do i = 1, nx
            df(i, 1, k) = coeff * (ff(i, 2, k) - ff(i, ny, k))
         end do
         do j = 2, ny - 1
            do i = 1, nx
               df(i, j, k) = coeff * (ff(i, j + 1, k) - ff(i, j - 1, k))
            end do
         end do
         do i = 1, nx
            df(i, ny, k) = coeff * (ff(i, 1, k) - ff(i, ny - 1, k))
         end do
      end do
      !$acc end kernels
#else
      do k = 1, nz
         do i = 1, nx
            df(i, 1, k) = coeff * (ff(i, 2, k) - ff(i, ny, k))
         end do
         do j = 2, ny - 1
            do i = 1, nx
               df(i, j, k) = coeff * (ff(i, j + 1, k) - ff(i, j - 1, k))
            end do
         end do
         do i = 1, nx
            df(i, ny, k) = coeff * (ff(i, 1, k) - ff(i, ny - 1, k))
         end do
      end do
#endif

   end subroutine dery
   !=====================================================================
   ! Calculate gradient in Z (data in Z-pencil)
   !=====================================================================
   subroutine derz(df, ff, delta, nx, ny, nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx, ny, nz) :: df
      real(mytype), intent(in), dimension(nx, ny, nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff / delta

#if defined(_OPENMP_GPU)
      !$omp target teams distribute parallel do collapse(2) map(present,alloc:df,ff)
#elif defined(_GPU)
      !$acc kernels default(present)
#endif
      do j = 1, ny
         do i = 1, nx
            df(i, j, 1) = coeff * (ff(i, j, 2) - ff(i, j, nz))
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
      !$omp target teams distribute parallel do collapse(3) map(present,alloc:df,ff)
#endif
      do k = 2, nz - 1
         do j = 1, ny
            do i = 1, nx
               df(i, j, k) = coeff * (ff(i, j, k + 1) - ff(i, j, k - 1))
            end do
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
      !$omp target teams distribute parallel do collapse(2) map(present,alloc:df,ff)
#endif
      do j = 1, ny
         do i = 1, nx
            df(i, j, nz) = coeff * (ff(i, j, 1) - ff(i, j, nz - 1))
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
#elif defined(_GPU)
      !$acc end kernels
#endif

   end subroutine derz
   !=====================================================================
   ! Test derivatives against analytical solution (data in X-pencil)
   !=====================================================================
   subroutine test_derX(df)

      implicit none
      ! Arguments
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: df

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      real(mytype) :: dphi, dphi_num
      real(mytype) :: error = 0._mytype
      real(mytype) :: err_all = 0._mytype
      real(mytype) :: dphi2 = 0._mytype
      real(mytype) :: sum_dphi2 = 0._mytype
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Compute the error against analytical solution
#if defined(_OPENMP_GPU)
      !$omp target teams distribute parallel do collapse(3) reduction(+:error,dphi2) &
      !$omp& private(x,y,z,dphi,dphi_num) map(present,alloc:df)
#elif defined(_GPU)
      !$acc parallel loop default(present) reduction(+:error)
#endif
      do k = 1, xe3
         do j = 1, xe2
            do i = 1, xe1
               z = (k + xs3 - 2) * dz
               y = (j + xs2 - 2) * dy
               x = (i + xs1 - 2) * dx
               dphi = -2._mytype * (twopi / lx) &
                      * sin(twopi * (x / lx)) * cos(twopi * (y / ly)) * sin(twopi * (z / lz))
               dphi2 = dphi2 + dphi * dphi
               dphi_num = df(i, j, k)
               error = error + (dphi - dphi_num) * (dphi - dphi_num)
            end do
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
#elif defined(_GPU)
      !$acc end loop
#endif
      call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(dphi2, sum_dphi2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      err_all = sqrt(err_all / sum_dphi2) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
      error_ref = err_all

      if (nrank == 0) then
         write (*, *) 'DX error / mesh point: ', err_all
      end if

   end subroutine test_derX

   !=====================================================================
   ! Test derivatives against analytical solution (data in Y-pencil)
   !=====================================================================
   subroutine test_derY(df)

      implicit none
      ! Arguments
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: df

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      real(mytype) :: dphi, dphi_num
      real(mytype) :: error = 0._mytype
      real(mytype) :: err_all = 0._mytype
      real(mytype) :: dphi2 = 0._mytype
      real(mytype) :: sum_dphi2 = 0._mytype
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Compute the error against analytical solution
#if defined(_OPENMP_GPU)
      !$omp target teams distribute parallel do collapse(3) reduction(+:error,dphi2) &
      !$omp& private(x,y,z,dphi,dphi_num) map(present,alloc:df)
#elif defined(_GPU)
      !$acc parallel loop default(present) reduction(+:error)
#endif
      do k = 1, xe3
         do j = 1, xe2
            do i = 1, xe1
               z = (k + xs3 - 2) * dz
               y = (j + xs2 - 2) * dy
               x = (i + xs1 - 2) * dx
               dphi = -2._mytype * (twopi / ly) &
                      * cos(twopi * (x / lx)) * sin(twopi * (y / ly)) * sin(twopi * (z / lz))
               dphi2 = dphi2 + dphi * dphi
               dphi_num = df(i, j, k)
               error = error + (dphi - dphi_num) * (dphi - dphi_num)
            end do
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
#elif defined(_GPU)
      !$acc end loop
#endif
      call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(dphi2, sum_dphi2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      err_all = sqrt(err_all / sum_dphi2) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

      if (nrank == 0) then
         write (*, *) 'DY error / mesh point: ', err_all
      end if

      if (abs(err_all - error_ref) > 1.0e-5_mytype) all_pass = .false.

   end subroutine test_derY

   !=====================================================================
   ! Test derivatives against analytical solution (data in Z-pencil)
   !=====================================================================
   subroutine test_derZ(df)

      implicit none
      ! Arguments
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: df

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      real(mytype) :: dphi, dphi_num
      real(mytype) :: error = 0._mytype
      real(mytype) :: err_all = 0._mytype
      real(mytype) :: dphi2 = 0._mytype
      real(mytype) :: sum_dphi2 = 0._mytype
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Compute the error against analytical solution
#if defined(_OPENMP_GPU)
      !$omp target teams distribute parallel do collapse(3) reduction(+:error,dphi2) &
      !$omp& private(x,y,z,dphi,dphi_num) map(present,alloc:df)
#elif defined(_GPU)
      !$acc parallel loop default(present) reduction(+:error)
#endif
      do k = 1, xe3
         do j = 1, xe2
            do i = 1, xe1
               z = (k + xs3 - 2) * dz
               y = (j + xs2 - 2) * dy
               x = (i + xs1 - 2) * dx
               dphi = 2._mytype * (twopi / lz) &
                      * cos(twopi * (x / lx)) * cos(twopi * (y / ly)) * cos(twopi * (z / lz))
               dphi2 = dphi2 + dphi * dphi
               dphi_num = df(i, j, k)
               error = error + (dphi - dphi_num) * (dphi - dphi_num)
            end do
         end do
      end do
#if defined(_OPENMP_GPU)
      !$omp end target teams distribute parallel do
#elif defined(_GPU)
      !$acc end loop
#endif
      call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(dphi2, sum_dphi2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      err_all = sqrt(err_all / sum_dphi2) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

      if (nrank == 0) then
         write (*, *) 'DZ error / mesh point: ', err_all
      end if

      if (abs(err_all - error_ref) > 1.0e-5_mytype) all_pass = .false.

   end subroutine test_derZ
   !=====================================================================
   ! Write of the results (all data are in X-pencil)
   !=====================================================================
   subroutine write_data()

      use decomp_2d_io

      implicit none

      logical :: dir_exists

#if defined(_OPENMP_GPU)
      !$omp target update from(phi1,dphiX,dphiY,dphiZ)
#elif defined(_GPU)
      !$acc update self (phi1)
      !$acc update self (dphiX)
      !$acc update self (dphiY)
      !$acc update self (dphiZ)
#endif
      if (nrank == 0) then
         inquire (file="out", exist=dir_exists)
         if (.not. dir_exists) then
            call execute_command_line("mkdir out 2> /dev/null")
         end if
      end if

      call decomp_2d_io_init()

      ! Standard MPI I/O pattern - 1 file per field
      call decomp_2d_write_one(1, phi1, 'phi1.dat', opt_dirname='out')
      call decomp_2d_write_one(1, dphiX, 'dphiX.dat', opt_dirname='out')
      call decomp_2d_write_one(1, dphiY, 'dphiY.dat', opt_dirname='out')
      call decomp_2d_write_one(1, dphiZ, 'dphiZ.dat', opt_dirname='out')

      call decomp_2d_io_fin

      call write_xdmf()

   end subroutine write_data
   !=====================================================================
   ! Write of the xdmf file to visualise in paraview
   !=====================================================================
   subroutine write_xdmf()
      ! This subroutine is based on the xdmf writers in Xcompact3d.
      ! Copyright (c) 2012-2022, Xcompact3d
      ! SPDX-License-Identifier: BSD 3-Clause
      use, intrinsic :: iso_fortran_env, only: real64

      integer :: ioxdmf, code

      character(len=:), allocatable :: fmt

      integer :: precision

      integer :: varctr
      character(len=16) :: filename
      character(len=5) :: varname

      if (nrank == 0) then
         OPEN (newunit=ioxdmf, file="./out.xdmf", iostat=code)
         if (code /= 0) call decomp_2d_abort(code, "Grad3D : error when opening the file ./out.xdmf")

         write (ioxdmf, '(A22)') '<?xml version="1.0" ?>'
         write (ioxdmf, *) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
         write (ioxdmf, *) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
         write (ioxdmf, *) '<Domain>'

         write (ioxdmf, '(A)') '    <Topology name="topo" TopologyType="3DCoRectMesh"'
         fmt = "(A, I0, A, I0, A, I0, A)"
         write (ioxdmf, fmt) '        Dimensions="', nz, " ", ny, " ", nx, '">'
         write (ioxdmf, '(A)') '    </Topology>'

         write (ioxdmf, *) '    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
         write (ioxdmf, *) '        <!-- Origin -->'
         write (ioxdmf, *) '        <DataItem Format="XML" Dimensions="3">'
         write (ioxdmf, *) '          0.0 0.0 0.0'
         write (ioxdmf, *) '        </DataItem>'
         write (ioxdmf, *) '        <!-- DxDyDz -->'
         write (ioxdmf, *) '        <DataItem Format="XML" Dimensions="3">'
         if (mytype == kind(0._real64)) then
            fmt = "(A, E24.17, A, E24.17, A, E24.17)"
         else
            fmt = "(A, E16.9, A, E16.9, A, E16.9)"
         end if
         write (ioxdmf, fmt) '        ', 1.0_mytype, " ", 1.0_mytype, " ", 1.0_mytype
         write (ioxdmf, *) '        </DataItem>'
         write (ioxdmf, *) '    </Geometry>'

         write (ioxdmf, *) '   <Grid Name="1" GridType="Uniform">'
         write (ioxdmf, *) '       <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
         write (ioxdmf, *) '       <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
         do varctr = 1, 4
            select case (varctr)
            case (1)
               write (varname, "(A)") "phi1"
               write (filename, '(A)') "./out/phi1.dat"
            case (2)
               write (varname, "(A)") "dphiX"
               write (filename, '(A)') "./out/dphiX.dat"
            case (3)
               write (varname, "(A)") "dphiY"
               write (filename, '(A)') "./out/dphiY.dat"
            case (4)
               write (varname, "(A)") "dphiZ"
               write (filename, '(A)') "./out/dphiZ.dat"
            end select
            write (ioxdmf, *) '       <Attribute Name="'//trim(varname)//'" Center="Node">'
            write (ioxdmf, *) '          <DataItem Format="Binary"'

#ifdef DOUBLE_PREC
            print *, "Double precision build"
#ifdef SAVE_SINGLE
            precision = 4
#else
            precision = 8
#endif
#else
            precision = 4
#endif
            write (ioxdmf, "(A,I0,A)") '            DataType="Float" Precision="', precision, '" Endian="little" Seek="0"'

            fmt = "(A, I0, A, I0, A, I0, A)"
            write (ioxdmf, fmt) '            Dimensions="', nz, " ", ny, " ", nx, '">'

            write (ioxdmf, *) '              '//trim(filename)

            write (ioxdmf, *) '           </DataItem>'
            write (ioxdmf, *) '        </Attribute>'
         end do
         write (ioxdmf, '(/)')
         write (ioxdmf, *) '    </Grid>'
         write (ioxdmf, *) '</Domain>'
         write (ioxdmf, '(A7)') '</Xdmf>'
         close (ioxdmf)
      end if

   end subroutine write_xdmf

end program grad3d
