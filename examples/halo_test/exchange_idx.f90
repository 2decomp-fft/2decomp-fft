!! SPDX-License-Identifier: BSD-3-Clause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example tests that the exchanges place data correctly by exchanging global indices on an
! extended (with halos) domain.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program exchange_idx

   use mpi
#if defined(_GPU)
   use cudafor
   use openacc
#endif

   use m_halo
   use decomp_2d

   use decomp_2d_testing

   implicit none

   integer :: nranks_tot
   integer :: nx, ny, nz
   integer, parameter :: nx_base = 5, ny_base = 6, nz_base = 7
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain

   integer :: ierror

   logical :: all_pass

   integer :: ipencil
   integer, dimension(3) :: isize, istart, iend

  !! Initialisation
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

   call decomp_2d_init(nx, ny, nz, p_row, p_col)
   call decomp_2d_testing_log()

   all_pass = .true.
   do ipencil = 1, 3
      all_pass = all_pass .and. run_test(ipencil, [nx, ny, nz])
   end do

   call decomp_2d_finalize()
   if (.not. all_pass) then
      call decomp_2d_abort(1, "Error in exchange indices")
   end if

   call MPI_Finalize(ierror)

contains

   logical function run_test(ipencil, ns) result(test_pass)
      integer, intent(in) :: ipencil
      integer, dimension(3), intent(in) :: ns

      real(mytype), dimension(:, :, :), allocatable :: a1
#if defined(_GPU)
      attributes(device) :: a1
#endif
      type(halo_extents_t) :: halo_extents

      integer :: i, j, k
      integer :: ii, jj, kk
      integer :: idx
      integer :: hx, hy, hz
      integer :: nx, ny, nz

      nx = ns(1)
      ny = ns(2)
      nz = ns(3)

      test_pass = .true.

      if (ipencil == 1) then
         isize(:) = xsize(:)
         istart(:) = xstart(:)
         iend(:) = xend(:)
      else if (ipencil == 2) then
         isize(:) = ysize(:)
         istart(:) = ystart(:)
         iend(:) = yend(:)
      else if (ipencil == 3) then
         isize(:) = zsize(:)
         istart(:) = zstart(:)
         iend(:) = zend(:)
      else
         call decomp_2d_abort(1, "Error in exchange indices")
      end if

      ! Create halo-extended arrays
      hx = 1; hy = 1; hz = 1

      halo_extents = init_halo_extents(ipencil, isize, decomp_main, [hx, hy, hz], .false.)

      if (ipencil == 1) then
         call alloc_x(a1, opt_global=.false., opt_levels=[hx, hy, hz])
      else if (ipencil == 2) then
         call alloc_y(a1, opt_global=.false., opt_levels=[hx, hy, hz])
      else if (ipencil == 3) then
         call alloc_z(a1, opt_global=.false., opt_levels=[hx, hy, hz])
      else
         call decomp_2d_abort(1, "Error in exchange indices")
      end if

      a1(:, :, :) = -1.0_mytype
      !$acc data copyin(istart, iend, halo_extents)
      !$acc kernels default(present)
      do k = halo_extents%zs + hz, halo_extents%ze - hz
         do j = halo_extents%ys + hy, halo_extents%ye - hy
            do i = halo_extents%xs + hx, halo_extents%xe - hx
               ! Extended indices
               ii = istart(1) + (i - halo_extents%xs)
               jj = istart(2) + (j - halo_extents%ys)
               kk = istart(3) + (k - halo_extents%zs)

               ! Global index
               idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)
               a1(i, j, k) = idx
            end do
            if (istart(1) == 1) then
               do i = halo_extents%xs, halo_extents%xs + (hx - 1)
                  ! Extended indices
                  ii = istart(1) + (i - halo_extents%xs)
                  jj = istart(2) + (j - halo_extents%ys)
                  kk = istart(3) + (k - halo_extents%zs)

                  ! Global index
                  idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)
                  a1(i, j, k) = idx
               end do
            end if
            if (iend(1) == nx) then
               do i = halo_extents%xe - hx + 1, halo_extents%xe
                  ! Extended indices
                  ii = istart(1) + (i - halo_extents%xs)
                  jj = istart(2) + (j - halo_extents%ys)
                  kk = istart(3) + (k - halo_extents%zs)

                  ! Global index
                  idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)
                  a1(i, j, k) = idx
               end do
            end if
         end do
         if (istart(2) == 1) then
            do j = halo_extents%ys, halo_extents%ys + (hy - 1)
               do i = halo_extents%xs + hx, halo_extents%xe - hx
                  ! Extended indices
                  ii = istart(1) + (i - halo_extents%xs)
                  jj = istart(2) + (j - halo_extents%ys)
                  kk = istart(3) + (k - halo_extents%zs)

                  ! Global index
                  idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)
                  a1(i, j, k) = idx
               end do
            end do
         end if
         if (iend(2) == ny) then
            do j = halo_extents%ye - hy + 1, halo_extents%ye
               do i = halo_extents%xs + hx, halo_extents%xe - hx
                  ! Extended indices
                  ii = istart(1) + (i - halo_extents%xs)
                  jj = istart(2) + (j - halo_extents%ys)
                  kk = istart(3) + (k - halo_extents%zs)

                  ! Global index
                  idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)
                  a1(i, j, k) = idx
               end do
            end do
         end if
      end do
      !$acc end kernels
      if (istart(3) == 1) then
         !$acc kernels default(present)
         do k = halo_extents%zs, halo_extents%zs + (hz - 1)
            do j = halo_extents%ys + hy, halo_extents%ye - hy
               do i = halo_extents%xs + hx, halo_extents%xe - hx
                  ! Extended indices
                  ii = istart(1) + (i - halo_extents%xs)
                  jj = istart(2) + (j - halo_extents%ys)
                  kk = istart(3) + (k - halo_extents%zs)

                  ! Global index
                  idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)
                  a1(i, j, k) = idx
               end do
            end do
         end do
         !$acc end kernels
      end if
      if (iend(3) == nz) then
         !$acc kernels default(present)
         do k = halo_extents%ze - hz + 1, halo_extents%ze
            do j = halo_extents%ys + hy, halo_extents%ye - hy
               do i = halo_extents%xs + hx, halo_extents%xe - hx
                  ! Extended indices
                  ii = istart(1) + (i - halo_extents%xs)
                  jj = istart(2) + (j - halo_extents%ys)
                  kk = istart(3) + (k - halo_extents%zs)

                  ! Global index
                  idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)
                  a1(i, j, k) = idx
               end do
            end do
         end do
         !$acc end kernels
      end if

      call halo_exchange(a1, ipencil, opt_levels=[hx, hy, hz])
      !$acc kernels default(present)
      do k = halo_extents%zs, halo_extents%ze
         do j = halo_extents%ys, halo_extents%ye
            do i = halo_extents%xs, halo_extents%xe
               ! Extended indices
               ii = istart(1) + (i - halo_extents%xs)
               jj = istart(2) + (j - halo_extents%ys)
               kk = istart(3) + (k - halo_extents%zs)

               ! Global index
               idx = ii + (jj - 1) * (nx + 2 * hx) + (kk - 1) * (nx + 2 * hx) * (ny + 2 * hy)

               ! Check all points are valid, excluding "corners"
               if (.not. is_global_corner(i, j, k, halo_extents, [nx, ny, nz], [hx, hy, hz], istart, iend)) then
                  if (a1(i, j, k) /= idx) then
                     test_pass = .false.
                  end if
               end if
            end do
         end do
      end do
      !$acc end kernels
      !$acc end data

   end function run_test

   ! Determine if we are in a global "corner" of the extended domain
   logical function is_global_corner(i, j, k, halo_extents, ns, hs, istart, iend)
      integer, intent(in) :: i, j, k
      type(halo_extents_t), intent(in) :: halo_extents
      integer, dimension(3), intent(in) :: ns
      integer, dimension(3), intent(in) :: hs
      integer, dimension(3), intent(in) :: istart, iend

      integer :: nx, ny, nz
      integer :: hx, hy, hz

      nx = ns(1)
      ny = ns(2)
      nz = ns(3)
      hx = hs(1)
      hy = hs(2)
      hz = hs(3)

      is_global_corner = .false.

      if (istart(1) == 1 .and. i < halo_extents%xs + hx) then
         if (istart(2) == 1) then
            is_global_corner = is_global_corner .or. (j < halo_extents%ys + hy)
         end if
         if (iend(2) == ny) then
            is_global_corner = is_global_corner .or. (j > halo_extents%ye - hy)
         end if
         if (istart(3) == 1) then
            is_global_corner = is_global_corner .or. (k < halo_extents%zs + hz)
         end if
         if (iend(3) == nz) then
            is_global_corner = is_global_corner .or. (k > halo_extents%ze - hz)
         end if
      end if
      if (iend(1) == nx .and. i > halo_extents%xe - hx) then
         if (istart(2) == 1) then
            is_global_corner = is_global_corner .or. (j < halo_extents%ys + hy)
         end if
         if (iend(2) == ny) then
            is_global_corner = is_global_corner .or. (j > halo_extents%ye - hy)
         end if
         if (istart(3) == 1) then
            is_global_corner = is_global_corner .or. (k < halo_extents%zs + hz)
         end if
         if (iend(3) == nz) then
            is_global_corner = is_global_corner .or. (k > halo_extents%ze - hz)
         end if
      end if

      if (istart(2) == 1 .and. j < halo_extents%ys + hy) then
         if (istart(3) == 1) then
            is_global_corner = is_global_corner .or. (k < halo_extents%zs + hz)
         end if
         if (iend(3) == nz) then
            is_global_corner = is_global_corner .or. (k > halo_extents%ze - hz)
         end if
      end if
      if (iend(2) == ny .and. j > halo_extents%ye - hy) then
         if (istart(3) == 1) then
            is_global_corner = is_global_corner .or. (k < halo_extents%zs + hz)
         end if
         if (iend(3) == nz) then
            is_global_corner = is_global_corner .or. (k > halo_extents%ze - hz)
         end if
      end if

   end function is_global_corner

end program exchange_idx

