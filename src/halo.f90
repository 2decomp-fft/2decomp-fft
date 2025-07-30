!! SPDX-License-Identifier: BSD-3-Clause

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Halo cell support for neighbouring pencils to exchange data
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module m_halo

   use mpi

   use decomp_2d_constants
   use decomp_2d_mpi
   use m_info

   implicit none

   interface update_halo
      module procedure update_halo_real
      module procedure update_halo_real_short
      module procedure update_halo_complex
      module procedure update_halo_complex_short
   end interface update_halo

   ! define neighboring blocks (to be used in halo-cell support)
   !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
   ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom
   integer, save, dimension(3, 6) :: neighbour

   private
   public :: update_halo
   public :: init_neighbour

contains

   !---------------------------------------------------------------------
   ! To support halo-cell exchange:
   !   find the MPI ranks of neighbouring pencils
   !---------------------------------------------------------------------
   subroutine init_neighbour

      integer :: ierror

      ! For X-pencil
      neighbour(1, 1) = MPI_PROC_NULL               ! east
      neighbour(1, 2) = MPI_PROC_NULL               ! west
      call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 0, 1, &
                          neighbour(1, 4), neighbour(1, 3), ierror) ! north & south
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
      call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X, 1, 1, &
                          neighbour(1, 6), neighbour(1, 5), ierror) ! top & bottom
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

      ! For Y-pencil
      call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 0, 1, &
                          neighbour(2, 2), neighbour(2, 1), ierror) ! east & west
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
      neighbour(2, 3) = MPI_PROC_NULL               ! north
      neighbour(2, 4) = MPI_PROC_NULL               ! south
      call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Y, 1, 1, &
                          neighbour(2, 6), neighbour(2, 5), ierror) ! top & bottom
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")

      ! For Z-pencil
      call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 0, 1, &
                          neighbour(3, 2), neighbour(3, 1), ierror) ! east & west
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
      call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_Z, 1, 1, &
                          neighbour(3, 4), neighbour(3, 3), ierror) ! north & south
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_CART_SHIFT")
      neighbour(3, 5) = MPI_PROC_NULL               ! top
      neighbour(3, 6) = MPI_PROC_NULL               ! bottom
      return

   end subroutine init_neighbour

   subroutine update_halo_real_short(in, out, level, opt_global, opt_pencil)

      implicit none

      integer, intent(IN) :: level      ! levels of halo cells required
      real(mytype), dimension(:, :, :), intent(IN) :: in
      real(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
      attributes(device) :: out
#endif
      logical, optional :: opt_global
      integer, intent(in), optional :: opt_pencil

      call update_halo(in, out, level, decomp_main, opt_global, opt_pencil)

   end subroutine update_halo_real_short

   subroutine update_halo_real(in, out, level, decomp, opt_global, opt_pencil)

      implicit none

      integer, intent(IN) :: level      ! levels of halo cells required
      real(mytype), dimension(:, :, :), intent(IN) :: in
      real(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
      attributes(device) :: out
#endif
      TYPE(DECOMP_INFO), intent(in) :: decomp
      logical, optional :: opt_global
      integer, intent(in), optional :: opt_pencil

      logical :: global

      ! starting/ending index of array with halo cells
      integer :: xs, ys, zs, xe, ye, ze
      ! additional start end
      integer :: ist, ien, jst, jen, kst, ken

      integer :: i, j, k, s1, s2, s3, ierror
      integer :: data_type

      integer :: icount, ilength, ijump
      integer :: halo12, halo21, halo31, halo32
      integer, dimension(4) :: requests
      integer, dimension(MPI_STATUS_SIZE, 4) :: status
      integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

      integer :: ipencil
      logical, save :: first_call_x = .true., first_call_y = .true., first_call_z = .true.

      data_type = real_type

#include "halo_common.f90"

      return
   end subroutine update_halo_real

   subroutine update_halo_complex_short(in, out, level, opt_global, opt_pencil)

      implicit none

      integer, intent(IN) :: level      ! levels of halo cells required
      complex(mytype), dimension(:, :, :), intent(IN) :: in
      complex(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
      attributes(device) :: out
#endif
      logical, optional :: opt_global
      integer, intent(in), optional :: opt_pencil

      call update_halo(in, out, level, decomp_main, opt_global, opt_pencil)

   end subroutine update_halo_complex_short

   subroutine update_halo_complex(in, out, level, decomp, opt_global, opt_pencil)

      implicit none

      integer, intent(IN) :: level      ! levels of halo cells required
      complex(mytype), dimension(:, :, :), intent(IN) :: in
      complex(mytype), allocatable, dimension(:, :, :), intent(OUT) :: out
#if defined(_GPU)
      attributes(device) :: out
#endif
      TYPE(DECOMP_INFO), intent(in) :: decomp
      logical, optional :: opt_global
      integer, intent(in), optional :: opt_pencil

      logical :: global

      ! starting/ending index of array with halo cells
      integer :: xs, ys, zs, xe, ye, ze
      ! additional start end
      integer :: ist, ien, jst, jen, kst, ken

      integer :: i, j, k, s1, s2, s3, ierror
      integer :: data_type

      integer :: icount, ilength, ijump
      integer :: halo12, halo21, halo31, halo32
      integer, dimension(4) :: requests
      integer, dimension(MPI_STATUS_SIZE, 4) :: status
      integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

      integer :: ipencil
      logical, save :: first_call_x = .true., first_call_y = .true., first_call_z = .true.

      data_type = complex_type

#include "halo_common.f90"

      return
   end subroutine update_halo_complex

end module m_halo
