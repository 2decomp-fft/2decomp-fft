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

   ! define neighboring blocks (to be used in halo-cell support)
   !  first dimension 1=X-pencil, 2=Y-pencil, 3=Z-pencil
   ! second dimension 1=east, 2=west, 3=north, 4=south, 5=top, 6=bottom
   integer, save, dimension(3, 6) :: neighbour

   type :: halo_extents_t
      integer :: xs, xe
      integer :: ys, ye
      integer :: zs, ze
   end type halo_extents_t

   interface halo_extents_t
      procedure init_halo_extents
   end interface halo_extents_t

   interface update_halo
      module procedure update_halo_real
      module procedure update_halo_real_short
      module procedure update_halo_complex
      module procedure update_halo_complex_short
   end interface update_halo

   interface halo_exchange
      procedure halo_exchange_real
      procedure halo_exchange_complex
   end interface halo_exchange

   private
   public :: update_halo 
   public :: init_neighbour
   public :: halo_exchange
   public :: halo_extents_t
  
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
      type(halo_extents_t) :: halo_extents
      
      ! additional start end
      integer :: ist, ien, jst, jen, kst, ken

      integer :: i, j, k, s1, s2, s3

      integer :: ipencil

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
      type(halo_extents_t) :: halo_extents
      
      ! additional start end
      integer :: ist, ien, jst, jen, kst, ken

      integer :: i, j, k, s1, s2, s3

      integer :: ipencil

#include "halo_common.f90"

      return
   end subroutine update_halo_complex

  integer function get_pencil(sizes, decomp, opt_pencil)
    integer, dimension(3), intent(in) :: sizes
    type(decomp_info), intent(in) :: decomp
    integer, intent(in), optional :: opt_pencil

    integer :: s1, s2, s3
    logical, save :: first_call_x = .true.
    logical, save :: first_call_y = .true.
    logical, save :: first_call_z = .true.

    s1 = sizes(1)
    s2 = sizes(2)
    s3 = sizes(3)

    if (present(opt_pencil)) then
       get_pencil = opt_pencil
    else
       ! Historic/default behaviour
       if (s1 == decomp%xsz(1)) then
          get_pencil = 1
          if (first_call_x) then
             first_call_x = .false.
             call decomp_2d_warning(__FILE__, __LINE__, &
                                    0, "Deprecated interface - calling halo in X without explicit pencil")
          end if
       else if (s2 == decomp%ysz(2)) then
          get_pencil = 2
          if (first_call_y) then
             first_call_y = .false.
             call decomp_2d_warning(__FILE__, __LINE__, &
                                    0, "Deprecated interface - calling halo in Y without explicit pencil")
          end if
       else if (s3 == decomp%zsz(3)) then
          get_pencil = 3
          if (first_call_z) then
             first_call_z = .false.
             call decomp_2d_warning(__FILE__, __LINE__, &
                                    0, "Deprecated interface - calling halo in Z without explicit pencil")
          end if
       else
          get_pencil = 0
          call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid decomposition size")
       end if
    end if

  end function get_pencil

  type(halo_extents_t) function init_halo_extents(ipencil, sizes, decomp, level, global) result(halo_extents)
    integer, intent(in) :: ipencil
    integer, dimension(3), intent(in) :: sizes
    type(decomp_info), intent(in) :: decomp
    integer, intent(in) :: level
    logical, intent(in) :: global

    integer :: s1, s2, s3

    s1 = sizes(1)
    s2 = sizes(2)
    s3 = sizes(3)
    
    ! Calculate the starting index and ending index of output
    if (ipencil == 1) then  ! X-pencil input
       if (global) then
          halo_extents%xs = decomp%xst(1)
          halo_extents%xe = decomp%xen(1)
          halo_extents%ys = decomp%xst(2) - level
          halo_extents%ye = decomp%xen(2) + level
          halo_extents%zs = decomp%xst(3) - level
          halo_extents%ze = decomp%xen(3) + level
       else
          halo_extents%xs = 1
          halo_extents%xe = s1
          halo_extents%ys = 1 - level
          halo_extents%ye = s2 + level
          halo_extents%zs = 1 - level
          halo_extents%ze = s3 + level
       end if
    else if (ipencil == 2) then  ! Y-pencil input
       if (global) then
          halo_extents%xs = decomp%yst(1) - level
          halo_extents%xe = decomp%yen(1) + level
          halo_extents%ys = decomp%yst(2)
          halo_extents%ye = decomp%yen(2)
          halo_extents%zs = decomp%yst(3) - level
          halo_extents%ze = decomp%yen(3) + level
       else
          halo_extents%xs = 1 - level
          halo_extents%xe = s1 + level
          halo_extents%ys = 1
          halo_extents%ye = s2
          halo_extents%zs = 1 - level
          halo_extents%ze = s3 + level
       end if
    else if (ipencil == 3) then  ! Z-pencil input
       if (global) then
          halo_extents%xs = decomp%zst(1) - level
          halo_extents%xe = decomp%zen(1) + level
          halo_extents%ys = decomp%zst(2) - level
          halo_extents%ye = decomp%zen(2) + level
          halo_extents%zs = decomp%zst(3)
          halo_extents%ze = decomp%zen(3)
       else
          halo_extents%xs = 1 - level
          halo_extents%xe = s1 + level
          halo_extents%ys = 1 - level
          halo_extents%ye = s2 + level
          halo_extents%zs = 1
          halo_extents%ze = s3
       end if
    else
       ! invalid input

       ! Set defaults to silence "uninitialised errors"
       halo_extents%xs = 1; halo_extents%xe = 1
       halo_extents%ys = 1; halo_extents%ye = 1
       halo_extents%zs = 1; halo_extents%ze = 1

       call decomp_2d_abort(__FILE__, __LINE__, 10, &
                            'Invalid data passed to update_halo')
    end if
  end function init_halo_extents

  subroutine halo_exchange_real(arr, ipencil, halo_extents, level, sizes)

    type(halo_extents_t), intent(in) :: halo_extents 
    real(mytype), dimension(halo_extents%xs:halo_extents%xe,halo_extents%ys:halo_extents%ye,halo_extents%zs:halo_extents%ze), intent(inout) :: arr
#if defined(_GPU)
    attributes(device) :: arr
#endif
    integer, intent(in) :: ipencil
    integer, intent(in) :: level
    integer, dimension(3), intent(in) :: sizes

    integer :: s1, s2, s3

    integer, parameter :: data_type = real_type

    integer :: icount, ilength, ijump
    integer :: halo12, halo21, halo31, halo32
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE, 4) :: status
    integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

    integer :: ierror

    s1 = sizes(1)
    s2 = sizes(2)
    s3 = sizes(3)
    
    !$acc enter data create(requests,neighbour)
    if (ipencil == 1) then
#include "halo_exchange_x_body.f90"
    else if (ipencil == 2) then
#include "halo_exchange_y_body.f90"
    else if (ipencil == 3) then
#include "halo_exchange_z_body.f90"
    end if  ! pencil
    !$acc exit data delete(neighbour,requests)

  end subroutine halo_exchange_real

  subroutine halo_exchange_complex(arr, ipencil, halo_extents, level, sizes)

    type(halo_extents_t), intent(in) :: halo_extents 
    complex(mytype), dimension(halo_extents%xs:halo_extents%xe,halo_extents%ys:halo_extents%ye,halo_extents%zs:halo_extents%ze), intent(inout) :: arr
#if defined(_GPU)
    attributes(device) :: arr
#endif
    integer, intent(in) :: ipencil
    integer, intent(in) :: level
    integer, dimension(3), intent(in) :: sizes

    integer :: s1, s2, s3

    integer, parameter :: data_type = complex_type

    integer :: icount, ilength, ijump
    integer :: halo12, halo21, halo31, halo32
    integer, dimension(4) :: requests
    integer, dimension(MPI_STATUS_SIZE, 4) :: status
    integer :: tag_e, tag_w, tag_n, tag_s, tag_t, tag_b

    integer :: ierror
#ifdef HALO_DEBUG
    integer :: i, j, k
#endif

    s1 = sizes(1)
    s2 = sizes(2)
    s3 = sizes(3)
    
    !$acc enter data create(requests,neighbour)
    if (ipencil == 1) then
#include "halo_exchange_x_body.f90"
    else if (ipencil == 2) then
#include "halo_exchange_y_body.f90"
    else if (ipencil == 3) then
#include "halo_exchange_z_body.f90"
    end if  ! pencil
    !$acc exit data delete(neighbour,requests)

  end subroutine halo_exchange_complex
  
end module m_halo
