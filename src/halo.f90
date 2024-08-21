!! SPDX-License-Identifier: BSD-3-Clause

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Halo cell support for neighbouring pencils to exchange data
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

     integer :: i, j, k, s1, s2, s3
     integer :: data_type

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

     integer :: i, j, k, s1, s2, s3
     integer :: data_type

     integer :: ipencil
     logical, save :: first_call_x = .true., first_call_y = .true., first_call_z = .true.

     data_type = complex_type

#include "halo_common.f90"

     return
  end subroutine update_halo_complex

  subroutine exchange_halo_x_real_short(inout, opt_xlevel)

     implicit none

     real(mytype), dimension(:, :, :), intent(INOUT) :: inout
     integer, dimension(3), optional :: opt_xlevel

     call exchange_halo_x_real(inout, decomp_main, opt_xlevel)

  end subroutine exchange_halo_x_real_short

  subroutine exchange_halo_x_real(inout, decomp, opt_xlevel)

     implicit none

     real(mytype), dimension(:, :, :), intent(INOUT) :: inout
     TYPE(DECOMP_INFO), intent(in) :: decomp
     integer, dimension(3), optional :: opt_xlevel

     integer :: level_x, level_y, level_z
     integer :: ierror
     integer :: icount, ilength, ijump
     integer :: halo12
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_n, tag_s, tag_t, tag_b
     integer :: data_type
     integer :: xs, ys, zs, ye, ze, s1, s2, s3

     data_type = real_type

#include "halo_comm_x.f90"

     return
  end subroutine exchange_halo_x_real

  subroutine exchange_halo_x_complex_short(inout, opt_xlevel)

     implicit none

     complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
     integer, dimension(3), optional :: opt_xlevel

     call exchange_halo_x_complex(inout, decomp_main, opt_xlevel)

  end subroutine exchange_halo_x_complex_short

  subroutine exchange_halo_x_complex(inout, decomp, opt_xlevel)

     implicit none

     complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
     TYPE(DECOMP_INFO), intent(in) :: decomp
     integer, dimension(3), optional :: opt_xlevel

     integer :: level_x, level_y, level_z
     integer :: ierror
     integer :: icount, ilength, ijump
     integer :: halo12
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_n, tag_s, tag_t, tag_b
     integer :: data_type
     integer :: xs, ys, zs, ye, ze, s1, s2, s3

     data_type = complex_type

#include "halo_comm_x.f90"

     return
  end subroutine exchange_halo_x_complex

  subroutine exchange_halo_y_real_short(inout, opt_ylevel)

     implicit none

     real(mytype), dimension(:, :, :), intent(INOUT) :: inout
     integer, dimension(3), optional :: opt_ylevel

     call exchange_halo_y_real(inout, decomp_main, opt_ylevel)

  end subroutine exchange_halo_y_real_short

  subroutine exchange_halo_y_real(inout, decomp, opt_ylevel)

     implicit none

     real(mytype), dimension(:, :, :), intent(INOUT) :: inout
     TYPE(DECOMP_INFO), intent(in) :: decomp
     integer, dimension(3), optional :: opt_ylevel

     integer :: level_x, level_y, level_z
     integer :: ierror
     integer :: icount, ilength, ijump
     integer :: halo21
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_e, tag_w, tag_t, tag_b
     integer :: data_type
     integer :: xs, ys, zs, xe, ze, s1, s2, s3

     data_type = real_type

#include "halo_comm_y.f90"

     return
  end subroutine exchange_halo_y_real

  subroutine exchange_halo_y_complex_short(inout, opt_ylevel)

     implicit none

     complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
     integer, dimension(3), optional :: opt_ylevel

     call exchange_halo_y_complex(inout, decomp_main, opt_ylevel)

  end subroutine exchange_halo_y_complex_short

  subroutine exchange_halo_y_complex(inout, decomp, opt_ylevel)

     implicit none

     complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
     TYPE(DECOMP_INFO), intent(in) :: decomp
     integer, dimension(3), optional :: opt_ylevel

     integer :: level_x, level_y, level_z
     integer :: ierror
     integer :: icount, ilength, ijump
     integer :: halo21
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_e, tag_w, tag_t, tag_b
     integer :: data_type
     integer :: xs, ys, zs, xe, ze, s1, s2, s3

     data_type = complex_type

#include "halo_comm_y.f90"

     return
  end subroutine exchange_halo_y_complex

  subroutine exchange_halo_z_real_short(inout, opt_zlevel)
     implicit none
     real(mytype), dimension(:, :, :), intent(INOUT) :: inout
     integer, dimension(3), optional :: opt_zlevel

     call exchange_halo_z_real(inout, decomp_main, opt_zlevel)
  end subroutine exchange_halo_z_real_short

  subroutine exchange_halo_z_real(inout, decomp, opt_zlevel)
     implicit none
     real(mytype), dimension(:, :, :), intent(INOUT) :: inout
     TYPE(DECOMP_INFO), intent(in) :: decomp
     integer, dimension(3), optional :: opt_zlevel

     integer :: level_x, level_y, level_z
     integer :: ierror
     integer :: icount, ilength, ijump
     integer :: halo31, halo32
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_e, tag_w, tag_n, tag_s
     integer :: data_type
     integer :: xs, ys, zs, xe, ye, s1, s2, s3

     data_type = real_type

#include "halo_comm_z.f90"

     return
  end subroutine exchange_halo_z_real

  subroutine exchange_halo_z_complex_short(inout, opt_zlevel)

     implicit none

     complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
     integer, dimension(3), optional :: opt_zlevel

     call exchange_halo_z_complex(inout, decomp_main, opt_zlevel)
  end subroutine exchange_halo_z_complex_short

  subroutine exchange_halo_z_complex(inout, decomp, opt_zlevel)

     implicit none

     complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
     TYPE(DECOMP_INFO), intent(in) :: decomp
     integer, dimension(3), optional :: opt_zlevel

     integer :: level_x, level_y, level_z
     integer :: ierror
     integer :: icount, ilength, ijump
     integer :: halo31, halo32
     integer, dimension(4) :: requests
     integer, dimension(MPI_STATUS_SIZE, 4) :: status
     integer :: tag_e, tag_w, tag_n, tag_s
     integer :: data_type
     integer :: xs, ys, zs, xe, ye, s1, s2, s3

     data_type = complex_type

#include "halo_comm_z.f90"

     return
  end subroutine exchange_halo_z_complex
