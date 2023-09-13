!! SPDX-License-Identifier: BSD-3-Clause
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routine to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! X-pencil real arrays
  subroutine alloc_x_freal_short(var, opt_global)

     implicit none

     real(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_freal_short

  subroutine alloc_x_dreal_short(var, opt_global)

     implicit none

     real(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_dreal_short

  subroutine alloc_x_freal(var, decomp, opt_global)

     implicit none

     real(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%xst(1):decomp%xen(1), &
                      decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_x_freal

  subroutine alloc_x_dreal(var, decomp, opt_global)

     implicit none

     real(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%xst(1):decomp%xen(1), &
                      decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_x_dreal

  ! X-pencil complex arrays
  subroutine alloc_x_fcplx_short(var, opt_global)

     implicit none

     complex(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_fcplx_short

  subroutine alloc_x_dcplx_short(var, opt_global)

     implicit none

     complex(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_dcplx_short

  subroutine alloc_x_fcplx(var, decomp, opt_global)

     implicit none

     complex(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%xst(1):decomp%xen(1), &
                      decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_x_fcplx

  subroutine alloc_x_dcplx(var, decomp, opt_global)

     implicit none

     complex(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%xst(1):decomp%xen(1), &
                      decomp%xst(2):decomp%xen(2), decomp%xst(3):decomp%xen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%xsz(1), decomp%xsz(2), decomp%xsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_x_dcplx

  ! Y-pencil real arrays
  subroutine alloc_y_freal_short(var, opt_global)

     implicit none

     real(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_freal_short

  subroutine alloc_y_dreal_short(var, opt_global)

     implicit none

     real(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_dreal_short

  subroutine alloc_y_freal(var, decomp, opt_global)

     implicit none

     real(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%yst(1):decomp%yen(1), &
                      decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_y_freal

  subroutine alloc_y_dreal(var, decomp, opt_global)

     implicit none

     real(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%yst(1):decomp%yen(1), &
                      decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_y_dreal

  ! Y-pencil complex arrays
  subroutine alloc_y_fcplx_short(var, opt_global)

     implicit none

     complex(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_fcplx_short

  subroutine alloc_y_dcplx_short(var, opt_global)

     implicit none

     complex(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_dcplx_short

  subroutine alloc_y_fcplx(var, decomp, opt_global)

     implicit none

     complex(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%yst(1):decomp%yen(1), &
                      decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_y_fcplx

  subroutine alloc_y_dcplx(var, decomp, opt_global)

     implicit none

     complex(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%yst(1):decomp%yen(1), &
                      decomp%yst(2):decomp%yen(2), decomp%yst(3):decomp%yen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%ysz(1), decomp%ysz(2), decomp%ysz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_y_dcplx

  ! Z-pencil real arrays
  subroutine alloc_z_freal_short(var, opt_global)

     implicit none

     real(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_freal_short

  subroutine alloc_z_dreal_short(var, opt_global)

     implicit none

     real(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_dreal_short

  subroutine alloc_z_freal(var, decomp, opt_global)

     implicit none

     real(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%zst(1):decomp%zen(1), &
                      decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_z_freal

  subroutine alloc_z_dreal(var, decomp, opt_global)

     implicit none

     real(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%zst(1):decomp%zen(1), &
                      decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_z_dreal

  ! Y-pencil complex arrays
  subroutine alloc_z_fcplx_short(var, opt_global)

     implicit none

     complex(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_fcplx_short

  subroutine alloc_z_dcplx_short(var, opt_global)

     implicit none

     complex(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_dcplx_short

  subroutine alloc_z_fcplx(var, decomp, opt_global)

     implicit none

     complex(kind(0._real32)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%zst(1):decomp%zen(1), &
                      decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_z_fcplx

  subroutine alloc_z_dcplx(var, decomp, opt_global)

     implicit none

     complex(kind(0._real64)), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     logical :: global
     integer :: alloc_stat, errorcode

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global) then
        allocate (var(decomp%zst(1):decomp%zen(1), &
                      decomp%zst(2):decomp%zen(2), decomp%zst(3):decomp%zen(3)), &
                  stat=alloc_stat)
     else
        allocate (var(decomp%zsz(1), decomp%zsz(2), decomp%zsz(3)), &
                  stat=alloc_stat)
     end if

     if (alloc_stat /= 0) then
        errorcode = 8
        call decomp_2d_abort(errorcode, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc_z_dcplx
