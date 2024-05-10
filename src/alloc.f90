!! SPDX-License-Identifier: BSD-3-Clause

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Utility routine to help allocate 3D arrays
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !
  ! Allocate a 3D array
  !
  subroutine alloc(ipencil, decomp, &
                   opt_global, &
                   freal, dreal, fcplx, dcplx, ints, logs)

     implicit none

     integer, intent(IN) :: ipencil
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global
     real(real32), allocatable, dimension(:, :, :), intent(INOUT), optional :: freal
     real(real64), allocatable, dimension(:, :, :), intent(INOUT), optional :: dreal
     complex(real32), allocatable, dimension(:, :, :), intent(INOUT), optional :: fcplx
     complex(real64), allocatable, dimension(:, :, :), intent(INOUT), optional :: dcplx
     integer, allocatable, dimension(:, :, :), intent(INOUT), optional :: ints
     logical, allocatable, dimension(:, :, :), intent(INOUT), optional :: logs

     logical :: global
     integer :: ierror
     integer :: st1, st2, st3
     integer :: en1, en2, en3

     if (present(opt_global)) then
        global = opt_global
     else
        global = default_opt_global
     end if

     if (global .and. ipencil == 1) then
        st1 = decomp%xst(1)
        st2 = decomp%xst(2)
        st3 = decomp%xst(3)
        en1 = decomp%xen(1)
        en2 = decomp%xen(2)
        en3 = decomp%xen(3)
     else if (global .and. ipencil == 2) then
        st1 = decomp%yst(1)
        st2 = decomp%yst(2)
        st3 = decomp%yst(3)
        en1 = decomp%yen(1)
        en2 = decomp%yen(2)
        en3 = decomp%yen(3)
     else if (global .and. ipencil == 3) then
        st1 = decomp%zst(1)
        st2 = decomp%zst(2)
        st3 = decomp%zst(3)
        en1 = decomp%zen(1)
        en2 = decomp%zen(2)
        en3 = decomp%zen(3)
     else if (ipencil == 1) then
        st1 = 1
        st2 = 1
        st3 = 1
        en1 = decomp%xsz(1)
        en2 = decomp%xsz(2)
        en3 = decomp%xsz(3)
     else if (ipencil == 2) then
        st1 = 1
        st2 = 1
        st3 = 1
        en1 = decomp%ysz(1)
        en2 = decomp%ysz(2)
        en3 = decomp%ysz(3)
     else if (ipencil == 3) then
        st1 = 1
        st2 = 1
        st3 = 1
        en1 = decomp%zsz(1)
        en2 = decomp%zsz(2)
        en3 = decomp%zsz(3)
     else
        st1 = 1
        st2 = 1
        st3 = 1
        en1 = 1
        en2 = 1
        en3 = 1
        call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Invalid pencil provided")
     end if

     if (present(freal)) then
        if (allocated(freal)) deallocate (freal)
     else if (present(dreal)) then
        if (allocated(dreal)) deallocate (dreal)
     else if (present(fcplx)) then
        if (allocated(fcplx)) deallocate (fcplx)
     else if (present(dcplx)) then
        if (allocated(dcplx)) deallocate (dcplx)
     else if (present(ints)) then
        if (allocated(ints)) deallocate (ints)
     else if (present(logs)) then
        if (allocated(logs)) deallocate (logs)
     else
        call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid arguments")
     end if

     if (present(freal)) then
        allocate (freal(st1:en1, st2:en2, st3:en3), stat=ierror)
     else if (present(dreal)) then
        allocate (dreal(st1:en1, st2:en2, st3:en3), stat=ierror)
     else if (present(fcplx)) then
        allocate (fcplx(st1:en1, st2:en2, st3:en3), stat=ierror)
     else if (present(dcplx)) then
        allocate (dcplx(st1:en1, st2:en2, st3:en3), stat=ierror)
     else if (present(ints)) then
        allocate (ints(st1:en1, st2:en2, st3:en3), stat=ierror)
     else if (present(logs)) then
        allocate (logs(st1:en1, st2:en2, st3:en3), stat=ierror)
     end if

     if (ierror /= 0) then
        call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                             'Memory allocation failed when creating new arrays')
     end if

  end subroutine alloc

  !
  ! X pencil
  !
  subroutine alloc_x_freal_short(var, opt_global)

     implicit none

     real(real32), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_freal_short

  subroutine alloc_x_dreal_short(var, opt_global)

     implicit none

     real(real64), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_dreal_short

  subroutine alloc_x_fcplx_short(var, opt_global)

     implicit none

     complex(real32), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_fcplx_short

  subroutine alloc_x_dcplx_short(var, opt_global)

     implicit none

     complex(real64), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_dcplx_short

  subroutine alloc_x_ints_short(var, opt_global)

     implicit none

     integer, allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_ints_short

  subroutine alloc_x_logs_short(var, opt_global)

     implicit none

     logical, allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_x(var, decomp_main, opt_global)

  end subroutine alloc_x_logs_short

  subroutine alloc_x_freal(var, decomp, opt_global)

     implicit none

     real(real32), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(1, decomp, opt_global, freal=var)

  end subroutine alloc_x_freal

  subroutine alloc_x_dreal(var, decomp, opt_global)

     implicit none

     real(real64), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(1, decomp, opt_global, dreal=var)

  end subroutine alloc_x_dreal

  subroutine alloc_x_fcplx(var, decomp, opt_global)

     implicit none

     complex(real32), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(1, decomp, opt_global, fcplx=var)

  end subroutine alloc_x_fcplx

  subroutine alloc_x_dcplx(var, decomp, opt_global)

     implicit none

     complex(real64), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(1, decomp, opt_global, dcplx=var)

  end subroutine alloc_x_dcplx

  subroutine alloc_x_ints(var, decomp, opt_global)

     implicit none

     integer, allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(1, decomp, opt_global, ints=var)

  end subroutine alloc_x_ints

  subroutine alloc_x_logs(var, decomp, opt_global)

     implicit none

     logical, allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(1, decomp, opt_global, logs=var)

  end subroutine alloc_x_logs

  !
  ! Y pencil
  !
  subroutine alloc_y_freal_short(var, opt_global)

     implicit none

     real(real32), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_freal_short

  subroutine alloc_y_dreal_short(var, opt_global)

     implicit none

     real(real64), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_dreal_short

  subroutine alloc_y_fcplx_short(var, opt_global)

     implicit none

     complex(real32), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_fcplx_short

  subroutine alloc_y_dcplx_short(var, opt_global)

     implicit none

     complex(real64), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_dcplx_short

  subroutine alloc_y_ints_short(var, opt_global)

     implicit none

     integer, allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_ints_short

  subroutine alloc_y_logs_short(var, opt_global)

     implicit none

     logical, allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_y(var, decomp_main, opt_global)

  end subroutine alloc_y_logs_short

  subroutine alloc_y_freal(var, decomp, opt_global)

     implicit none

     real(real32), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(2, decomp, opt_global, freal=var)

  end subroutine alloc_y_freal

  subroutine alloc_y_dreal(var, decomp, opt_global)

     implicit none

     real(real64), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(2, decomp, opt_global, dreal=var)

  end subroutine alloc_y_dreal

  subroutine alloc_y_fcplx(var, decomp, opt_global)

     implicit none

     complex(real32), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(2, decomp, opt_global, fcplx=var)

  end subroutine alloc_y_fcplx

  subroutine alloc_y_dcplx(var, decomp, opt_global)

     implicit none

     complex(real64), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(2, decomp, opt_global, dcplx=var)

  end subroutine alloc_y_dcplx

  subroutine alloc_y_ints(var, decomp, opt_global)

     implicit none

     integer, allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(2, decomp, opt_global, ints=var)

  end subroutine alloc_y_ints

  subroutine alloc_y_logs(var, decomp, opt_global)

     implicit none

     logical, allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(2, decomp, opt_global, logs=var)

  end subroutine alloc_y_logs

  !
  ! Z pencil
  !
  subroutine alloc_z_freal_short(var, opt_global)

     implicit none

     real(real32), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_freal_short

  subroutine alloc_z_dreal_short(var, opt_global)

     implicit none

     real(real64), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_dreal_short

  subroutine alloc_z_fcplx_short(var, opt_global)

     implicit none

     complex(real32), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_fcplx_short

  subroutine alloc_z_dcplx_short(var, opt_global)

     implicit none

     complex(real64), allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_dcplx_short

  subroutine alloc_z_ints_short(var, opt_global)

     implicit none

     integer, allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_ints_short

  subroutine alloc_z_logs_short(var, opt_global)

     implicit none

     logical, allocatable, dimension(:, :, :) :: var
     logical, intent(IN), optional :: opt_global

     call alloc_z(var, decomp_main, opt_global)

  end subroutine alloc_z_logs_short

  subroutine alloc_z_freal(var, decomp, opt_global)

     implicit none

     real(real32), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(3, decomp, opt_global, freal=var)

  end subroutine alloc_z_freal

  subroutine alloc_z_dreal(var, decomp, opt_global)

     implicit none

     real(real64), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(3, decomp, opt_global, dreal=var)

  end subroutine alloc_z_dreal

  subroutine alloc_z_fcplx(var, decomp, opt_global)

     implicit none

     complex(real32), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(3, decomp, opt_global, fcplx=var)

  end subroutine alloc_z_fcplx

  subroutine alloc_z_dcplx(var, decomp, opt_global)

     implicit none

     complex(real64), allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(3, decomp, opt_global, dcplx=var)

  end subroutine alloc_z_dcplx

  subroutine alloc_z_ints(var, decomp, opt_global)

     implicit none

     integer, allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(3, decomp, opt_global, ints=var)

  end subroutine alloc_z_ints

  subroutine alloc_z_logs(var, decomp, opt_global)

     implicit none

     logical, allocatable, dimension(:, :, :) :: var
     TYPE(DECOMP_INFO), intent(IN) :: decomp
     logical, intent(IN), optional :: opt_global

     call alloc(3, decomp, opt_global, logs=var)

  end subroutine alloc_z_logs
