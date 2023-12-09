!! SPDX-License-Identifier: BSD-3-Clause

!
! Module for the caliper profiler
!
module decomp_2d_profiler

   use caliper_mod, only: ConfigManager, ConfigManager_new, &
                          cali_begin_region, cali_end_region
   use decomp_2d_constants, only: decomp_profiler_caliper, &
                                  decomp_profiler_none
   use decomp_2d_mpi, only: decomp_2d_warning

   implicit none

   type(ConfigManager), save :: mgr

   private

   ! public user routines
   public :: decomp_profiler_init, &
             decomp_profiler_fin, &
             decomp_profiler_prep, &
             decomp_profiler_log, &
             decomp_profiler_start, &
             decomp_profiler_end

   ! Generic interface to initialize the profiler
   interface decomp_profiler_init
      module procedure decomp_profiler_init_noarg
   end interface decomp_profiler_init

   ! Generic interface to finalize the profiler
   interface decomp_profiler_fin
      module procedure decomp_profiler_fin_int
   end interface decomp_profiler_fin

   ! Generic interface for the profiler to log setup
   interface decomp_profiler_log
      module procedure decomp_profiler_log_int
   end interface decomp_profiler_log

   ! Generic interface to prepare the profiler before init.
   interface decomp_profiler_prep
      module procedure decomp_profiler_prep_bool
   end interface decomp_profiler_prep

   ! Generic interface for the profiler to start a given timer
   interface decomp_profiler_start
      module procedure decomp_profiler_start_char
   end interface decomp_profiler_start

   ! Generic interface for the profiler to end a given timer
   interface decomp_profiler_end
      module procedure decomp_profiler_end_char
   end interface decomp_profiler_end

contains

   !
   ! Initialize
   !
   subroutine decomp_profiler_init_noarg()

      implicit none

      ! Create the config manager with basic reporting
      mgr = ConfigManager_new()
      call manager_error(mgr)
      call mgr%add("runtime-report")
      call manager_error(mgr)

      ! Start the manager
      call mgr%start
      call manager_error(mgr)

   end subroutine decomp_profiler_init_noarg

   !
   ! Finalize
   !
   subroutine decomp_profiler_fin_int(decomp_profiler)

      implicit none

      integer, intent(out) :: decomp_profiler

      call mgr%flush()
      call manager_error(mgr)
      call mgr%stop()
      call manager_error(mgr)
      call mgr%delete()
      decomp_profiler = decomp_profiler_none

   end subroutine decomp_profiler_fin_int

   !
   ! Log setup, currently only one profiling option
   !
   subroutine decomp_profiler_log_int(io_unit)

      implicit none

      integer, intent(in) :: io_unit

      write (io_unit, *) "Caliper profiling active"
      write (io_unit, *) "   Profiling mode : runtime-report."

   end subroutine decomp_profiler_log_int

   !
   ! The setup of the profiler can be modified before init
   !
   ! TODO Current setup is limited to "runtime-report"
   !      Advance setup is possible using
   !         - extra optional argument of type char at prepare stage
   !         - env. variable CALI_CONFIG
   subroutine decomp_profiler_prep_bool(decomp_profiler, b_transpose, b_io, b_fft, b_2d, profiler_setup)

      implicit none

      integer, intent(out) :: decomp_profiler
      logical, intent(out) :: b_transpose, b_io, b_fft, b_2d
      logical, dimension(4), intent(in), optional :: profiler_setup

      ! Set the profiler id
      decomp_profiler = decomp_profiler_caliper

      ! Change the setup if provided
      if (present(profiler_setup)) then
         b_transpose = profiler_setup(1)
         b_io = profiler_setup(2)
         b_fft = profiler_setup(3)
         b_2d = profiler_setup(4)
      end if

   end subroutine decomp_profiler_prep_bool

   !
   ! Start a timer
   !
   subroutine decomp_profiler_start_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      call cali_begin_region(timer_name)

   end subroutine decomp_profiler_start_char

   !
   ! Stop a timer
   !
   subroutine decomp_profiler_end_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      call cali_end_region(timer_name)

   end subroutine decomp_profiler_end_char

   !
   ! Not frequent, but the manager can produce an error
   !
   subroutine manager_error(manager)

      implicit none

      ! Argument
      type(ConfigManager), intent(in) :: manager

      ! Local variables
      logical :: ret
      character(len=:), allocatable :: errmsg

      ret = manager%error()
      if (ret) then
         errmsg = manager%error_msg()
         call decomp_2d_warning(1609, "Caliper manager error: "//trim(errmsg))
         deallocate (errmsg)
      end if

   end subroutine manager_error

end module decomp_2d_profiler
