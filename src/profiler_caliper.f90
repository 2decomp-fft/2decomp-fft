!
! Submodule for the caliper profiler
!
#ifdef PROFILER
submodule (decomp_2d) d2d_profiler_caliper

   use caliper_mod, only : ConfigManager, ConfigManager_new, &
                           cali_begin_region, cali_end_region

   implicit none

   type(ConfigManager), save :: mgr

   contains

   !
   ! Initialize
   !
   module subroutine decomp_profiler_init_noarg

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
   module subroutine decomp_profiler_fin_noarg
      call mgr%flush()
      call manager_error(mgr)
      call mgr%stop()
      call manager_error(mgr)
      call mgr%delete()
      decomp_profiler = decomp_profiler_none
   end subroutine decomp_profiler_fin_noarg

   !
   ! Log setup, currently only one profiling option
   !
   module subroutine decomp_profiler_log_int(io_unit)

      implicit none

      integer, intent(in) :: io_unit

      write(io_unit, *) "Caliper profiling active"
      write(io_unit, *) "   Profiling mode : runtime-report."
      write(io_unit, *) "   Profiling transpose : ", decomp_profiler_transpose
      write(io_unit, *) "   Profiling IO : ", decomp_profiler_io
      write(io_unit, *) "   Profiling FFT : ", decomp_profiler_fft
      write(io_unit, *) "   Profiling decomp_2d : ", decomp_profiler_d2d
      write(io_unit, *) ""

   end subroutine decomp_profiler_log_int

   !
   ! The setup of the profiler can be modified before init
   !
   module subroutine decomp_profiler_prep_bool(profiler_setup)

      implicit none

      logical, dimension(4), intent(in), optional :: profiler_setup

      ! Set the profiler id
      decomp_profiler = decomp_profiler_caliper

      ! Change the setup if provided
      if (present(profiler_setup)) then
         decomp_profiler_transpose = profiler_setup(1)
         decomp_profiler_io = profiler_setup(2)
         decomp_profiler_fft = profiler_setup(3)
         decomp_profiler_d2d = profiler_setup(4)
      endif

   end subroutine decomp_profiler_prep_bool

   !
   ! Start a timer
   !
   module subroutine decomp_profiler_start_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      call cali_begin_region(timer_name)

   end subroutine decomp_profiler_start_char

   !
   ! Stop a timer
   !
   module subroutine decomp_profiler_end_char(timer_name)

      implicit none

      character(len=*), intent(in) :: timer_name

      call cali_end_region(timer_name)

   end subroutine decomp_profiler_end_char

   !
   ! Not frequent, but the manager can produce an error
   !
   subroutine manager_error(mmm)

      implicit none

      ! Argument
      type(ConfigManager), intent(in) :: mmm

      ! Local variables
      logical :: ret
      character(len=:), allocatable :: errmsg

      ret = mmm%error()
      if (ret) then
         errmsg = mmm%error_msg()
         call decomp_2d_warning(1609, "Caliper manager error: "//trim(errmsg))
         deallocate(errmsg)
      endif

   end subroutine manager_error

end submodule d2d_profiler_caliper
#endif
