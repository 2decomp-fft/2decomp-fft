!=======================================================================
! This is part of the 2DECOMP&FFT library
!
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil)
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2012 Ning Li, the Numerical Algorithms Group (NAG)
! Copyright (C) 2021               the University of Edinburgh (UoE)
!
!=======================================================================

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

   end subroutine decomp_profiler_log_int

   !
   ! The setup of the profiler can be modified before init
   !
   ! TODO Current setup is limited to "runtime-report"
   !      Advance setup is possible using
   !         - extra optional argument of type char at prepare stage
   !         - env. variable CALI_CONFIG
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
         deallocate(errmsg)
      endif

   end subroutine manager_error

end submodule d2d_profiler_caliper
#endif
