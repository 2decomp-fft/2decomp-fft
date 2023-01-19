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

! Constants for the 2decomp&fft library

module decomp_2d_constants

   implicit none

   !
   ! Output for the log can be changed by the external code before calling decomp_2d_init
   !
   !    0 => No log output
   !    1 => Master rank log output to stdout
   !    2 => Master rank log output to the file "decomp_2d_setup.log"
   !    3 => All ranks log output to a dedicated file
   !
   ! The default value is 2 (3 for debug builds)
   !
   integer, parameter, public :: D2D_LOG_QUIET = 0
   integer, parameter, public :: D2D_LOG_STDOUT = 1
   integer, parameter, public :: D2D_LOG_TOFILE = 2
   integer, parameter, public :: D2D_LOG_TOFILE_FULL = 3

   !
   ! Debug level can be changed by the external code before calling decomp_2d_init
   !
   ! The environment variable "DECOMP_2D_DEBUG" can be used to change the debug level
   !
   ! Debug checks are performed only when the preprocessor variable DEBUG is defined
   !
   enum, bind(c)
      enumerator :: D2D_DEBUG_LEVEL_OFF = 0
      enumerator :: D2D_DEBUG_LEVEL_CRITICAL = 1
      enumerator :: D2D_DEBUG_LEVEL_ERROR = 2
      enumerator :: D2D_DEBUG_LEVEL_WARN = 3
      enumerator :: D2D_DEBUG_LEVEL_INFO = 4
      enumerator :: D2D_DEBUG_LEVEL_DEBUG = 5
      enumerator :: D2D_DEBUG_LEVEL_TRACE = 6
   end enum

   !
   ! Profiler section
   !
   ! Integer to select the profiling tool
   !    0 => no profiling, default
   !    1 => Caliper (https://github.com/LLNL/Caliper)
   !
   enum, bind(c)
      enumerator :: decomp_profiler_none = 0
      enumerator :: decomp_profiler_caliper = 1
   end enum

   !
   ! Supported FFT backends
   !
   integer, parameter, public :: D2D_FFT_BACKEND_GENERIC = 0
   integer, parameter, public :: D2D_FFT_BACKEND_FFTW3 = 1
   integer, parameter, public :: D2D_FFT_BACKEND_FFTW3_F03 = 2
   integer, parameter, public :: D2D_FFT_BACKEND_MKL = 3
   integer, parameter, public :: D2D_FFT_BACKEND_CUFFT = 4

   !
   ! Complex-to-complex FFT can be forward or backward
   !
   integer, parameter, public :: DECOMP_2D_FFT_FORWARD = -1
   integer, parameter, public :: DECOMP_2D_FFT_BACKWARD = 1

   !
   ! Input / output of the FFT are distributed along different pencils
   !
   integer, parameter, public :: PHYSICAL_IN_X = 1 ! Forward is input in X, output in Z
   integer, parameter, public :: PHYSICAL_IN_Z = 3 ! Forward is input in Z, output in X

end module decomp_2d_constants

