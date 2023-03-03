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

   use mpi
   use, intrinsic :: iso_fortran_env, only: real32, real64
#if defined(_GPU) && defined(_NCCL)
   use nccl 
#endif

   implicit none

   !private 
 
#ifdef DOUBLE_PREC
   integer, parameter, public :: mytype = KIND(0._real64)
   integer, parameter, public :: real_type = MPI_DOUBLE_PRECISION
   integer, parameter, public :: real2_type = MPI_2DOUBLE_PRECISION
   integer, parameter, public :: complex_type = MPI_DOUBLE_COMPLEX
#ifdef SAVE_SINGLE
   integer, parameter, public :: mytype_single = KIND(0._real32)
   integer, parameter, public :: real_type_single = MPI_REAL
#else
   integer, parameter, public :: mytype_single = KIND(0._real64)
   integer, parameter, public :: real_type_single = MPI_DOUBLE_PRECISION
#endif
#if defined(_GPU) && defined(_NCCL)
   type(ncclDataType) :: ncclType = ncclDouble
#endif
#else
   integer, parameter, public :: mytype = KIND(0._real32)
   integer, parameter, public :: real_type = MPI_REAL
   integer, parameter, public :: real2_type = MPI_2REAL
   integer, parameter, public :: complex_type = MPI_COMPLEX
   integer, parameter, public :: mytype_single = KIND(0._real32)
   integer, parameter, public :: real_type_single = MPI_REAL
#if defined(_GPU) && defined(_NCCL)
   type(ncclDataType) :: ncclType = ncclFloat
#endif
#endif

   integer, save, public :: mytype_bytes
   ! Global MPI parameters
   integer, save, public :: nrank = -1 ! local MPI rank
   integer, save, public :: nproc = -1 ! total number of processors
   integer, save, public :: decomp_2d_comm = MPI_COMM_NULL ! MPI communicator
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

   public :: decomp_2d_mpi_init, &
             decomp_2d_mpi_fin, &
             decomp_2d_mpi_comm_free, &
             decomp_2d_abort, &
             decomp_2d_warning 

   interface decomp_2d_abort
      module procedure decomp_2d_abort_basic
      module procedure decomp_2d_abort_file_line
#if defined(_GPU) && defined(_NCCL)
      module procedure decomp_2d_abort_nccl_basic
      module procedure decomp_2d_abort_nccl_file_line
#endif
   end interface decomp_2d_abort

   interface decomp_2d_warning
      module procedure decomp_2d_warning_basic
      module procedure decomp_2d_warning_file_line
   end interface decomp_2d_warning

contains
   !
   ! Set the main MPI comunicator 
   ! 
   subroutine decomp_2d_mpi_init(comm)
      
      implicit none
      
      integer, intent(in), optional :: comm

      integer :: ierror

      ! Use the provided MPI communicator if present
      if (present(comm)) then
         decomp_2d_comm = comm
      else
         decomp_2d_comm = MPI_COMM_WORLD
      end if

      ! If the external code has not set nrank and nproc
      if (nrank == -1) then
         call MPI_COMM_RANK(decomp_2d_comm, nrank, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_RANK")
         end if
      end if
      if (nproc == -1) then
         call MPI_COMM_SIZE(decomp_2d_comm, nproc, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SIZE")
         end if
      end if
   
   end subroutine decomp_2d_mpi_init
   !
   ! Reset the MPI
   ! 
   subroutine decomp_2d_mpi_fin
      
      implicit none

      nrank = -1
      nproc = -1

   end subroutine decomp_2d_mpi_fin
   !
   ! Small wrapper to free a MPI communicator
   !
   subroutine decomp_2d_mpi_comm_free(mpi_comm)

      implicit none

      integer, intent(inout) :: mpi_comm
      integer :: ierror

      ! Return if no MPI comm to free
      if (mpi_comm == MPI_COMM_NULL) return

      ! Free the provided MPI communicator
      call MPI_COMM_FREE(mpi_comm, ierror)
      if (ierror /= 0) call decomp_2d_warning(__FILE__, __LINE__, ierror, "MPI_COMM_FREE")
      mpi_comm = MPI_COMM_NULL

   end subroutine decomp_2d_mpi_comm_free

   subroutine decomp_2d_abort_basic(errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode
      character(len=*), intent(IN) :: msg

      integer :: ierror

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT ERROR - errorcode: ', errorcode
         write (*, *) 'ERROR MESSAGE: '//msg
         write (error_unit, *) '2DECOMP&FFT ERROR - errorcode: ', errorcode
         write (error_unit, *) 'ERROR MESSAGE: '//msg
      end if
      call MPI_ABORT(decomp_2d_comm, errorcode, ierror)

   end subroutine decomp_2d_abort_basic

   subroutine decomp_2d_abort_file_line(file, line, errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode, line
      character(len=*), intent(IN) :: msg, file

      integer :: ierror

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT ERROR'
         write (*, *) '  errorcode:     ', errorcode
         write (*, *) '  error in file  '//file
         write (*, *) '           line  ', line
         write (*, *) '  error message: '//msg
         write (error_unit, *) '2DECOMP&FFT ERROR'
         write (error_unit, *) '  errorcode:     ', errorcode
         write (error_unit, *) '  error in file  '//file
         write (error_unit, *) '           line  ', line
         write (error_unit, *) '  error message: '//msg
      end if
      call MPI_ABORT(decomp_2d_comm, errorcode, ierror)

   end subroutine decomp_2d_abort_file_line

#if defined(_GPU) && defined(_NCCL)
   !
   ! This is based on the file "nccl.h" in nvhpc 22.1
   !
   function _ncclresult_to_integer(errorcode)

      implicit none

      type(ncclresult), intent(IN) :: errorcode
      integer :: _ncclresult_to_integer

      if (errorcode == ncclSuccess) then
         _ncclresult_to_integer = 0
      elseif (errorcode == ncclUnhandledCudaError) then
         _ncclresult_to_integer = 1
      elseif (errorcode == ncclSystemError) then
         _ncclresult_to_integer = 2
      elseif (errorcode == ncclInternalError) then
         _ncclresult_to_integer = 3
      elseif (errorcode == ncclInvalidArgument) then
         _ncclresult_to_integer = 4
      elseif (errorcode == ncclInvalidUsage) then
         _ncclresult_to_integer = 5
      elseif (errorcode == ncclNumResults) then
         _ncclresult_to_integer = 6
      else
         _ncclresult_to_integer = -1
         call decomp_2d_warning(__FILE__, __LINE__, _ncclresult_to_integer, &
                                "NCCL error handling needs some update")
      end if

   end function _ncclresult_to_integer
   !
   ! Small wrapper for basic NCCL errors
   !
   subroutine decomp_2d_abort_nccl_basic(errorcode, msg)

      implicit none

      type(ncclresult), intent(IN) :: errorcode
      character(len=*), intent(IN) :: msg

      call decomp_2d_abort(_ncclresult_to_integer(errorcode), &
                           msg//" "//ncclGetErrorString(errorcode))

   end subroutine decomp_2d_abort_nccl_basic

   !
   ! Small wrapper for NCCL errors
   !
   subroutine decomp_2d_abort_nccl_file_line(file, line, errorcode, msg)

      implicit none

      type(ncclresult), intent(IN) :: errorcode
      integer, intent(in) :: line
      character(len=*), intent(IN) :: msg, file

      call decomp_2d_abort(file, &
                           line, &
                           _ncclresult_to_integer(errorcode), &
                           msg//" "//ncclGetErrorString(errorcode))

   end subroutine decomp_2d_abort_nccl_file_line
#endif

   subroutine decomp_2d_warning_basic(errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode
      character(len=*), intent(IN) :: msg

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT WARNING - errorcode: ', errorcode
         write (*, *) 'ERROR MESSAGE: '//msg
         write (error_unit, *) '2DECOMP&FFT WARNING - errorcode: ', errorcode
         write (error_unit, *) 'ERROR MESSAGE: '//msg
      end if

   end subroutine decomp_2d_warning_basic

   subroutine decomp_2d_warning_file_line(file, line, errorcode, msg)

      use iso_fortran_env, only: error_unit

      implicit none

      integer, intent(IN) :: errorcode, line
      character(len=*), intent(IN) :: msg, file

      if (nrank == 0) then
         write (*, *) '2DECOMP&FFT WARNING'
         write (*, *) '  errorcode:     ', errorcode
         write (*, *) '  error in file  '//file
         write (*, *) '           line  ', line
         write (*, *) '  error message: '//msg
         write (error_unit, *) '2DECOMP&FFT WARNING'
         write (error_unit, *) '  errorcode:     ', errorcode
         write (error_unit, *) '  error in file  '//file
         write (error_unit, *) '           line  ', line
         write (error_unit, *) '  error message: '//msg
      end if

   end subroutine decomp_2d_warning_file_line


end module decomp_2d_constants

