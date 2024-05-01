!! SPDX-License-Identifier: BSD-3-Clause

! This module provides MPI-IO facilities for applications based on
! 2D decomposition.

!
! The external code writing a 3D array to a file will call decomp_2d_write_one
!
! Arguments :
!   - ipencil : pencil orientation of the input data
!   - var : 3D array (real / complex / integer / logical)
!   - varname : name of the file
!   - opt_dirname : optional, name of the folder containing the file
!   - opt_mpi_file_open_info : optional, MPI hints provided to MPI_FILE_OPEN
!   - opt_mpi_file_set_view_info : optional, MPI hints provided to MPI_FILE_SET_VIEW
!   - opt_reduce_prec : optional, write double precision real / complex arrays in single precision
!   - opt_decomp : optional, decomp_info object describing the array. decomp_main is used when this is not provided
!   - opt_nb_req : optional, MPI_REQUEST associated with the non-blocking MPI-IO operation
!   - opt_io : optional, use it to close the file when the request is completed (non-blocking MPI-IO)
!

!
! The external code reading a 3D array from a file will call decomp_2d_read_one
!
! Same arguments as decomp_2d_write_one
!

!
! The external code writing a 3D array to a file already opened will call decomp_2d_write_var
!
! Arguments :
!   - io : d2d_io_mpi object obtained when opening the file (mode=decomp_2d_write_mode or decomp_2d_append_mode)
!   - ipencil : pencil orientation of the input data
!   - var : 3D array (real / complex / integer / logical)
!   - opt_reduce_prec : optional, write double precision real / complex arrays in single precision
!   - opt_decomp : optional, decomp_info object describing the array. decomp_main is used when this is not provided
!   - opt_nb_req : optional, MPI_REQUEST associated with the non-blocking MPI-IO operation
!

!
! The external code reading a 3D array from a file already opened will call decomp_2d_read_var
!
! Same arguments as decomp_2d_write_var except for the IO object, it should be opened with mode=decomp_2d_read_mode
!

!
! The external code writing planes to a file will call decomp_2d_write_plane
!
! Arguments :
!   - ipencil : pencil orientation of the input data
!   - var : array (real / complex / integer / logical)
!   - varname : name of the file
!   - opt_nplanes : optional, number of planes stacked together (mandatory is stacked planes are provided)
!   - opt_iplane : optional, location of the plane in the 3D array (mandatory if a 3D array is provided)
!   - opt_dirname : optional, name of the folder containing the file
!   - opt_mpi_file_open_info : optional, MPI hints provided to MPI_FILE_OPEN
!   - opt_mpi_file_set_view_info : optional, MPI hints provided to MPI_FILE_SET_VIEW
!   - opt_reduce_prec : optional, write double precision real / complex arrays in single precision
!   - opt_decomp : optional, decomp_info object describing the array. decomp_main is used when this is not provided
!   - opt_nb_req : optional, MPI_REQUEST associated with the non-blocking MPI-IO operation
!   - opt_io : optional, use it to close the file when the request is completed (non-blocking MPI-IO)
!

!
! The external code reading planes from a file will call decomp_2d_read_plane
!
! Arguments :
!   - ipencil : pencil orientation of the input data
!   - var : array (real / complex / integer / logical)
!   - varname : name of the file
!   - nplanes : number of planes stacked together
!   - opt_dirname : optional, name of the folder containing the file
!   - opt_mpi_file_open_info : optional, MPI hints provided to MPI_FILE_OPEN
!   - opt_mpi_file_set_view_info : optional, MPI hints provided to MPI_FILE_SET_VIEW
!   - opt_reduce_prec : optional, write double precision real / complex arrays in single precision
!   - opt_decomp : optional, decomp_info object describing the array. decomp_main is used when this is not provided
!   - opt_nb_req : optional, MPI_REQUEST associated with the non-blocking MPI-IO operation
!   - opt_io : optional, use it to close the file when the request is completed (non-blocking MPI-IO)
!

!
! The external code writing a 1D array to a file will call decomp_2d_write_scalar
!
! Arguments :
!   - io : d2d_io_mpi object obtained when opening the file (mode=decomp_2d_write_mode or decomp_2d_append_mode)
!   - n : size of the array
!   - var : 3D array (real / complex / integer / logical)
!

!
! The external code reading a 1D array from a file will call decomp_2d_read_scalar
!
! Same arguments as decomp_2d_write_scalar except for the IO object, it should be opened with mode=decomp_2d_read_mode
!

module decomp_2d_io

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io_object_mpi
   use decomp_2d_io_utilities
   use decomp_2d_mpi
   use decomp_2d_profiler
   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

   private

   public :: decomp_2d_io_init, &
             decomp_2d_io_fin, &
             decomp_2d_write_one, &
             decomp_2d_read_one, &
             decomp_2d_write_var, &
             decomp_2d_read_var, &
             decomp_2d_write_plane, &
             decomp_2d_read_plane, &
             decomp_2d_write_scalar, &
             decomp_2d_read_scalar

   ! Generic interface to handle multiple data types

   interface decomp_2d_write_one
      module procedure write_one_freal
      module procedure write_one_fcplx
      module procedure write_one_dreal
      module procedure write_one_dcplx
      module procedure write_one_ints
      module procedure write_one_logs
   end interface decomp_2d_write_one

   interface decomp_2d_read_one
      module procedure read_one_freal
      module procedure read_one_fcplx
      module procedure read_one_dreal
      module procedure read_one_dcplx
      module procedure read_one_ints
      module procedure read_one_logs
   end interface decomp_2d_read_one

   interface decomp_2d_write_var
      module procedure write_var_freal
      module procedure write_var_fcplx
      module procedure write_var_dreal
      module procedure write_var_dcplx
      module procedure write_var_ints
      module procedure write_var_logs
   end interface decomp_2d_write_var

   interface decomp_2d_read_var
      module procedure read_var_freal
      module procedure read_var_fcplx
      module procedure read_var_dreal
      module procedure read_var_dcplx
      module procedure read_var_ints
      module procedure read_var_logs
   end interface decomp_2d_read_var

   interface decomp_2d_write_plane
      module procedure write_plane_freal
      module procedure write_plane_fcplx
      module procedure write_plane_dreal
      module procedure write_plane_dcplx
      module procedure write_plane_ints
      module procedure write_plane_logs
   end interface decomp_2d_write_plane

   interface decomp_2d_read_plane
      module procedure read_plane_freal
      module procedure read_plane_fcplx
      module procedure read_plane_dreal
      module procedure read_plane_dcplx
      module procedure read_plane_ints
      module procedure read_plane_logs
   end interface decomp_2d_read_plane

   interface decomp_2d_write_scalar
      module procedure write_scalar_freal
      module procedure write_scalar_fcplx
      module procedure write_scalar_dreal
      module procedure write_scalar_dcplx
      module procedure write_scalar_ints
      module procedure write_scalar_logs
   end interface decomp_2d_write_scalar

   interface decomp_2d_read_scalar
      module procedure read_scalar_freal
      module procedure read_scalar_fcplx
      module procedure read_scalar_dreal
      module procedure read_scalar_dcplx
      module procedure read_scalar_ints
      module procedure read_scalar_logs
   end interface decomp_2d_read_scalar

contains

   !
   ! Initialize the MPI IO module
   !
   !    The external code can specify default MPI_INFO values for
   !    calls to MPI_FILE_OPEN and MPI_FILE_SET_VIEW
   !
   subroutine decomp_2d_io_init(file_open_info, file_set_view_info)

      implicit none

      integer, intent(in), optional :: file_open_info, file_set_view_info

      if (decomp_profiler_io) call decomp_profiler_start("io_mpi_init")

      ! Initialize the MPI IO object module
      call decomp_2d_io_object_mpi_init(file_open_info, file_set_view_info)

      if (decomp_profiler_io) call decomp_profiler_end("io_mpi_init")

   end subroutine decomp_2d_io_init

   !
   ! Finalize the MPI IO module
   !
   subroutine decomp_2d_io_fin()

      implicit none

      if (decomp_profiler_io) call decomp_profiler_start("io_mpi_fin")

      ! Finalize the MPI IO object module
      call decomp_2d_io_object_mpi_fin()

      if (decomp_profiler_io) call decomp_profiler_end("io_mpi_fin")

   end subroutine decomp_2d_io_fin

   !
   ! Low-level. Using MPI-IO to write a 3D array to a file
   !
   ! Inputs
   !   - ipencil : pencil orientation of the variable
   !   - varname : name of the variable
   !   - decomp : decomp_info for the variable
   !   - opt_dirname : This is mandatory if no IO reader / writer was provided
   !   - opt_mpi_xxx_info : Hints for MPI IO operations
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   !   If opt_dirname is present, MPI IO will write to the file opt_dirname/varname
   !   Otherwise, MPI IO will write to the file varname
   !
   subroutine write_one(ipencil, varname, decomp, &
                        opt_dirname, &
                        opt_mpi_file_open_info, &
                        opt_mpi_file_set_view_info, &
                        opt_nb_req, &
                        opt_io, &
                        freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      integer, intent(IN) :: ipencil
      character(len=*), intent(in) :: varname
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      integer, optional :: opt_nb_req
      type(d2d_io_mpi), optional :: opt_io
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx
      integer, contiguous, dimension(:, :, :), intent(IN), optional :: ints
      logical, contiguous, dimension(:, :, :), intent(IN), optional :: logs

      type(d2d_io_mpi) :: io

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

      ! MPI-IO
      if (present(opt_nb_req) .and. present(opt_io)) then
         call opt_io%open(varname, decomp_2d_write_mode, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call write_var(opt_io, ipencil, decomp, &
                        opt_nb_req=opt_nb_req, &
                        freal=freal, &
                        dreal=dreal, &
                        fcplx=fcplx, &
                        dcplx=dcplx, &
                        ints=ints, &
                        logs=logs)
      else
         call io%open(varname, decomp_2d_write_mode, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call write_var(io, ipencil, decomp, &
                        freal=freal, &
                        dreal=dreal, &
                        fcplx=fcplx, &
                        dcplx=dcplx, &
                        ints=ints, &
                        logs=logs)
         call io%close()
      end if

   end subroutine write_one

   !
   ! Low-level. Using MPI-IO to read a 3D array from a file
   !
   ! Inputs
   !   - ipencil : pencil orientation of the variable
   !   - varname : name of the variable
   !   - decomp : decomp_info for the variable
   !   - opt_dirname : This is mandatory if no IO reader / writer was provided
   !   - opt_mpi_xxx_info : Hints for MPI IO operations
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   !   If opt_dirname is present, MPI IO will read from the file opt_dirname/varname
   !   Otherwise, MPI IO will read from the file varname
   !
   subroutine read_one(ipencil, varname, decomp, &
                       opt_dirname, &
                       opt_mpi_file_open_info, &
                       opt_mpi_file_set_view_info, &
                       opt_nb_req, &
                       opt_io, &
                       freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      integer, intent(IN) :: ipencil
      character(len=*), intent(in) :: varname
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      integer, optional :: opt_nb_req
      type(d2d_io_mpi), optional :: opt_io
      real(real32), contiguous, dimension(:, :, :), intent(OUT), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(OUT), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(OUT), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(OUT), optional :: dcplx
      integer, contiguous, dimension(:, :, :), intent(OUT), optional :: ints
      logical, contiguous, dimension(:, :, :), intent(OUT), optional :: logs

      type(d2d_io_mpi) :: io

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

      ! MPI-IO
      if (present(opt_nb_req) .and. present(opt_io)) then
         call opt_io%open(varname, decomp_2d_read_mode, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call read_var(opt_io, ipencil, decomp, &
                       opt_nb_req=opt_nb_req, &
                       freal=freal, &
                       dreal=dreal, &
                       fcplx=fcplx, &
                       dcplx=dcplx, &
                       ints=ints, &
                       logs=logs)
      else
         call io%open(varname, decomp_2d_read_mode, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call read_var(io, ipencil, decomp, &
                       freal=freal, &
                       dreal=dreal, &
                       fcplx=fcplx, &
                       dcplx=dcplx, &
                       ints=ints, &
                       logs=logs)
         call io%close()
      end if

   end subroutine read_one

   !
   ! Low-level. Using MPI-IO to add a 3D array inside a file already opened
   !
   ! Inputs
   !   - io : d2d_io_mpi object obtained when opening the file
   !   - ipencil : pencil orientation of the variable
   !   - decomp : decomp_info for the variable
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   subroutine write_var(io, ipencil, decomp, &
                        opt_nb_req, &
                        freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, optional :: opt_nb_req
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx
      integer, contiguous, dimension(:, :, :), intent(IN), optional :: ints
      logical, contiguous, dimension(:, :, :), intent(IN), optional :: logs

      integer, dimension(3) :: sizes, starts, subsizes

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

      ! Safety check
      if (present(opt_nb_req)) then
         if (opt_nb_req /= MPI_REQUEST_NULL) then
            call decomp_2d_warning(__FILE__, __LINE__, -1, "Provided MPI request was not finished "//io%label)
         end if
         ! Wait if needed and set opt_nb_req to MPI_REQUEST_NULL
         call decomp_2d_mpi_wait(opt_nb_req)
      end if

      ! Process the size
      call io_get_size(ipencil, decomp, sizes, starts, subsizes)

      ! MPI IO
      call read_or_write(.false., io, sizes, subsizes, starts, &
                         opt_nb_req=opt_nb_req, &
                         freal=freal, &
                         dreal=dreal, &
                         fcplx=fcplx, &
                         dcplx=dcplx, &
                         ints=ints, &
                         logs=logs)

   end subroutine write_var

   !
   ! Low-level. Using MPI-IO to read a 3D array from a file already opened
   !
   ! Inputs
   !   - io : d2d_io_mpi object obtained when opening the file
   !   - ipencil : pencil orientation of the variable
   !   - decomp : decomp_info for the variable
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   subroutine read_var(io, ipencil, decomp, &
                       opt_nb_req, &
                       freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, optional :: opt_nb_req
      real(real32), contiguous, dimension(:, :, :), intent(OUT), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(OUT), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(OUT), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(OUT), optional :: dcplx
      integer, contiguous, dimension(:, :, :), intent(OUT), optional :: ints
      logical, contiguous, dimension(:, :, :), intent(OUT), optional :: logs

      integer, dimension(3) :: sizes, starts, subsizes

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

      ! Safety check
      if (present(opt_nb_req)) then
         if (opt_nb_req /= MPI_REQUEST_NULL) then
            call decomp_2d_warning(__FILE__, __LINE__, -1, "Provided MPI request was not finished "//io%label)
         end if
         ! Wait if needed and set opt_nb_req to MPI_REQUEST_NULL
         call decomp_2d_mpi_wait(opt_nb_req)
      end if

      ! Process the size
      call io_get_size(ipencil, decomp, sizes, starts, subsizes)

      ! MPI IO
      call read_or_write(.true., io, sizes, subsizes, starts, &
                         opt_nb_req=opt_nb_req, &
                         freal=freal, &
                         dreal=dreal, &
                         fcplx=fcplx, &
                         dcplx=dcplx, &
                         ints=ints, &
                         logs=logs)

   end subroutine read_var

   !
   ! Low-level. Using MPI-IO to write planes to a file
   !
   ! Inputs
   !   - ipencil : pencil orientation of the planes
   !   - varname : name of the variable
   !   - decomp : decomp_info for the variable
   !   - nplanes : number of planes stacked in the array
   !   - opt_dirname : This is mandatory if no IO reader / writer was provided
   !   - opt_mpi_xxx_info : Hints for MPI IO operations
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   !   If opt_dirname is present, MPI IO will write to the file opt_dirname/varname
   !   Otherwise, MPI IO will write to the file varname
   !
   subroutine write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname, &
                          opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info, &
                          opt_nb_req, &
                          opt_io, &
                          freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      integer, intent(IN) :: ipencil
      character(len=*), intent(in) :: varname
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      integer, optional :: opt_nb_req
      type(d2d_io_mpi), optional :: opt_io
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx
      integer, contiguous, dimension(:, :, :), intent(IN), optional :: ints
      logical, contiguous, dimension(:, :, :), intent(IN), optional :: logs

      type(d2d_io_mpi) :: io
      integer, dimension(3) :: sizes, starts, subsizes

      ! Safety check
      if (nplanes < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, nplanes, "Error invalid value of nplanes ")
      end if
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

      ! Safety check
      if (present(opt_nb_req)) then
         if (opt_nb_req /= MPI_REQUEST_NULL) then
            call decomp_2d_warning(__FILE__, __LINE__, -1, "Provided MPI request was not finished "//io%label)
         end if
         ! Wait if needed and set opt_nb_req to MPI_REQUEST_NULL
         call decomp_2d_mpi_wait(opt_nb_req)
      end if

      ! Process the size
      call io_get_size(ipencil, decomp, sizes, starts, subsizes, nplanes)

      ! MPI-IO
      if (present(opt_nb_req) .and. present(opt_io)) then
         call opt_io%open(varname, decomp_2d_write_mode, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call read_or_write(.false., opt_io, sizes, subsizes, starts, &
                            opt_nb_req=opt_nb_req, &
                            freal=freal, &
                            dreal=dreal, &
                            fcplx=fcplx, &
                            dcplx=dcplx, &
                            ints=ints, &
                            logs=logs)
      else
         call io%open(varname, decomp_2d_write_mode, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call read_or_write(.false., io, sizes, subsizes, starts, &
                            freal=freal, &
                            dreal=dreal, &
                            fcplx=fcplx, &
                            dcplx=dcplx, &
                            ints=ints, &
                            logs=logs)
         call io%close()
      end if

   end subroutine write_plane

   !
   ! Low-level. Using MPI-IO to read planes from a file
   !
   ! Inputs
   !   - ipencil : pencil orientation of the planes
   !   - varname : name of the variable
   !   - decomp : decomp_info for the variable
   !   - nplanes : number of planes stacked in the array
   !   - opt_dirname : This is mandatory if no IO reader / writer was provided
   !   - opt_mpi_xxx_info : Hints for MPI IO operations
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   !   If opt_dirname is present, MPI IO will write read from the file opt_dirname/varname
   !   Otherwise, MPI IO will write from the file varname
   !
   subroutine read_plane(ipencil, varname, decomp, nplanes, &
                         opt_dirname, &
                         opt_mpi_file_open_info, &
                         opt_mpi_file_set_view_info, &
                         opt_nb_req, &
                         opt_io, &
                         freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      integer, intent(IN) :: ipencil
      character(len=*), intent(in) :: varname
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      integer, optional :: opt_nb_req
      type(d2d_io_mpi), optional :: opt_io
      real(real32), contiguous, dimension(:, :, :), intent(OUT), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(OUT), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(OUT), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(OUT), optional :: dcplx
      integer, contiguous, dimension(:, :, :), intent(OUT), optional :: ints
      logical, contiguous, dimension(:, :, :), intent(OUT), optional :: logs

      type(d2d_io_mpi) :: io
      integer, dimension(3) :: sizes, starts, subsizes

      ! Safety check
      if (nplanes < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, nplanes, "Error invalid value of nplanes ")
      end if
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

      ! Safety check
      if (present(opt_nb_req)) then
         if (opt_nb_req /= MPI_REQUEST_NULL) then
            call decomp_2d_warning(__FILE__, __LINE__, -1, "Provided MPI request was not finished "//io%label)
         end if
         ! Wait if needed and set opt_nb_req to MPI_REQUEST_NULL
         call decomp_2d_mpi_wait(opt_nb_req)
      end if

      ! Process the size
      call io_get_size(ipencil, decomp, sizes, starts, subsizes, nplanes)

      ! MPI-IO
      if (present(opt_nb_req) .and. present(opt_io)) then
         call opt_io%open(varname, decomp_2d_read_mode, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call read_or_write(.true., opt_io, sizes, subsizes, starts, &
                            opt_nb_req=opt_nb_req, &
                            freal=freal, &
                            dreal=dreal, &
                            fcplx=fcplx, &
                            dcplx=dcplx, &
                            ints=ints, &
                            logs=logs)
      else
         call io%open(varname, decomp_2d_read_mode, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info)
         call read_or_write(.true., io, sizes, subsizes, starts, &
                            opt_nb_req=opt_nb_req, &
                            freal=freal, &
                            dreal=dreal, &
                            fcplx=fcplx, &
                            dcplx=dcplx, &
                            ints=ints, &
                            logs=logs)
         call io%close()
      end if

   end subroutine read_plane

   !
   ! Low-level. Using MPI-IO to write a 1D array to a file already opened
   !
   ! Inputs
   !   - io : d2d_io_mpi object obtained when opening the file
   !   - n : size for the array
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   subroutine write_scalar(io, n, &
                           freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(in) :: n
      real(real32), contiguous, dimension(:), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:), intent(IN), optional :: dcplx
      integer, contiguous, dimension(:), intent(IN), optional :: ints
      logical, contiguous, dimension(:), intent(IN), optional :: logs

      integer :: size, subsize

      ! Process the size
      size = n
      if (nrank == 0) then
         subsize = n
      else
         subsize = 0
      end if

      ! MPI IO
      call read_or_write_scalar(.false., io, size, subsize, &
                                freal=freal, &
                                dreal=dreal, &
                                fcplx=fcplx, &
                                dcplx=dcplx, &
                                ints=ints, &
                                logs=logs)

   end subroutine write_scalar

   !
   ! Low-level. Using MPI-IO to read a 1D array from a file already opened
   !
   ! Inputs
   !   - io : d2d_io_mpi object obtained when opening the file
   !   - n : size for the array
   !   - freal / dreal / fcplx / dcplx / ints / logs : array
   !
   subroutine read_scalar(io, n, &
                          freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(in) :: n
      real(real32), contiguous, dimension(:), intent(OUT), optional :: freal
      real(real64), contiguous, dimension(:), intent(OUT), optional :: dreal
      complex(real32), contiguous, dimension(:), intent(OUT), optional :: fcplx
      complex(real64), contiguous, dimension(:), intent(OUT), optional :: dcplx
      integer, contiguous, dimension(:), intent(OUT), optional :: ints
      logical, contiguous, dimension(:), intent(OUT), optional :: logs

      integer :: size, subsize

      ! Process the size
      size = n
      subsize = n

      ! MPI IO
      call read_or_write_scalar(.true., io, size, subsize, &
                                freal=freal, &
                                dreal=dreal, &
                                fcplx=fcplx, &
                                dcplx=dcplx, &
                                ints=ints, &
                                logs=logs)

   end subroutine read_scalar

   !
   ! Low-level MPI-IO subroutine performing the read / write operations for 3D arrays
   !
   subroutine read_or_write(flag_read, io, sizes, subsizes, starts, &
                            opt_nb_req, &
                            freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      logical, intent(in) :: flag_read
      type(d2d_io_mpi), intent(inout) :: io
      integer, dimension(3), intent(in) :: sizes, subsizes, starts
      integer, optional :: opt_nb_req
      real(real32), contiguous, dimension(:, :, :), optional :: freal
      real(real64), contiguous, dimension(:, :, :), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), optional :: dcplx
      integer, contiguous, dimension(:, :, :), optional :: ints
      logical, contiguous, dimension(:, :, :), optional :: logs

      logical :: non_blocking
      integer :: ierror, data_type, newtype, type_bytes

      ! Safety check
      if (.not. io%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "IO reader / writer was not opened "//io%label)
      end if

      ! Safety check
      if (present(opt_nb_req)) then
         if (opt_nb_req /= MPI_REQUEST_NULL) then
            call decomp_2d_abort(__FILE__, __LINE__, -1, "Provided MPI request was not finished "//io%label)
         end if
      end if

      ! Allow non-blocking MPI IO
      if (present(opt_nb_req)) then
         non_blocking = .true.
      else
         non_blocking = .false.
      end if

      ! Select the MPI data type
      if (present(freal)) then
         data_type = MPI_REAL
      else if (present(dreal)) then
         data_type = MPI_DOUBLE_PRECISION
      else if (present(fcplx)) then
         data_type = MPI_COMPLEX
      else if (present(dcplx)) then
         data_type = MPI_DOUBLE_COMPLEX
      else if (present(ints)) then
         data_type = MPI_INTEGER
      else if (present(logs)) then
         data_type = MPI_LOGICAL
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Invalid inputs for "//io%label)
      end if

      ! Get the corresponding record size
      call MPI_TYPE_SIZE(data_type, type_bytes, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")

      ! Do the MPI IO
      call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, &
                                    MPI_ORDER_FORTRAN, data_type, newtype, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_CREATE_SUBARRAY")
      call MPI_TYPE_COMMIT(newtype, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
      call MPI_FILE_SET_VIEW(io%fh, io%disp, data_type, newtype, 'native', io%mpi_file_set_view_info, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      if (flag_read) then
         if (non_blocking) then
            if (present(freal)) then
               call MPI_FILE_IREAD_ALL(io%fh, freal, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(dreal)) then
               call MPI_FILE_IREAD_ALL(io%fh, dreal, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_IREAD_ALL(io%fh, fcplx, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_IREAD_ALL(io%fh, dcplx, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(ints)) then
               call MPI_FILE_IREAD_ALL(io%fh, ints, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(logs)) then
               call MPI_FILE_IREAD_ALL(io%fh, logs, product(subsizes), data_type, opt_nb_req, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_IREAD_ALL")
         else
            if (present(freal)) then
               call MPI_FILE_READ_ALL(io%fh, freal, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dreal)) then
               call MPI_FILE_READ_ALL(io%fh, dreal, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_READ_ALL(io%fh, fcplx, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_READ_ALL(io%fh, dcplx, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(ints)) then
               call MPI_FILE_READ_ALL(io%fh, ints, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(logs)) then
               call MPI_FILE_READ_ALL(io%fh, logs, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_READ_ALL")
            if (present(opt_nb_req)) opt_nb_req = MPI_REQUEST_NULL
         end if
      else
         if (non_blocking) then
            if (present(freal)) then
               call MPI_FILE_IWRITE_ALL(io%fh, freal, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(dreal)) then
               call MPI_FILE_IWRITE_ALL(io%fh, dreal, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_IWRITE_ALL(io%fh, fcplx, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_IWRITE_ALL(io%fh, dcplx, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(ints)) then
               call MPI_FILE_IWRITE_ALL(io%fh, ints, product(subsizes), data_type, opt_nb_req, ierror)
            else if (present(logs)) then
               call MPI_FILE_IWRITE_ALL(io%fh, logs, product(subsizes), data_type, opt_nb_req, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_IWRITE_ALL")
         else
            if (present(freal)) then
               call MPI_FILE_WRITE_ALL(io%fh, freal, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dreal)) then
               call MPI_FILE_WRITE_ALL(io%fh, dreal, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_WRITE_ALL(io%fh, fcplx, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_WRITE_ALL(io%fh, dcplx, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(ints)) then
               call MPI_FILE_WRITE_ALL(io%fh, ints, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(logs)) then
               call MPI_FILE_WRITE_ALL(io%fh, logs, product(subsizes), data_type, MPI_STATUS_IGNORE, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
            if (present(opt_nb_req)) opt_nb_req = MPI_REQUEST_NULL
         end if
      end if
      call MPI_TYPE_FREE(newtype, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")

      ! Update displacement for the next write operation
      call decomp_2d_io_update_disp(io%disp, sizes, type_bytes)

   end subroutine read_or_write

   !
   ! Low-level MPI-IO subroutine performing the read / write operations for 1D arrays
   !
   subroutine read_or_write_scalar(flag_read, io, size, subsize, &
                                   opt_nb_req, &
                                   freal, dreal, fcplx, dcplx, ints, logs)

      implicit none

      logical, intent(in) :: flag_read
      type(d2d_io_mpi), intent(inout) :: io
      integer, intent(in) :: size, subsize
      integer, optional :: opt_nb_req
      real(real32), contiguous, dimension(:), optional :: freal
      real(real64), contiguous, dimension(:), optional :: dreal
      complex(real32), contiguous, dimension(:), optional :: fcplx
      complex(real64), contiguous, dimension(:), optional :: dcplx
      integer, contiguous, dimension(:), optional :: ints
      logical, contiguous, dimension(:), optional :: logs

      logical :: non_blocking
      integer :: ierror, data_type, type_bytes

      ! Safety check
      if (.not. io%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "IO reader / writer was not opened "//io%label)
      end if

      ! Safety check
      if (present(opt_nb_req)) then
         if (opt_nb_req /= MPI_REQUEST_NULL) then
            call decomp_2d_abort(__FILE__, __LINE__, -1, "Provided MPI request was not finished "//io%label)
         end if
      end if

      ! Allow non-blocking MPI IO
      if (present(opt_nb_req)) then
         non_blocking = .true.
      else
         non_blocking = .false.
      end if

      ! Select the MPI data type
      if (present(freal)) then
         data_type = MPI_REAL
      else if (present(dreal)) then
         data_type = MPI_DOUBLE_PRECISION
      else if (present(fcplx)) then
         data_type = MPI_COMPLEX
      else if (present(dcplx)) then
         data_type = MPI_DOUBLE_COMPLEX
      else if (present(ints)) then
         data_type = MPI_INTEGER
      else if (present(logs)) then
         data_type = MPI_LOGICAL
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Invalid inputs for "//io%label)
      end if

      ! Get the corresponding record size
      call MPI_TYPE_SIZE(data_type, type_bytes, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")

      ! Do the MPI IO
      call MPI_FILE_SET_VIEW(io%fh, io%disp, data_type, data_type, 'native', io%mpi_file_set_view_info, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      if (flag_read) then
         if (non_blocking) then
            if (present(freal)) then
               call MPI_FILE_IREAD_ALL(io%fh, freal, subsize, data_type, opt_nb_req, ierror)
            else if (present(dreal)) then
               call MPI_FILE_IREAD_ALL(io%fh, dreal, subsize, data_type, opt_nb_req, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_IREAD_ALL(io%fh, fcplx, subsize, data_type, opt_nb_req, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_IREAD_ALL(io%fh, dcplx, subsize, data_type, opt_nb_req, ierror)
            else if (present(ints)) then
               call MPI_FILE_IREAD_ALL(io%fh, ints, subsize, data_type, opt_nb_req, ierror)
            else if (present(logs)) then
               call MPI_FILE_IREAD_ALL(io%fh, logs, subsize, data_type, opt_nb_req, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_IREAD_ALL")
         else
            if (present(freal)) then
               call MPI_FILE_READ_ALL(io%fh, freal, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dreal)) then
               call MPI_FILE_READ_ALL(io%fh, dreal, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_READ_ALL(io%fh, fcplx, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_READ_ALL(io%fh, dcplx, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(ints)) then
               call MPI_FILE_READ_ALL(io%fh, ints, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(logs)) then
               call MPI_FILE_READ_ALL(io%fh, logs, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_READ_ALL")
            if (present(opt_nb_req)) opt_nb_req = MPI_REQUEST_NULL
         end if
      else
         if (non_blocking) then
            if (present(freal)) then
               call MPI_FILE_IWRITE_ALL(io%fh, freal, subsize, data_type, opt_nb_req, ierror)
            else if (present(dreal)) then
               call MPI_FILE_IWRITE_ALL(io%fh, dreal, subsize, data_type, opt_nb_req, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_IWRITE_ALL(io%fh, fcplx, subsize, data_type, opt_nb_req, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_IWRITE_ALL(io%fh, dcplx, subsize, data_type, opt_nb_req, ierror)
            else if (present(ints)) then
               call MPI_FILE_IWRITE_ALL(io%fh, ints, subsize, data_type, opt_nb_req, ierror)
            else if (present(logs)) then
               call MPI_FILE_IWRITE_ALL(io%fh, logs, subsize, data_type, opt_nb_req, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_IWRITE_ALL")
         else
            if (present(freal)) then
               call MPI_FILE_WRITE_ALL(io%fh, freal, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dreal)) then
               call MPI_FILE_WRITE_ALL(io%fh, dreal, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(fcplx)) then
               call MPI_FILE_WRITE_ALL(io%fh, fcplx, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(dcplx)) then
               call MPI_FILE_WRITE_ALL(io%fh, dcplx, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(ints)) then
               call MPI_FILE_WRITE_ALL(io%fh, ints, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            else if (present(logs)) then
               call MPI_FILE_WRITE_ALL(io%fh, logs, subsize, data_type, MPI_STATUS_IGNORE, ierror)
            end if
            if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
            if (present(opt_nb_req)) opt_nb_req = MPI_REQUEST_NULL
         end if
      end if

      ! Update displacement for the next write operation
      io%disp = io%disp + &
                int(size, kind=MPI_OFFSET_KIND) * int(type_bytes, kind=MPI_OFFSET_KIND)

   end subroutine read_or_write_scalar

   !
   !
   !
   ! The code below was generated with the script "write_one.f90" located in the folder scripts
   !
   !
   !
   subroutine write_one_freal(ipencil, var, varname, &
                              opt_dirname, &
                              opt_mpi_file_open_info, &
                              opt_mpi_file_set_view_info, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req, &
                              opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_one(ipencil, varname, decomp, &
                     opt_dirname=opt_dirname, &
                     opt_mpi_file_open_info=opt_mpi_file_open_info, &
                     opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                     opt_nb_req=opt_nb_req, &
                     opt_io=opt_io, &
                     freal=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_freal
   !
   subroutine write_one_fcplx(ipencil, var, varname, &
                              opt_dirname, &
                              opt_mpi_file_open_info, &
                              opt_mpi_file_set_view_info, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req, &
                              opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_one(ipencil, varname, decomp, &
                     opt_dirname=opt_dirname, &
                     opt_mpi_file_open_info=opt_mpi_file_open_info, &
                     opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                     opt_nb_req=opt_nb_req, &
                     opt_io=opt_io, &
                     fcplx=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_fcplx
   !
   subroutine write_one_dreal(ipencil, var, varname, &
                              opt_dirname, &
                              opt_mpi_file_open_info, &
                              opt_mpi_file_set_view_info, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req, &
                              opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         call write_one(ipencil, varname, decomp, &
                        opt_dirname=opt_dirname, &
                        opt_mpi_file_open_info=opt_mpi_file_open_info, &
                        opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                        freal=real(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call write_one(ipencil, varname, decomp, &
                        opt_dirname=opt_dirname, &
                        opt_mpi_file_open_info=opt_mpi_file_open_info, &
                        opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                        opt_nb_req=opt_nb_req, &
                        opt_io=opt_io, &
                        dreal=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_dreal
   !
   subroutine write_one_dcplx(ipencil, var, varname, &
                              opt_dirname, &
                              opt_mpi_file_open_info, &
                              opt_mpi_file_set_view_info, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req, &
                              opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         call write_one(ipencil, varname, decomp, &
                        opt_dirname=opt_dirname, &
                        opt_mpi_file_open_info=opt_mpi_file_open_info, &
                        opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                        fcplx=cmplx(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call write_one(ipencil, varname, decomp, &
                        opt_dirname=opt_dirname, &
                        opt_mpi_file_open_info=opt_mpi_file_open_info, &
                        opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                        opt_nb_req=opt_nb_req, &
                        opt_io=opt_io, &
                        dcplx=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_dcplx
   !
   subroutine write_one_ints(ipencil, var, varname, &
                              opt_dirname, &
                              opt_mpi_file_open_info, &
                              opt_mpi_file_set_view_info, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req, &
                              opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      integer, contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_one(ipencil, varname, decomp, &
                     opt_dirname=opt_dirname, &
                     opt_mpi_file_open_info=opt_mpi_file_open_info, &
                     opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                     opt_nb_req=opt_nb_req, &
                     opt_io=opt_io, &
                     ints=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_ints
   !
   subroutine write_one_logs(ipencil, var, varname, &
                              opt_dirname, &
                              opt_mpi_file_open_info, &
                              opt_mpi_file_set_view_info, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req, &
                              opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      logical, contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_one(ipencil, varname, decomp, &
                     opt_dirname=opt_dirname, &
                     opt_mpi_file_open_info=opt_mpi_file_open_info, &
                     opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                     opt_nb_req=opt_nb_req, &
                     opt_io=opt_io, &
                     logs=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_logs
   !

   !
   !
   !
   ! The code below was generated with the script "read_one.f90" located in the folder scripts
   !
   !
   !
   subroutine read_one_freal(ipencil, var, varname, &
                             opt_dirname, &
                             opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req, &
                             opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_one(ipencil, varname, decomp, &
                    opt_dirname=opt_dirname, &
                    opt_mpi_file_open_info=opt_mpi_file_open_info, &
                    opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                    opt_nb_req=opt_nb_req, &
                    opt_io=opt_io, &
                    freal=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_freal
   !
   subroutine read_one_fcplx(ipencil, var, varname, &
                             opt_dirname, &
                             opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req, &
                             opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_one(ipencil, varname, decomp, &
                    opt_dirname=opt_dirname, &
                    opt_mpi_file_open_info=opt_mpi_file_open_info, &
                    opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                    opt_nb_req=opt_nb_req, &
                    opt_io=opt_io, &
                    fcplx=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_fcplx
   !
   subroutine read_one_dreal(ipencil, var, varname, &
                             opt_dirname, &
                             opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req, &
                             opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp
      real(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can read a file in single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then

         if (ipencil == 1) then
            call alloc_x(tmp, decomp)
         else if (ipencil == 2) then
            call alloc_y(tmp, decomp)
         else
            call alloc_z(tmp, decomp)
         end if
         call read_one(ipencil, varname, decomp, &
                       opt_dirname=opt_dirname, &
                       opt_mpi_file_open_info=opt_mpi_file_open_info, &
                       opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                       freal=tmp)
         var = tmp
         deallocate (tmp)

      else

         call read_one(ipencil, varname, decomp, &
                       opt_dirname=opt_dirname, &
                       opt_mpi_file_open_info=opt_mpi_file_open_info, &
                       opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                       opt_nb_req=opt_nb_req, &
                       opt_io=opt_io, &
                       dreal=var)

      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_dreal
   !
   subroutine read_one_dcplx(ipencil, var, varname, &
                             opt_dirname, &
                             opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req, &
                             opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp
      complex(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can read a file in single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then

         if (ipencil == 1) then
            call alloc_x(tmp, decomp)
         else if (ipencil == 2) then
            call alloc_y(tmp, decomp)
         else
            call alloc_z(tmp, decomp)
         end if
         call read_one(ipencil, varname, decomp, &
                       opt_dirname=opt_dirname, &
                       opt_mpi_file_open_info=opt_mpi_file_open_info, &
                       opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                       fcplx=tmp)
         var = tmp
         deallocate (tmp)

      else

         call read_one(ipencil, varname, decomp, &
                       opt_dirname=opt_dirname, &
                       opt_mpi_file_open_info=opt_mpi_file_open_info, &
                       opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                       opt_nb_req=opt_nb_req, &
                       opt_io=opt_io, &
                       dcplx=var)

      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_dcplx
   !
   subroutine read_one_ints(ipencil, var, varname, &
                             opt_dirname, &
                             opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req, &
                             opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      integer, contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_one(ipencil, varname, decomp, &
                    opt_dirname=opt_dirname, &
                    opt_mpi_file_open_info=opt_mpi_file_open_info, &
                    opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                    opt_nb_req=opt_nb_req, &
                    opt_io=opt_io, &
                    ints=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_ints
   !
   subroutine read_one_logs(ipencil, var, varname, &
                             opt_dirname, &
                             opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req, &
                             opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      logical, contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_one")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_one(ipencil, varname, decomp, &
                    opt_dirname=opt_dirname, &
                    opt_mpi_file_open_info=opt_mpi_file_open_info, &
                    opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                    opt_nb_req=opt_nb_req, &
                    opt_io=opt_io, &
                    logs=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_logs
   !

   !
   !
   !
   ! The code below was generated with the script "write_var.f90" located in the folder scripts
   !
   !
   !
   subroutine write_var_freal(io, ipencil, var, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_var(io, ipencil, decomp, &
                     opt_nb_req=opt_nb_req, &
                     freal=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_var")

   end subroutine write_var_freal
   !
   subroutine write_var_fcplx(io, ipencil, var, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_var(io, ipencil, decomp, &
                     opt_nb_req=opt_nb_req, &
                     fcplx=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_var")

   end subroutine write_var_fcplx
   !
   subroutine write_var_dreal(io, ipencil, var, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         call write_var(io, ipencil, decomp, &
                        freal=real(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call write_var(io, ipencil, decomp, &
                        opt_nb_req=opt_nb_req, &
                        dreal=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_var")

   end subroutine write_var_dreal
   !
   subroutine write_var_dcplx(io, ipencil, var, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         call write_var(io, ipencil, decomp, &
                        fcplx=cmplx(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call write_var(io, ipencil, decomp, &
                        opt_nb_req=opt_nb_req, &
                        dcplx=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_var")

   end subroutine write_var_dcplx
   !
   subroutine write_var_ints(io, ipencil, var, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      integer, contiguous, dimension(:, :, :), intent(IN) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_var(io, ipencil, decomp, &
                     opt_nb_req=opt_nb_req, &
                     ints=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_var")

   end subroutine write_var_ints
   !
   subroutine write_var_logs(io, ipencil, var, &
                              opt_reduce_prec, &
                              opt_decomp, &
                              opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      logical, contiguous, dimension(:, :, :), intent(IN) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_write_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call write_var(io, ipencil, decomp, &
                     opt_nb_req=opt_nb_req, &
                     logs=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_var")

   end subroutine write_var_logs
   !

   !
   !
   !
   ! The code below was generated with the script "read_var.f90" located in the folder scripts
   !
   !
   !
   subroutine read_var_freal(io, ipencil, var, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_var(io, ipencil, decomp, &
                    opt_nb_req=opt_nb_req, &
                    freal=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_var")

   end subroutine read_var_freal
   !
   subroutine read_var_fcplx(io, ipencil, var, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_var(io, ipencil, decomp, &
                    opt_nb_req=opt_nb_req, &
                    fcplx=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_var")

   end subroutine read_var_fcplx
   !
   subroutine read_var_dreal(io, ipencil, var, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp
      real(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         if (ipencil == 1) then
            call alloc_x(tmp, decomp)
         else if (ipencil == 2) then
            call alloc_y(tmp, decomp)
         else
            call alloc_z(tmp, decomp)
         end if
         call read_var(io, ipencil, decomp, &
                       opt_nb_req=opt_nb_req, &
                      freal=tmp)
         var = tmp
         deallocate(tmp)
      else
         call read_var(io, ipencil, decomp, &
                       opt_nb_req=opt_nb_req, &
                           dreal=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_var")

   end subroutine read_var_dreal
   !
   subroutine read_var_dcplx(io, ipencil, var, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp
      complex(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         if (ipencil == 1) then
            call alloc_x(tmp, decomp)
         else if (ipencil == 2) then
            call alloc_y(tmp, decomp)
         else
            call alloc_z(tmp, decomp)
         end if
         call read_var(io, ipencil, decomp, &
                       opt_nb_req=opt_nb_req, &
                      fcplx=tmp)
         var = tmp
         deallocate(tmp)
      else
         call read_var(io, ipencil, decomp, &
                       opt_nb_req=opt_nb_req, &
                           dcplx=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_var")

   end subroutine read_var_dcplx
   !
   subroutine read_var_ints(io, ipencil, var, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      integer, contiguous, dimension(:, :, :), intent(OUT) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_var(io, ipencil, decomp, &
                    opt_nb_req=opt_nb_req, &
                    ints=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_var")

   end subroutine read_var_ints
   !
   subroutine read_var_logs(io, ipencil, var, &
                             opt_reduce_prec, &
                             opt_decomp, &
                             opt_nb_req)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: ipencil
      logical, contiguous, dimension(:, :, :), intent(OUT) :: var
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req

      ! Local variable(s)
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_var")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_var(io, ipencil, decomp, &
                    opt_nb_req=opt_nb_req, &
                    logs=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_var")

   end subroutine read_var_logs
   !

   !
   !
   !
   ! The code below was generated with the script "write_plane.f90" located in the folder scripts
   !
   !
   !
   subroutine write_plane_freal(ipencil, var, varname, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_dirname, &
                                opt_mpi_file_open_info, &
                                opt_mpi_file_set_view_info, &
                                opt_reduce_prec, &
                                opt_decomp, &
                                opt_nb_req, &
                                opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      real(real32), allocatable, dimension(:, :, :) :: var2d
      TYPE(DECOMP_INFO), pointer :: decomp
      integer :: nplanes, iplane

      if (decomp_profiler_io) call decomp_profiler_start("io_write_plane")

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if
      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Error invalid value of iplane")
      end if

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      if (present(opt_nplanes)) then
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          opt_nb_req=opt_nb_req, &
                          opt_io=opt_io, &
                          freal=var)
      else
         if (ipencil == 1) then
            allocate (var2d(1, decomp%xsz(2), decomp%xsz(3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (ipencil == 2) then
            allocate (var2d(decomp%ysz(1), 1, decomp%ysz(3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (ipencil == 3) then
            allocate (var2d(decomp%zsz(1), decomp%zsz(2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          freal=var2d)
         deallocate (var2d)
      end if

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_freal
   !
   subroutine write_plane_fcplx(ipencil, var, varname, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_dirname, &
                                opt_mpi_file_open_info, &
                                opt_mpi_file_set_view_info, &
                                opt_reduce_prec, &
                                opt_decomp, &
                                opt_nb_req, &
                                opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      complex(real32), allocatable, dimension(:, :, :) :: var2d
      TYPE(DECOMP_INFO), pointer :: decomp
      integer :: nplanes, iplane

      if (decomp_profiler_io) call decomp_profiler_start("io_write_plane")

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if
      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Error invalid value of iplane")
      end if

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      if (present(opt_nplanes)) then
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          opt_nb_req=opt_nb_req, &
                          opt_io=opt_io, &
                          fcplx=var)
      else
         if (ipencil == 1) then
            allocate (var2d(1, decomp%xsz(2), decomp%xsz(3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (ipencil == 2) then
            allocate (var2d(decomp%ysz(1), 1, decomp%ysz(3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (ipencil == 3) then
            allocate (var2d(decomp%zsz(1), decomp%zsz(2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          fcplx=var2d)
         deallocate (var2d)
      end if

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_fcplx
   !
   subroutine write_plane_dreal(ipencil, var, varname, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_dirname, &
                                opt_mpi_file_open_info, &
                                opt_mpi_file_set_view_info, &
                                opt_reduce_prec, &
                                opt_decomp, &
                                opt_nb_req, &
                                opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      logical :: reduce
      real(real64), allocatable, dimension(:, :, :) :: var2d
      TYPE(DECOMP_INFO), pointer :: decomp
      integer :: nplanes, iplane

      if (decomp_profiler_io) call decomp_profiler_start("io_write_plane")

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if
      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Error invalid value of iplane")
      end if

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (present(opt_nplanes)) then
         if (reduce) then
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             freal=real(var, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             opt_nb_req=opt_nb_req, &
                             opt_io=opt_io, &
                             dreal=var)
         end if
      else
         if (ipencil == 1) then
            allocate (var2d(1, decomp%xsz(2), decomp%xsz(3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (ipencil == 2) then
            allocate (var2d(decomp%ysz(1), 1, decomp%ysz(3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (ipencil == 3) then
            allocate (var2d(decomp%zsz(1), decomp%zsz(2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         if (reduce) then
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             freal=real(var2d, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             dreal=var2d)
         end if
         deallocate (var2d)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_dreal
   !
   subroutine write_plane_dcplx(ipencil, var, varname, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_dirname, &
                                opt_mpi_file_open_info, &
                                opt_mpi_file_set_view_info, &
                                opt_reduce_prec, &
                                opt_decomp, &
                                opt_nb_req, &
                                opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      logical :: reduce
      complex(real64), allocatable, dimension(:, :, :) :: var2d
      TYPE(DECOMP_INFO), pointer :: decomp
      integer :: nplanes, iplane

      if (decomp_profiler_io) call decomp_profiler_start("io_write_plane")

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if
      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Error invalid value of iplane")
      end if

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (present(opt_nplanes)) then
         if (reduce) then
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             fcplx=cmplx(var, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             opt_nb_req=opt_nb_req, &
                             opt_io=opt_io, &
                             dcplx=var)
         end if
      else
         if (ipencil == 1) then
            allocate (var2d(1, decomp%xsz(2), decomp%xsz(3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (ipencil == 2) then
            allocate (var2d(decomp%ysz(1), 1, decomp%ysz(3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (ipencil == 3) then
            allocate (var2d(decomp%zsz(1), decomp%zsz(2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         if (reduce) then
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             fcplx=cmplx(var2d, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_dirname=opt_dirname, &
                             opt_mpi_file_open_info=opt_mpi_file_open_info, &
                             opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                             dcplx=var2d)
         end if
         deallocate (var2d)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_dcplx
   !
   subroutine write_plane_ints(ipencil, var, varname, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_dirname, &
                                opt_mpi_file_open_info, &
                                opt_mpi_file_set_view_info, &
                                opt_reduce_prec, &
                                opt_decomp, &
                                opt_nb_req, &
                                opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      integer, contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      integer, allocatable, dimension(:, :, :) :: var2d
      TYPE(DECOMP_INFO), pointer :: decomp
      integer :: nplanes, iplane

      if (decomp_profiler_io) call decomp_profiler_start("io_write_plane")

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if
      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Error invalid value of iplane")
      end if

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      if (present(opt_nplanes)) then
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          opt_nb_req=opt_nb_req, &
                          opt_io=opt_io, &
                          ints=var)
      else
         if (ipencil == 1) then
            allocate (var2d(1, decomp%xsz(2), decomp%xsz(3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (ipencil == 2) then
            allocate (var2d(decomp%ysz(1), 1, decomp%ysz(3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (ipencil == 3) then
            allocate (var2d(decomp%zsz(1), decomp%zsz(2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          ints=var2d)
         deallocate (var2d)
      end if

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_ints
   !
   subroutine write_plane_logs(ipencil, var, varname, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_dirname, &
                                opt_mpi_file_open_info, &
                                opt_mpi_file_set_view_info, &
                                opt_reduce_prec, &
                                opt_decomp, &
                                opt_nb_req, &
                                opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      logical, contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      logical, allocatable, dimension(:, :, :) :: var2d
      TYPE(DECOMP_INFO), pointer :: decomp
      integer :: nplanes, iplane

      if (decomp_profiler_io) call decomp_profiler_start("io_write_plane")

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if
      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Error invalid value of iplane")
      end if

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      if (present(opt_nplanes)) then
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          opt_nb_req=opt_nb_req, &
                          opt_io=opt_io, &
                          logs=var)
      else
         if (ipencil == 1) then
            allocate (var2d(1, decomp%xsz(2), decomp%xsz(3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (ipencil == 2) then
            allocate (var2d(decomp%ysz(1), 1, decomp%ysz(3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (ipencil == 3) then
            allocate (var2d(decomp%zsz(1), decomp%zsz(2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_dirname=opt_dirname, &
                          opt_mpi_file_open_info=opt_mpi_file_open_info, &
                          opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                          logs=var2d)
         deallocate (var2d)
      end if

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_logs
   !

   !
   !
   !
   ! The code below was generated with the script "read_plane.f90" located in the folder scripts
   !
   !
   !
   subroutine read_plane_freal(ipencil, var, varname, nplanes, &
                               opt_dirname, &
                               opt_mpi_file_open_info, &
                               opt_mpi_file_set_view_info, &
                               opt_reduce_prec, &
                               opt_decomp, &
                               opt_nb_req, &
                               opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_plane")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_plane(ipencil, varname, decomp, nplanes, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                      opt_nb_req=opt_nb_req, &
                      opt_io=opt_io, &
                      freal=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_plane")

   end subroutine read_plane_freal
   !
   subroutine read_plane_fcplx(ipencil, var, varname, nplanes, &
                               opt_dirname, &
                               opt_mpi_file_open_info, &
                               opt_mpi_file_set_view_info, &
                               opt_reduce_prec, &
                               opt_decomp, &
                               opt_nb_req, &
                               opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_plane")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_plane(ipencil, varname, decomp, nplanes, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                      opt_nb_req=opt_nb_req, &
                      opt_io=opt_io, &
                      fcplx=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_plane")

   end subroutine read_plane_fcplx
   !
   subroutine read_plane_dreal(ipencil, var, varname, nplanes, &
                               opt_dirname, &
                               opt_mpi_file_open_info, &
                               opt_mpi_file_set_view_info, &
                               opt_reduce_prec, &
                               opt_decomp, &
                               opt_nb_req, &
                               opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp
      real(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_plane")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then

         if (ipencil == 1) then
            allocate (tmp(nplanes, decomp%xsz(2), decomp%xsz(3)))
         else if (ipencil == 2) then
            allocate (tmp(decomp%ysz(1), nplanes, decomp%ysz(3)))
         else if (ipencil == 3) then
            allocate (tmp(decomp%zsz(1), decomp%zsz(2), nplanes))
         end if
         call read_plane(ipencil, varname, decomp, nplanes, &
                         opt_dirname=opt_dirname, &
                         opt_mpi_file_open_info=opt_mpi_file_open_info, &
                         opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                         freal=tmp)
         var = tmp
         deallocate (tmp)
      else
         call read_plane(ipencil, varname, decomp, nplanes, &
                         opt_dirname=opt_dirname, &
                         opt_mpi_file_open_info=opt_mpi_file_open_info, &
                         opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                         opt_nb_req=opt_nb_req, &
                         opt_io=opt_io, &
                         dreal=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_plane")

   end subroutine read_plane_dreal
   !
   subroutine read_plane_dcplx(ipencil, var, varname, nplanes, &
                               opt_dirname, &
                               opt_mpi_file_open_info, &
                               opt_mpi_file_set_view_info, &
                               opt_reduce_prec, &
                               opt_decomp, &
                               opt_nb_req, &
                               opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      logical :: reduce
      TYPE(DECOMP_INFO), pointer :: decomp
      complex(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_plane")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then

         if (ipencil == 1) then
            allocate (tmp(nplanes, decomp%xsz(2), decomp%xsz(3)))
         else if (ipencil == 2) then
            allocate (tmp(decomp%ysz(1), nplanes, decomp%ysz(3)))
         else if (ipencil == 3) then
            allocate (tmp(decomp%zsz(1), decomp%zsz(2), nplanes))
         end if
         call read_plane(ipencil, varname, decomp, nplanes, &
                         opt_dirname=opt_dirname, &
                         opt_mpi_file_open_info=opt_mpi_file_open_info, &
                         opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                         fcplx=tmp)
         var = tmp
         deallocate (tmp)
      else
         call read_plane(ipencil, varname, decomp, nplanes, &
                         opt_dirname=opt_dirname, &
                         opt_mpi_file_open_info=opt_mpi_file_open_info, &
                         opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                         opt_nb_req=opt_nb_req, &
                         opt_io=opt_io, &
                         dcplx=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_plane")

   end subroutine read_plane_dcplx
   !
   subroutine read_plane_ints(ipencil, var, varname, nplanes, &
                               opt_dirname, &
                               opt_mpi_file_open_info, &
                               opt_mpi_file_set_view_info, &
                               opt_reduce_prec, &
                               opt_decomp, &
                               opt_nb_req, &
                               opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      integer, contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_plane")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_plane(ipencil, varname, decomp, nplanes, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                      opt_nb_req=opt_nb_req, &
                      opt_io=opt_io, &
                      ints=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_plane")

   end subroutine read_plane_ints
   !
   subroutine read_plane_logs(ipencil, var, varname, nplanes, &
                               opt_dirname, &
                               opt_mpi_file_open_info, &
                               opt_mpi_file_set_view_info, &
                               opt_reduce_prec, &
                               opt_decomp, &
                               opt_nb_req, &
                               opt_io)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      logical, contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in) :: nplanes
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info
      integer, intent(in), optional :: opt_mpi_file_set_view_info
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(inout), optional :: opt_nb_req
      type(d2d_io_mpi), intent(inout), optional :: opt_io
      ! Local variables
      TYPE(DECOMP_INFO), pointer :: decomp

      if (decomp_profiler_io) call decomp_profiler_start("io_read_plane")

      ! Use the provided decomp_info or the default one
      if (present(opt_decomp)) then
         decomp => opt_decomp
      else
         decomp => decomp_main
      end if

      call read_plane(ipencil, varname, decomp, nplanes, &
                      opt_dirname=opt_dirname, &
                      opt_mpi_file_open_info=opt_mpi_file_open_info, &
                      opt_mpi_file_set_view_info=opt_mpi_file_set_view_info, &
                      opt_nb_req=opt_nb_req, &
                      opt_io=opt_io, &
                      logs=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_plane")

   end subroutine read_plane_logs
   !

   !
   !
   !
   ! The code below was generated with the script "write_scalar.f90" located in the folder scripts
   !
   !
   !
   subroutine write_scalar_freal(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      real(real32), contiguous, dimension(:), intent(IN) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_write_scalar")

      call write_scalar(io, n, &
                        freal=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_write_scalar")

   end subroutine write_scalar_freal
   !
   subroutine write_scalar_fcplx(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      complex(real32), contiguous, dimension(:), intent(IN) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_write_scalar")

      call write_scalar(io, n, &
                        fcplx=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_write_scalar")

   end subroutine write_scalar_fcplx
   !
   subroutine write_scalar_dreal(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      real(real64), contiguous, dimension(:), intent(IN) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_write_scalar")

      call write_scalar(io, n, &
                        dreal=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_write_scalar")

   end subroutine write_scalar_dreal
   !
   subroutine write_scalar_dcplx(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      complex(real64), contiguous, dimension(:), intent(IN) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_write_scalar")

      call write_scalar(io, n, &
                        dcplx=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_write_scalar")

   end subroutine write_scalar_dcplx
   !
   subroutine write_scalar_ints(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      integer, contiguous, dimension(:), intent(IN) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_write_scalar")

      call write_scalar(io, n, &
                        ints=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_write_scalar")

   end subroutine write_scalar_ints
   !
   subroutine write_scalar_logs(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      logical, contiguous, dimension(:), intent(IN) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_write_scalar")

      call write_scalar(io, n, &
                        logs=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_write_scalar")

   end subroutine write_scalar_logs
   !

   !
   !
   !
   ! The code below was generated with the script "read_scalar.f90" located in the folder scripts
   !
   !
   !
   subroutine read_scalar_freal(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      real(real32), contiguous, dimension(:), intent(OUT) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_read_scalar")

      call read_scalar(io, n, &
                       freal=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_read_scalar")

   end subroutine read_scalar_freal
   !
   subroutine read_scalar_fcplx(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      complex(real32), contiguous, dimension(:), intent(OUT) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_read_scalar")

      call read_scalar(io, n, &
                       fcplx=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_read_scalar")

   end subroutine read_scalar_fcplx
   !
   subroutine read_scalar_dreal(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      real(real64), contiguous, dimension(:), intent(OUT) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_read_scalar")

      call read_scalar(io, n, &
                       dreal=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_read_scalar")

   end subroutine read_scalar_dreal
   !
   subroutine read_scalar_dcplx(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      complex(real64), contiguous, dimension(:), intent(OUT) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_read_scalar")

      call read_scalar(io, n, &
                       dcplx=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_read_scalar")

   end subroutine read_scalar_dcplx
   !
   subroutine read_scalar_ints(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      integer, contiguous, dimension(:), intent(OUT) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_read_scalar")

      call read_scalar(io, n, &
                       ints=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_read_scalar")

   end subroutine read_scalar_ints
   !
   subroutine read_scalar_logs(io, n, var)

      implicit none

      ! Arguments
      type(d2d_io_mpi), intent(INOUT) :: io
      integer, intent(IN) :: n
      logical, contiguous, dimension(:), intent(OUT) :: var

      if (decomp_profiler_io) call decomp_profiler_start("io_read_scalar")

      call read_scalar(io, n, &
                       logs=var)
      if (decomp_profiler_io) call decomp_profiler_end("io_read_scalar")

   end subroutine read_scalar_logs
   !

end module decomp_2d_io
