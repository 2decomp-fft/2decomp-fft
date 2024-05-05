!! SPDX-License-Identifier: BSD-3-Clause

! This module provides parallel IO facilities for applications based on
! 2D decomposition.

module decomp_2d_io_adios

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io_object_adios
   use decomp_2d_io_utilities
   use decomp_2d_mpi
   use decomp_2d_profiler
   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64
#ifdef ADIOS2
   use adios2
#endif

   implicit none

   ! Default family for ADIOS2 IO
   type(d2d_io_family), save, pointer :: default_family => null()

   private

   public :: decomp_2d_io_adios_init, &
             decomp_2d_io_adios_fin, &
             decomp_2d_adios_write_var, &
             decomp_2d_adios_read_var, &
             decomp_2d_adios_write_plane, &
             decomp_2d_adios_read_plane, &
             decomp_2d_register_var, &
             decomp_2d_register_plane

   ! Generic interface to handle multiple data types

   interface decomp_2d_adios_write_var
      module procedure write_var_freal
      module procedure write_var_fcplx
      module procedure write_var_dreal
      module procedure write_var_dcplx
   end interface decomp_2d_adios_write_var

   interface decomp_2d_adios_read_var
      module procedure read_var_freal
      module procedure read_var_fcplx
      module procedure read_var_dreal
      module procedure read_var_dcplx
   end interface decomp_2d_adios_read_var

   interface decomp_2d_adios_write_plane
      module procedure write_plane_freal
      module procedure write_plane_fcplx
      module procedure write_plane_dreal
      module procedure write_plane_dcplx
   end interface decomp_2d_adios_write_plane

   interface decomp_2d_adios_read_plane
      module procedure read_plane_freal
      module procedure read_plane_fcplx
      module procedure read_plane_dreal
      module procedure read_plane_dcplx
   end interface decomp_2d_adios_read_plane

contains

   !
   ! High-level. Register a 3D variable for the default family of readers / writers.
   !
   subroutine decomp_2d_register_var(varname, &
                                     ipencil, &
                                     type, &
                                     opt_reduce_prec, &
                                     opt_decomp)

      implicit none

      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil ! (x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(in) :: type
      logical, intent(in), optional :: opt_reduce_prec
      type(decomp_info), intent(in), optional :: opt_decomp

      call default_family%register_var(varname, ipencil, type, &
                                       opt_reduce_prec=opt_reduce_prec, &
                                       opt_decomp=opt_decomp)

   end subroutine decomp_2d_register_var

   !
   ! High-level. Register planes for the default family of readers / writers.
   !
   subroutine decomp_2d_register_plane(varname, &
                                       ipencil, &
                                       type, &
                                       opt_reduce_prec, &
                                       opt_decomp, &
                                       opt_nplanes)

      implicit none

      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil ! (x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(in) :: type
      logical, intent(in), optional :: opt_reduce_prec
      type(decomp_info), intent(in), optional :: opt_decomp
      integer, intent(in), optional :: opt_nplanes

      call default_family%register_plane(varname, ipencil, type, &
                                         opt_reduce_prec=opt_reduce_prec, &
                                         opt_decomp=opt_decomp, &
                                         opt_nplanes=opt_nplanes)

   end subroutine decomp_2d_register_plane

   !
   ! Initialize the ADIOS2 IO module
   !
   subroutine decomp_2d_io_adios_init(adios_xml)

      implicit none

      character(len=*), intent(in), optional :: adios_xml

      if (decomp_profiler_io) call decomp_profiler_start("decomp_2d_io_adios_init")

      ! Initialize the ADIOS2 IO objects module
      call decomp_2d_adios_object_init(adios_xml)

      ! Get a pointer to the default IO family
      default_family => decomp_2d_adios_get_default_family()

      if (decomp_profiler_io) call decomp_profiler_end("decomp_2d_io_adios_init")

   end subroutine decomp_2d_io_adios_init

   !
   ! Finalize the ADIOS2 IO module
   !
   subroutine decomp_2d_io_adios_fin()

      implicit none

      if (decomp_profiler_io) call decomp_profiler_start("decomp_2d_io_adios_fin")

      ! Release the pointer to the default IO family
      nullify (default_family)

      ! Finalize the ADIOS2 IO objects module
      call decomp_2d_adios_object_fin()

      if (decomp_profiler_io) call decomp_profiler_end("decomp_2d_io_adios_fin")

   end subroutine decomp_2d_io_adios_fin

   !
   ! Low-level. Using ADIOS2 to write the provided array to a file
   !
   ! Inputs
   !   - io : IO reader / writer
   !   - varname : name of the variable
   !   - opt_mode : writing mode
   !   - opt_family : family of IO readers / writers
   !   - freal / dreal / fcplx / dcplx : array
   !
   subroutine adios_write(io, varname, &
                          opt_mode, &
                          opt_family, &
                          freal, dreal, fcplx, dcplx)

      implicit none

      type(d2d_io_adios), intent(inout) :: io
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(inout), optional :: opt_family
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx

      integer :: mode

      if (.not. io%is_open) then
         if (present(opt_family)) then
            call io%open(decomp_2d_write_mode, opt_family)
         else
            call io%open(decomp_2d_write_mode, default_family)
         end if
      end if
      if (.not. io%is_active) then
         call io%start()
      end if

      ! Default write mode is deferred
      if (present(opt_mode)) then
         mode = opt_mode
      else
         mode = decomp_2d_io_deferred
      end if

      ! Perform IO with ADIOS2
      call io%write(varname, mode, freal, dreal, fcplx, dcplx)

   end subroutine adios_write

   !
   ! Low-level. Using ADIOS2 to read the provided array from a file
   !
   ! Inputs
   !   - io : IO reader / writer
   !   - varname : name of the variable
   !   - opt_family : family of IO readers / writers
   !   - freal / dreal / fcplx / dcplx : array
   !
   subroutine adios_read(io, varname, &
                         opt_family, &
                         freal, dreal, fcplx, dcplx)

      implicit none

      type(d2d_io_adios), intent(inout) :: io
      character(len=*), intent(in) :: varname
      type(d2d_io_family), target, intent(inout), optional :: opt_family
      real(real32), contiguous, dimension(:, :, :), intent(out), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(out), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(out), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(out), optional :: dcplx

      if (.not. io%is_open) then
         if (present(opt_family)) then
            call io%open(decomp_2d_read_mode, opt_family)
         else
            call io%open(decomp_2d_read_mode, default_family)
         end if
      end if
      if (.not. io%is_active) then
         call io%start()
      end if

      ! Perform IO with ADIOS2
      call io%read(varname, freal, dreal, fcplx, dcplx)

   end subroutine adios_read

   !
   !
   !
   ! The code below was generated with the script adios_write_var located in the folder scripts
   !
   !
   !
   subroutine write_var_freal(io, var, varname, &
                              opt_mode, &
                              opt_family, &
                              opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_var")

      call adios_write(io, varname, &
                       opt_mode=opt_mode, &
                       opt_family=opt_family, &
                       freal=var)

      associate (p => opt_reduce_prec); end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_var")

   end subroutine write_var_freal
   !
   subroutine write_var_fcplx(io, var, varname, &
                              opt_mode, &
                              opt_family, &
                              opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_var")

      call adios_write(io, varname, &
                       opt_mode=opt_mode, &
                       opt_family=opt_family, &
                       fcplx=var)

      associate (p => opt_reduce_prec); end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_var")

   end subroutine write_var_fcplx
   !
   subroutine write_var_dreal(io, var, varname, &
                              opt_mode, &
                              opt_family, &
                              opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      ! Local variable(s)
      logical :: reduce

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_var")

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         call adios_write(io, varname, &
                          opt_mode=decomp_2d_io_sync, &
                          opt_family=opt_family, &
                          freal=real(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call adios_write(io, varname, &
                          opt_mode=opt_mode, &
                          opt_family=opt_family, &
                          dreal=var)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_var")

   end subroutine write_var_dreal
   !
   subroutine write_var_dcplx(io, var, varname, &
                              opt_mode, &
                              opt_family, &
                              opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      ! Local variable(s)
      logical :: reduce

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_var")

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         call adios_write(io, varname, &
                          opt_mode=decomp_2d_io_sync, &
                          opt_family=opt_family, &
                          fcplx=cmplx(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call adios_write(io, varname, &
                          opt_mode=opt_mode, &
                          opt_family=opt_family, &
                          dcplx=var)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_var")

   end subroutine write_var_dcplx
   !

   !
   !
   !
   ! The code below was generated with the script adios_read_var located in the folder scripts
   !
   !
   !
   subroutine read_var_freal(io, var, varname, &
                             opt_family, &
                             opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_var")

      call adios_read(io, varname, &
                      opt_family=opt_family, &
                      freal=var)

      associate (p => opt_reduce_prec); end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_read_var")

   end subroutine read_var_freal
   !
   subroutine read_var_fcplx(io, var, varname, &
                             opt_family, &
                             opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_var")

      call adios_read(io, varname, &
                      opt_family=opt_family, &
                      fcplx=var)

      associate (p => opt_reduce_prec); end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_read_var")

   end subroutine read_var_fcplx
   !
   subroutine read_var_dreal(io, var, varname, &
                             opt_family, &
                             opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      ! Local variable(s)
      logical :: reduce
      real(real32), dimension(:, :, :), allocatable :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_var")

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         allocate (tmp(size(var, 1), size(var, 2), size(var, 3)))
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         freal=tmp)
         var = real(tmp, kind=real64)
         deallocate (tmp)
      else
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         dreal=var)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("adios_read_var")

   end subroutine read_var_dreal
   !
   subroutine read_var_dcplx(io, var, varname, &
                             opt_family, &
                             opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      ! Local variable(s)
      logical :: reduce
      complex(real32), dimension(:, :, :), allocatable :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_var")

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then
         allocate (tmp(size(var, 1), size(var, 2), size(var, 3)))
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         fcplx=tmp)
         var = cmplx(tmp, kind=real64)
         deallocate (tmp)
      else
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         dcplx=var)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("adios_read_var")

   end subroutine read_var_dcplx
   !

   !
   !
   !
   ! The code below was generated with the script adios_write_plane located in the folder scripts
   !
   !
   !
   subroutine write_plane_freal(io, var, varname, &
                                opt_mode, &
                                opt_family, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_reduce_prec, &
                                opt_ipencil)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      logical, intent(in), optional :: opt_reduce_prec
      integer, intent(in), optional :: opt_ipencil

      ! Local variables
      real(real32), allocatable, dimension(:, :, :) :: var2d
      integer :: iplane

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_plane")

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if

      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Invalid value")
      end if
      if (.not. (present(opt_nplanes) .or. present(opt_ipencil))) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Invalid arguments")
      end if

      if (present(opt_nplanes)) then
         call adios_write(io, varname, &
                          opt_mode=opt_mode, &
                          opt_family=opt_family, &
                          freal=var)
      else
         if (opt_ipencil == 1) then
            allocate (var2d(1, size(var, 2), size(var, 3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (opt_ipencil == 2) then
            allocate (var2d(size(var, 1), 1, size(var, 3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (opt_ipencil == 3) then
            allocate (var2d(size(var, 1), size(var, 2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         call adios_write(io, varname, &
                          opt_mode=decomp_2d_io_sync, &
                          opt_family=opt_family, &
                          freal=var2d)
         deallocate (var2d)
      end if

      associate (p => opt_reduce_prec); end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_plane")

   end subroutine write_plane_freal
   !
   subroutine write_plane_fcplx(io, var, varname, &
                                opt_mode, &
                                opt_family, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_reduce_prec, &
                                opt_ipencil)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      logical, intent(in), optional :: opt_reduce_prec
      integer, intent(in), optional :: opt_ipencil

      ! Local variables
      complex(real32), allocatable, dimension(:, :, :) :: var2d
      integer :: iplane

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_plane")

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if

      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Invalid value")
      end if
      if (.not. (present(opt_nplanes) .or. present(opt_ipencil))) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Invalid arguments")
      end if

      if (present(opt_nplanes)) then
         call adios_write(io, varname, &
                          opt_mode=opt_mode, &
                          opt_family=opt_family, &
                          fcplx=var)
      else
         if (opt_ipencil == 1) then
            allocate (var2d(1, size(var, 2), size(var, 3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (opt_ipencil == 2) then
            allocate (var2d(size(var, 1), 1, size(var, 3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (opt_ipencil == 3) then
            allocate (var2d(size(var, 1), size(var, 2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         call adios_write(io, varname, &
                          opt_mode=decomp_2d_io_sync, &
                          opt_family=opt_family, &
                          fcplx=var2d)
         deallocate (var2d)
      end if

      associate (p => opt_reduce_prec); end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_plane")

   end subroutine write_plane_fcplx
   !
   subroutine write_plane_dreal(io, var, varname, &
                                opt_mode, &
                                opt_family, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_reduce_prec, &
                                opt_ipencil)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      logical, intent(in), optional :: opt_reduce_prec
      integer, intent(in), optional :: opt_ipencil

      ! Local variables
      logical :: reduce
      real(real64), allocatable, dimension(:, :, :) :: var2d
      real(real32), allocatable, dimension(:, :, :) :: var2dbis
      integer :: iplane

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_plane")

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if

      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Invalid value")
      end if
      if (.not. (present(opt_nplanes) .or. present(opt_ipencil))) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Invalid arguments")
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (present(opt_nplanes)) then
         if (reduce) then
            allocate (var2dbis(size(var, 1), &
                               size(var, 2), &
                               size(var, 3)))
            var2dbis = real(var, kind=real32)
            call adios_write(io, varname, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             freal=var2dbis)
         else
            call adios_write(io, varname, &
                             opt_mode=opt_mode, &
                             opt_family=opt_family, &
                             dreal=var)
         end if
      else
         if (reduce .and. opt_ipencil == 1) then
            allocate (var2dbis(1, size(var, 2), size(var, 3)))
            var2dbis(1, :, :) = real(var(iplane, :, :), kind=real32)
         else if (reduce .and. opt_ipencil == 2) then
            allocate (var2dbis(size(var, 1), 1, size(var, 3)))
            var2dbis(:, 1, :) = real(var(:, iplane, :), kind=real32)
         else if (reduce .and. opt_ipencil == 3) then
            allocate (var2dbis(size(var, 1), size(var, 2), 1))
            var2dbis(:, :, 1) = real(var(:, :, iplane), kind=real32)
         else if (opt_ipencil == 1) then
            allocate (var2d(1, size(var, 2), size(var, 3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (opt_ipencil == 2) then
            allocate (var2d(size(var, 1), 1, size(var, 3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (opt_ipencil == 3) then
            allocate (var2d(size(var, 1), size(var, 2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         if (reduce) then
            call adios_write(io, varname, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             freal=var2dbis)
         else
            call adios_write(io, varname, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             dreal=var2d)
         end if
      end if
      if (allocated(var2d)) deallocate (var2d)
      if (allocated(var2dbis)) deallocate (var2dbis)

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_plane")

   end subroutine write_plane_dreal
   !
   subroutine write_plane_dcplx(io, var, varname, &
                                opt_mode, &
                                opt_family, &
                                opt_nplanes, &
                                opt_iplane, &
                                opt_reduce_prec, &
                                opt_ipencil)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), intent(inout), optional :: opt_family
      integer, intent(in), optional :: opt_nplanes
      integer, intent(in), optional :: opt_iplane
      logical, intent(in), optional :: opt_reduce_prec
      integer, intent(in), optional :: opt_ipencil

      ! Local variables
      logical :: reduce
      complex(real64), allocatable, dimension(:, :, :) :: var2d
      complex(real32), allocatable, dimension(:, :, :) :: var2dbis
      integer :: iplane

      if (decomp_profiler_io) call decomp_profiler_start("adios_write_plane")

      if (present(opt_iplane)) then
         iplane = opt_iplane
      else
         iplane = 1
      end if

      ! Safety check
      if (iplane < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, iplane, "Invalid value")
      end if
      if (.not. (present(opt_nplanes) .or. present(opt_ipencil))) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Invalid arguments")
      end if

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (present(opt_nplanes)) then
         if (reduce) then
            allocate (var2dbis(size(var, 1), &
                               size(var, 2), &
                               size(var, 3)))
            var2dbis = cmplx(var, kind=real32)
            call adios_write(io, varname, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             fcplx=var2dbis)
         else
            call adios_write(io, varname, &
                             opt_mode=opt_mode, &
                             opt_family=opt_family, &
                             dcplx=var)
         end if
      else
         if (reduce .and. opt_ipencil == 1) then
            allocate (var2dbis(1, size(var, 2), size(var, 3)))
            var2dbis(1, :, :) = cmplx(var(iplane, :, :), kind=real32)
         else if (reduce .and. opt_ipencil == 2) then
            allocate (var2dbis(size(var, 1), 1, size(var, 3)))
            var2dbis(:, 1, :) = cmplx(var(:, iplane, :), kind=real32)
         else if (reduce .and. opt_ipencil == 3) then
            allocate (var2dbis(size(var, 1), size(var, 2), 1))
            var2dbis(:, :, 1) = cmplx(var(:, :, iplane), kind=real32)
         else if (opt_ipencil == 1) then
            allocate (var2d(1, size(var, 2), size(var, 3)))
            var2d(1, :, :) = var(iplane, :, :)
         else if (opt_ipencil == 2) then
            allocate (var2d(size(var, 1), 1, size(var, 3)))
            var2d(:, 1, :) = var(:, iplane, :)
         else if (opt_ipencil == 3) then
            allocate (var2d(size(var, 1), size(var, 2), 1))
            var2d(:, :, 1) = var(:, :, iplane)
         end if
         if (reduce) then
            call adios_write(io, varname, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             fcplx=var2dbis)
         else
            call adios_write(io, varname, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             dcplx=var2d)
         end if
      end if
      if (allocated(var2d)) deallocate (var2d)
      if (allocated(var2dbis)) deallocate (var2dbis)

      if (decomp_profiler_io) call decomp_profiler_end("adios_write_plane")

   end subroutine write_plane_dcplx
   !

   !
   !
   !
   ! The code below was generated with the script adios_read_plane located in the folder scripts
   !
   !
   !
   subroutine read_plane_freal(io, var, varname, &
                               opt_family, &
                               opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_plane")

      call adios_read(io, varname, &
                      opt_family=opt_family, &
                      freal=var)
      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_read_plane")

   end subroutine read_plane_freal
   !
   subroutine read_plane_fcplx(io, var, varname, &
                               opt_family, &
                               opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_plane")

      call adios_read(io, varname, &
                      opt_family=opt_family, &
                      fcplx=var)
      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("adios_read_plane")

   end subroutine read_plane_fcplx
   !
   subroutine read_plane_dreal(io, var, varname, &
                               opt_family, &
                               opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      ! Local variables
      logical :: reduce
      real(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_plane")

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then

         allocate (tmp(size(var, 1), &
                       size(var, 2), &
                       size(var, 3)))
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         freal=tmp)
         var = real(tmp, kind=real64)
         deallocate (tmp)
      else
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         dreal=var)
      end if
      if (decomp_profiler_io) call decomp_profiler_end("adios_read_plane")

   end subroutine read_plane_dreal
   !
   subroutine read_plane_dcplx(io, var, varname, &
                               opt_family, &
                               opt_reduce_prec)

      implicit none

      ! Arguments
      type(d2d_io_adios), intent(inout) :: io
      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var
      character(len=*), intent(in) :: varname
      type(d2d_io_family), intent(inout), optional :: opt_family
      logical, intent(in), optional :: opt_reduce_prec

      ! Local variables
      logical :: reduce
      complex(real32), allocatable, dimension(:, :, :) :: tmp

      if (decomp_profiler_io) call decomp_profiler_start("adios_read_plane")

      ! One can write to single precision using opt_reduce_prec
      if (present(opt_reduce_prec)) then
         reduce = opt_reduce_prec
      else
         reduce = default_opt_reduce_prec
      end if

      if (reduce) then

         allocate (tmp(size(var, 1), &
                       size(var, 2), &
                       size(var, 3)))
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         fcplx=tmp)
         var = cmplx(tmp, kind=real64)
         deallocate (tmp)
      else
         call adios_read(io, varname, &
                         opt_family=opt_family, &
                         dcplx=var)
      end if
      if (decomp_profiler_io) call decomp_profiler_end("adios_read_plane")

   end subroutine read_plane_dcplx
   !

end module decomp_2d_io_adios
