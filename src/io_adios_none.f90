!! SPDX-License-Identifier: BSD-3-Clause

! Preprocessor macro to deal with unused variables
#define unused(x) associate(tmp => x); end associate

!
! Dummy module when there is no ADIOS2 IO available
!
module decomp_2d_io_object_adios

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_io_utilities, only: io_get_size
   use decomp_2d_profiler
   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

   !
   ! derived type to store info for a family of readers / writers
   !
   type, public :: d2d_io_family
   contains
      procedure :: init => d2d_io_family_init                     ! Init the family of readers / writers
      procedure :: fin => d2d_io_family_fin                       ! Clear the family of readers / writers
      procedure :: register_var => d2d_io_family_register_var     ! Register a 3D variable
      procedure :: register_plane => d2d_io_family_register_plane ! Register a 2D variable or stacked 2D variables
      procedure :: log => d2d_io_family_log                       ! Print some information about the object
   end type d2d_io_family

   ! derived type to store info for a reader / writer
   type, public :: d2d_io_adios
   contains
      procedure :: open => d2d_io_adios_open              ! Open the IO
      procedure :: start => d2d_io_adios_start            ! Start the IO
      procedure :: open_start => d2d_io_adios_open_start  ! Open and start the IO
      procedure :: end => d2d_io_adios_end                ! End the IO
      procedure :: close => d2d_io_adios_close            ! Close the IO
      procedure :: end_close => d2d_io_adios_end_close    ! End and close the IO
      procedure :: write => d2d_io_adios_write            ! Write the provided array
      procedure :: read => d2d_io_adios_read              ! Read the provided array
      procedure :: log => d2d_io_adios_log                ! Print some information about the object
   end type d2d_io_adios

   private

   public :: decomp_2d_adios_object_init, &
             decomp_2d_adios_object_fin, &
             decomp_2d_adios_get_default_family

contains

   subroutine d2d_io_family_init(family, label)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: label

      unused(family)
      unused(label)

   end subroutine d2d_io_family_init

   subroutine d2d_io_family_fin(family)

      implicit none

      class(d2d_io_family), intent(inout) :: family

      unused(family)

   end subroutine d2d_io_family_fin

   subroutine d2d_io_family_register_var(family, &
                                         varname, &
                                         ipencil, &
                                         type, &
                                         opt_reduce_prec, &
                                         opt_decomp)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil ! (x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(in) :: type
      logical, intent(in), optional :: opt_reduce_prec
      type(decomp_info), intent(in), optional :: opt_decomp

      unused(family)
      unused(varname)
      unused(ipencil)
      unused(type)
      unused(opt_reduce_prec)
      unused(opt_decomp)

   end subroutine d2d_io_family_register_var

   subroutine d2d_io_family_register_plane(family, &
                                           varname, &
                                           ipencil, &
                                           type, &
                                           opt_reduce_prec, &
                                           opt_decomp, &
                                           opt_nplanes)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil ! (x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(in) :: type
      logical, intent(in), optional :: opt_reduce_prec
      type(decomp_info), intent(in), optional :: opt_decomp
      integer, intent(in), optional :: opt_nplanes

      unused(family)
      unused(varname)
      unused(ipencil)
      unused(type)
      unused(opt_reduce_prec)
      unused(opt_decomp)
      unused(opt_nplanes)

   end subroutine d2d_io_family_register_plane

   subroutine register_var(family, varname, ipencil, typex, reduce_prec, decomp, nplanes)

      implicit none

      type(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil
      integer, intent(in) :: typex
      logical, intent(in) :: reduce_prec
      type(decomp_info), intent(in) :: decomp
      integer, intent(in) :: nplanes

      unused(family)
      unused(varname)
      unused(ipencil)
      unused(typex)
      unused(reduce_prec)
      unused(decomp)
      unused(nplanes)

   end subroutine register_var

   subroutine d2d_io_family_log(family, opt_iounit)

      implicit none

      class(d2d_io_family), intent(in) :: family
      integer, intent(in), optional :: opt_iounit

      unused(family)
      unused(opt_iounit)

   end subroutine d2d_io_family_log

   subroutine set_default_adios(adios_xml)

      implicit none

      character(len=*), intent(in) :: adios_xml

      unused(adios_xml)

   end subroutine set_default_adios

   subroutine d2d_io_adios_open(writer, mode, family)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer
      integer, intent(in) :: mode
      type(d2d_io_family), target, intent(inout) :: family

      unused(writer)
      unused(mode)
      unused(family)

   end subroutine d2d_io_adios_open

   subroutine d2d_io_adios_start(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      unused(writer)

   end subroutine d2d_io_adios_start

   subroutine d2d_io_adios_open_start(writer, mode, opt_family)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer
      integer, intent(in) :: mode
      type(d2d_io_family), target, intent(inout), optional :: opt_family

      unused(writer)
      unused(mode)
      unused(opt_family)

   end subroutine d2d_io_adios_open_start

   subroutine d2d_io_adios_end(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      unused(writer)

   end subroutine d2d_io_adios_end

   subroutine d2d_io_adios_close(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      unused(writer)

   end subroutine d2d_io_adios_close

   subroutine d2d_io_adios_end_close(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      unused(writer)

   end subroutine d2d_io_adios_end_close

   subroutine d2d_io_adios_write(writer, varname, mode, &
                                 freal, dreal, fcplx, dcplx)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer
      character(len=*), intent(in) :: varname
      integer, intent(in) :: mode
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx

      unused(writer)
      unused(varname)
      unused(mode)
      unused(freal)
      unused(dreal)
      unused(fcplx)
      unused(dcplx)

   end subroutine d2d_io_adios_write

   subroutine d2d_io_adios_read(writer, varname, freal, dreal, fcplx, dcplx)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer
      character(len=*), intent(in) :: varname
      real(real32), contiguous, dimension(:, :, :), intent(out), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(out), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(out), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(out), optional :: dcplx

      unused(writer)
      unused(varname)
      unused(freal)
      unused(dreal)
      unused(fcplx)
      unused(dcplx)

   end subroutine d2d_io_adios_read

   subroutine d2d_io_adios_log(writer, opt_iounit)

      implicit none

      class(d2d_io_adios), intent(in) :: writer
      integer, intent(in), optional :: opt_iounit

      unused(writer)
      unused(opt_iounit)

   end subroutine d2d_io_adios_log

   function decomp_2d_adios_get_default_family()

      implicit none

      type(d2d_io_family), pointer :: decomp_2d_adios_get_default_family

      decomp_2d_adios_get_default_family => null()

   end function decomp_2d_adios_get_default_family

   subroutine decomp_2d_adios_object_init(adios_xml)

      implicit none

      character(len=*), intent(in), optional :: adios_xml

      unused(adios_xml)

   end subroutine decomp_2d_adios_object_init

   !
   ! This should be called at the end
   !
   subroutine decomp_2d_adios_object_fin

      implicit none

   end subroutine decomp_2d_adios_object_fin

end module decomp_2d_io_object_adios

!
! Dummy module when there is no ADIOS2 IO available
!
module decomp_2d_io_adios

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io_utilities
   use decomp_2d_io_object_adios
   use decomp_2d_mpi
   use decomp_2d_profiler
   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

   ! The external code can use this variable to check if adios2 is available
   logical, parameter, public :: decomp_2d_with_adios2 = .false.

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
   ! Register a 3D variable
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

      unused(varname)
      unused(ipencil)
      unused(type)
      unused(opt_reduce_prec)
      unused(opt_decomp)

   end subroutine decomp_2d_register_var

   !
   ! Register planes
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

      unused(varname)
      unused(ipencil)
      unused(type)
      unused(opt_reduce_prec)
      unused(opt_decomp)
      unused(opt_nplanes)

   end subroutine decomp_2d_register_plane

   !
   ! Initialize the dummy ADIOS2 IO module
   !
   subroutine decomp_2d_io_adios_init(adios_xml)

      implicit none

      character(len=*), intent(in), optional :: adios_xml

      unused(adios_xml)

   end subroutine decomp_2d_io_adios_init

   !
   ! Finalize the dummy ADIOS2 IO module
   !
   subroutine decomp_2d_io_adios_fin()

      implicit none

   end subroutine decomp_2d_io_adios_fin

   !
   ! Dummy routines for reading / writing
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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

   end subroutine read_var_dcplx
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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_nplanes)
      unused(opt_iplane)
      unused(opt_reduce_prec)
      unused(opt_ipencil)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_nplanes)
      unused(opt_iplane)
      unused(opt_reduce_prec)
      unused(opt_ipencil)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_nplanes)
      unused(opt_iplane)
      unused(opt_reduce_prec)
      unused(opt_ipencil)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_mode)
      unused(opt_family)
      unused(opt_nplanes)
      unused(opt_iplane)
      unused(opt_reduce_prec)
      unused(opt_ipencil)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

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

      unused(io)
      unused(var)
      unused(varname)
      unused(opt_family)
      unused(opt_reduce_prec)

   end subroutine read_plane_dcplx

end module decomp_2d_io_adios
