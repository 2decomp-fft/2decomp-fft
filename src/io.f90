!! SPDX-License-Identifier: BSD-3-Clause

! This module provides parallel IO facilities for applications based on
! 2D decomposition.

module decomp_2d_io

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io_family
   use decomp_2d_io_object
   use decomp_2d_io_utilities
   use decomp_2d_mpi
   use decomp_2d_profiler
   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64
#ifdef ADIOS2
   use adios2
#endif

   implicit none

   ! Default IO family of readers / writers
   type(d2d_io_family), target, save :: default_io_family
   type(d2d_io_family), pointer, save :: default_mpi_io_family => null()

   integer, parameter :: MAX_IOH = 10 ! How many live IO things should we handle?
   character(len=*), parameter :: io_sep = "::"
#ifndef ADIOS2
   integer, dimension(MAX_IOH), save :: fh_registry
   character(len=1024), dimension(MAX_IOH), target, save :: fh_names
   integer(kind=MPI_OFFSET_KIND), dimension(MAX_IOH), save :: fh_disp
#else
   type(adios2_adios), target :: adios
   character(len=1024), dimension(MAX_IOH), target, save :: engine_names
   logical, dimension(MAX_IOH), target, save :: engine_live
   type(adios2_engine), dimension(MAX_IOH), save :: engine_registry
#endif

   private        ! Make everything private unless declared public

   public :: decomp_2d_write_one, decomp_2d_read_one, &
             decomp_2d_write_var, decomp_2d_read_var, &
             decomp_2d_write_scalar, decomp_2d_read_scalar, &
             decomp_2d_write_plane, decomp_2d_write_every, &
             decomp_2d_write_subdomain, &
             decomp_2d_write_outflow, decomp_2d_read_inflow, &
             decomp_2d_io_init, decomp_2d_io_fin, & ! XXX: initialise/finalise 2decomp&fft IO module
             decomp_2d_io_register_var3d, &
             decomp_2d_io_register_var2d, &
             gen_iodir_name

   ! Generic interface to handle multiple data types

   interface decomp_2d_write_one
      module procedure write_one_freal
      module procedure write_one_fcplx
      module procedure write_one_dreal
      module procedure write_one_dcplx
   end interface decomp_2d_write_one

   interface decomp_2d_read_one
      module procedure read_one_freal
      module procedure read_one_fcplx
      module procedure read_one_dreal
      module procedure read_one_dcplx
   end interface decomp_2d_read_one

   interface decomp_2d_write_var
      module procedure write_var_real
      module procedure write_var_complex
   end interface decomp_2d_write_var

   interface decomp_2d_read_var
      module procedure read_var_real
      module procedure read_var_complex
   end interface decomp_2d_read_var

   interface decomp_2d_write_scalar
      module procedure write_scalar_real
      module procedure write_scalar_complex
      module procedure write_scalar_integer
      module procedure write_scalar_logical
   end interface decomp_2d_write_scalar

   interface decomp_2d_read_scalar
      module procedure read_scalar_real
      module procedure read_scalar_complex
      module procedure read_scalar_integer
      module procedure read_scalar_logical
   end interface decomp_2d_read_scalar

   interface decomp_2d_write_plane
      module procedure write_plane_freal
      module procedure write_plane_fcplx
      module procedure write_plane_dreal
      module procedure write_plane_dcplx
   end interface decomp_2d_write_plane

   interface decomp_2d_write_every
      module procedure write_every_real
      module procedure write_every_complex
   end interface decomp_2d_write_every

   interface decomp_2d_write_subdomain
      module procedure write_subdomain
   end interface decomp_2d_write_subdomain

   interface decomp_2d_write_outflow
      module procedure write_outflow
   end interface decomp_2d_write_outflow

   interface decomp_2d_read_inflow
      module procedure read_inflow
   end interface decomp_2d_read_inflow

contains

   !
   ! High-level. Register one 3D variable for the default family of readers / writers.
   !
   !    See io_family.f90
   !
   subroutine decomp_2d_io_register_var3d(varname, &
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

      call default_io_family%register_var3d(varname, ipencil, type, &
                                            opt_reduce_prec=opt_reduce_prec, &
                                            opt_decomp=opt_decomp)

   end subroutine decomp_2d_io_register_var3d

   !
   ! High-level. Register planes for the default family of readers / writers.
   !
   !    See io_family.f90
   !
   subroutine decomp_2d_io_register_var2d(varname, &
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

      call default_io_family%register_var2d(varname, ipencil, type, &
                                            opt_reduce_prec=opt_reduce_prec, &
                                            opt_decomp=opt_decomp, &
                                            opt_nplanes=opt_nplanes)

   end subroutine decomp_2d_io_register_var2d

   !
   ! Initialize the IO module
   !
   subroutine decomp_2d_io_init()

#ifdef ADIOS2
      integer :: ierror
      character(len=80) :: config_file = "adios2_config.xml"
#endif

      if (decomp_profiler_io) call decomp_profiler_start("io_init")

      ! Initialize the IO family module
      call decomp_2d_io_family_init()

      ! Initialize the IO object module
      call decomp_2d_io_object_init()

#ifdef ADIOS2
      ! The ADIOS2 module is always initialized with the file "adios2_config.xml"
      ! This can be improved
      call adios2_init(adios, trim(config_file), decomp_2d_comm, ierror)
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                              "Error initialising ADIOS2 - is "//trim(config_file)//" present and valid?")
      end if
      call decomp_2d_io_family_set_default_adios(adios)
      engine_live(:) = .false.
#endif

      ! Initialize the default families of readers / writers
      call default_io_family%init("default")
#ifdef ADIOS2
      ! ADIOS2 does not support all IO operations currently
      allocate (default_mpi_io_family)
      call default_mpi_io_family%mpi_init("default_mpi")
#else
      default_mpi_io_family => default_io_family
#endif

      ! Send the default family to the IO object module
      call decomp_2d_io_object_set_default_family(default_io_family)

      if (decomp_profiler_io) call decomp_profiler_end("io_init")

   end subroutine decomp_2d_io_init

   !
   ! Finalize the IO module
   !
   subroutine decomp_2d_io_fin()

#ifdef ADIOS2
      use adios2
#endif

      implicit none

#ifdef ADIOS2
      integer :: ierror
#endif

      if (decomp_profiler_io) call decomp_profiler_start("io_fin")

#ifdef ADIOS2
      call adios2_finalize(adios, ierror)
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_finalize")
      end if
      call decomp_2d_io_family_set_default_adios()
#endif

      ! Finalize the default IO families
      call decomp_2d_io_object_set_default_family()
      call default_io_family%fin()
#ifdef ADIOS2
      call default_mpi_io_family%fin()
      deallocate (default_mpi_io_family)
#endif
      nullify (default_mpi_io_family)

      ! Finalize the IO object module
      call decomp_2d_io_object_fin()

      ! Finalize the IO family module
      call decomp_2d_io_family_fin()

      if (decomp_profiler_io) call decomp_profiler_end("io_fin")

   end subroutine decomp_2d_io_fin

   !
   ! High-level. Using MPI-IO / ADIOS2 to write a 3D array to a file
   !
   ! The code below below was generated with the script gen_io_write_one.py
   subroutine write_one_freal(ipencil, var, varname, opt_mode, opt_family, &
                              opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                     opt_mode=opt_mode, &
                     opt_family=opt_family, &
                     opt_io=opt_io, &
                     opt_dirname=opt_dirname, &
                     opt_nb_req=opt_nb_req, &
                     freal=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_freal
   !
   subroutine write_one_fcplx(ipencil, var, varname, opt_mode, opt_family, &
                              opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                     opt_mode=opt_mode, &
                     opt_family=opt_family, &
                     opt_io=opt_io, &
                     opt_dirname=opt_dirname, &
                     opt_nb_req=opt_nb_req, &
                     fcplx=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_fcplx
   !
   subroutine write_one_dreal(ipencil, var, varname, opt_mode, opt_family, &
                              opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                        opt_mode=decomp_2d_io_sync, &
                        opt_family=opt_family, &
                        opt_io=opt_io, &
                        opt_dirname=opt_dirname, &
                        opt_nb_req=opt_nb_req, &
                        freal=real(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call write_one(ipencil, varname, decomp, &
                        opt_mode=opt_mode, &
                        opt_family=opt_family, &
                        opt_io=opt_io, &
                        opt_dirname=opt_dirname, &
                        opt_nb_req=opt_nb_req, &
                        dreal=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_dreal
   !
   subroutine write_one_dcplx(ipencil, var, varname, opt_mode, opt_family, &
                              opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                        opt_mode=decomp_2d_io_sync, &
                        opt_family=opt_family, &
                        opt_io=opt_io, &
                        opt_dirname=opt_dirname, &
                        opt_nb_req=opt_nb_req, &
                        fcplx=cmplx(var, kind=real32)) ! Warning, implicit memory allocation
      else
         call write_one(ipencil, varname, decomp, &
                        opt_mode=opt_mode, &
                        opt_family=opt_family, &
                        opt_io=opt_io, &
                        opt_dirname=opt_dirname, &
                        opt_nb_req=opt_nb_req, &
                        dcplx=var)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_one")

   end subroutine write_one_dcplx

   !
   ! High-level. Using MPI-IO / ADIOS2 to read a 3D array from a file
   !
   ! The code below below was generated with the script gen_io_read_one.py
   subroutine read_one_freal(ipencil, var, varname, opt_mode, opt_family, &
                             opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      integer, intent(inout), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                    opt_mode=opt_mode, &
                    opt_family=opt_family, &
                    opt_io=opt_io, &
                    opt_dirname=opt_dirname, &
                    opt_nb_req=opt_nb_req, &
                    freal=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_freal
   !
   subroutine read_one_fcplx(ipencil, var, varname, opt_mode, opt_family, &
                             opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      integer, intent(inout), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                    opt_mode=opt_mode, &
                    opt_family=opt_family, &
                    opt_io=opt_io, &
                    opt_dirname=opt_dirname, &
                    opt_nb_req=opt_nb_req, &
                    fcplx=var)

      nullify (decomp)

      associate (p => opt_reduce_prec)
      end associate

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_fcplx
   !
   subroutine read_one_dreal(ipencil, var, varname, opt_mode, opt_family, &
                             opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      integer, intent(inout), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                       opt_mode=decomp_2d_io_sync, &
                       opt_family=opt_family, &
                       opt_io=opt_io, &
                       opt_dirname=opt_dirname, &
                       opt_nb_req=opt_nb_req, &
                       freal=tmp)
         var = tmp
         deallocate (tmp)

      else

         call read_one(ipencil, varname, decomp, &
                       opt_mode=opt_mode, &
                       opt_family=opt_family, &
                       opt_io=opt_io, &
                       opt_dirname=opt_dirname, &
                       opt_nb_req=opt_nb_req, &
                       dreal=var)

      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_dreal
   !
   subroutine read_one_dcplx(ipencil, var, varname, opt_mode, opt_family, &
                             opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(out) :: var
      character(len=*), intent(in) :: varname
      integer, intent(inout), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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
                       opt_mode=decomp_2d_io_sync, &
                       opt_family=opt_family, &
                       opt_io=opt_io, &
                       opt_dirname=opt_dirname, &
                       opt_nb_req=opt_nb_req, &
                       fcplx=tmp)
         var = tmp
         deallocate (tmp)

      else

         call read_one(ipencil, varname, decomp, &
                       opt_mode=opt_mode, &
                       opt_family=opt_family, &
                       opt_io=opt_io, &
                       opt_dirname=opt_dirname, &
                       opt_nb_req=opt_nb_req, &
                       dcplx=var)

      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_read_one")

   end subroutine read_one_dcplx

   !
   ! High-level. Using MPI-IO / ADIOS2 to write planes to a file
   !
   ! The code below below was generated with the script gen_io_write_plane.py
   subroutine write_plane_freal(ipencil, var, varname, &
                                opt_nplanes, opt_iplane, opt_mode, opt_family, &
                                opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes, opt_iplane
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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

      if (nplanes > 1) then
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_mode=opt_mode, &
                          opt_family=opt_family, &
                          opt_io=opt_io, &
                          opt_dirname=opt_dirname, &
                          opt_nb_req=opt_nb_req, &
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
                          opt_mode=decomp_2d_io_sync, &
                          opt_family=opt_family, &
                          opt_io=opt_io, &
                          opt_dirname=opt_dirname, &
                          opt_nb_req=opt_nb_req, &
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
                                opt_nplanes, opt_iplane, opt_mode, opt_family, &
                                opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes, opt_iplane
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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

      if (nplanes > 1) then
         call write_plane(ipencil, varname, decomp, nplanes, &
                          opt_mode=opt_mode, &
                          opt_family=opt_family, &
                          opt_io=opt_io, &
                          opt_dirname=opt_dirname, &
                          opt_nb_req=opt_nb_req, &
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
                          opt_mode=decomp_2d_io_sync, &
                          opt_family=opt_family, &
                          opt_io=opt_io, &
                          opt_dirname=opt_dirname, &
                          opt_nb_req=opt_nb_req, &
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
                                opt_nplanes, opt_iplane, opt_mode, opt_family, &
                                opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes, opt_iplane
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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

      if (nplanes > 1) then
         if (reduce) then
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
                             freal=real(var, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_mode=opt_mode, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
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
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
                             freal=real(var2d, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
                             dreal=var2d)
         end if
         deallocate (var2d)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_dreal
   !
   subroutine write_plane_dcplx(ipencil, var, varname, &
                                opt_nplanes, opt_iplane, opt_mode, opt_family, &
                                opt_io, opt_dirname, opt_reduce_prec, opt_decomp, opt_nb_req)

      implicit none

      ! Arguments
      integer, intent(IN) :: ipencil
      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var
      character(len=*), intent(in) :: varname
      integer, intent(in), optional :: opt_nplanes, opt_iplane
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      logical, intent(in), optional :: opt_reduce_prec
      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp
      integer, intent(out), optional :: opt_nb_req

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

      if (nplanes > 1) then
         if (reduce) then
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
                             fcplx=cmplx(var, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_mode=opt_mode, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
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
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
                             fcplx=cmplx(var2d, kind=real32)) ! Warning, implicit memory allocation
         else
            call write_plane(ipencil, varname, decomp, nplanes, &
                             opt_mode=decomp_2d_io_sync, &
                             opt_family=opt_family, &
                             opt_io=opt_io, &
                             opt_dirname=opt_dirname, &
                             opt_nb_req=opt_nb_req, &
                             dcplx=var2d)
         end if
         deallocate (var2d)
      end if

      nullify (decomp)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_plane")

   end subroutine write_plane_dcplx

   !
   ! Low-level. Using MPI-IO / ADIOS2 to write a 3D array to a file
   !
   ! Inputs
   !   - ipencil : pencil orientation of the variable
   !   - varname : name of the variable
   !   - decomp : decomp_info for the variable
   !   - opt_mode : writing mode
   !   - opt_family : IO family can be provided to avoid the default one
   !   - opt_io : IO reader / writer can be provided
   !   - opt_dirname : This is mandatory if no IO reader / writer was provided
   !   - freal / dreal / fcplx / dcplx : array
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !
   ! If opt_io is provided and open
   !    - MPI IO will write to the file handle opt_io%fh
   !    - ADIOS2 IO will write to the engine opt_io%engine
   !
   ! Otherwise
   !    - opt_dirname is mandatory
   !    - MPI IO will write to the file opt_dirname/varname
   !    - ADIOS2 IO will write to the folder opt_dirname
   !
   subroutine write_one(ipencil, varname, decomp, opt_mode, &
                        opt_family, opt_io, opt_dirname, &
                        freal, dreal, fcplx, dcplx, opt_nb_req)

      implicit none

      integer, intent(IN) :: ipencil
      character(len=*), intent(in) :: varname
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx
      integer, intent(out), optional :: opt_nb_req

      logical :: use_opt_io, opt_nb
      integer :: mode
      type(d2d_io) :: io

      ! Use opt_io only if present and open
      if (present(opt_io)) then
         if (opt_io%is_open) then
            use_opt_io = .true.
         else
            use_opt_io = .false.
         end if
      else
         use_opt_io = .false.
      end if

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if
      if ((.not. use_opt_io) .and. (.not. present(opt_dirname))) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, "Invalid arguments")
      end if

      ! Default write mode is deferred
      if (present(opt_mode)) then
         mode = opt_mode
      else
         mode = decomp_2d_io_deferred
      end if

      ! Default MPI IO for writing is blocking
      if (use_opt_io .and. mode == decomp_2d_io_deferred .and. present(opt_nb_req)) then
         opt_nb = .true.
      else
         opt_nb = .false.
      end if

      ! Default value
      if (present(opt_nb_req)) opt_nb_req = MPI_REQUEST_NULL

      if (use_opt_io) then

         if (opt_io%family%type == DECOMP_2D_IO_MPI) then

            ! MPI-IO
            call mpi_write_var(ipencil, opt_io, decomp, freal, dreal, fcplx, dcplx, opt_nb, opt_nb_req)

         else if (opt_io%family%type == DECOMP_2D_IO_ADIOS2) then

            ! ADIOS2
            call adios2_write_var(opt_io, varname, mode, freal, dreal, fcplx, dcplx)

         end if

      else

         call decomp_2d_io_object_open_and_start(io, &
                                                 opt_dirname, &
                                                 varname, &
                                                 decomp_2d_write_mode, &
                                                 opt_family=opt_family)

         if (io%family%type == DECOMP_2D_IO_MPI) then

            ! MPI-IO
            call mpi_write_var(ipencil, io, decomp, freal, dreal, fcplx, dcplx)

         else if (io%family%type == DECOMP_2D_IO_ADIOS2) then

            ! ADIOS2
            call adios2_write_var(io, varname, mode, freal, dreal, fcplx, dcplx)

         end if

         call io%end_close()

      end if

   end subroutine write_one

   !
   ! Low-level. Using MPI-IO / ADIOS2 to read a 3D array from a file
   !
   ! Inputs
   !   - ipencil : pencil orientation of the variable
   !   - varname : name of the variable
   !   - decomp : decomp_info for the variable
   !   - opt_mode : reading mode, no impact on ADIOS2
   !   - opt_family : IO family can be provided to avoid the default one
   !   - opt_io : IO reader / writer can be provided
   !   - opt_dirname : This is mandatory if no IO reader / writer was provided
   !   - freal / dreal / fcplx / dcplx : array
   !
   ! If opt_io is provided and open
   !    - MPI IO will read from the file handle opt_io%fh
   !    - ADIOS2 IO will read from the engine opt_io%engine
   !
   ! Otherwise
   !    - opt_dirname is mandatory
   !    - MPI IO will read from the file opt_dirname/varname
   !    - ADIOS2 IO will read from the folder opt_dirname
   !
   subroutine read_one(ipencil, varname, decomp, opt_mode, &
                       opt_family, opt_io, opt_dirname, &
                       freal, dreal, fcplx, dcplx, opt_nb_req)

      implicit none

      integer, intent(IN) :: ipencil
      character(len=*), intent(in) :: varname
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      real(real32), contiguous, dimension(:, :, :), intent(out), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(out), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(out), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(out), optional :: dcplx
      integer, intent(out), optional :: opt_nb_req

      logical :: use_opt_io, opt_nb
      type(d2d_io) :: io

      ! Use opt_io only if present and open
      if (present(opt_io)) then
         if (opt_io%is_open) then
            use_opt_io = .true.
         else
            use_opt_io = .false.
         end if
      else
         use_opt_io = .false.
      end if

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if
      if ((.not. use_opt_io) .and. (.not. present(opt_dirname))) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, "Invalid arguments")
      end if

      ! Default MPI IO for reading is blocking
      if (present(opt_mode)) then
         if (use_opt_io .and. opt_mode == decomp_2d_io_deferred .and. present(opt_nb_req)) then
            opt_nb = .true.
         else
            opt_nb = .false.
         end if
      else
         opt_nb = .false.
      end if

      ! Default value
      if (present(opt_nb_req)) opt_nb_req = MPI_REQUEST_NULL

      if (use_opt_io) then

         if (opt_io%family%type == DECOMP_2D_IO_MPI) then

            ! MPI-IO
            call mpi_read_var(ipencil, opt_io, decomp, freal, dreal, fcplx, dcplx, opt_nb, opt_nb_req)

         else if (opt_io%family%type == DECOMP_2D_IO_ADIOS2) then

            ! ADIOS2
            call adios2_read_var(opt_io, varname, freal, dreal, fcplx, dcplx)

         end if

      else

         call decomp_2d_io_object_open_and_start(io, &
                                                 opt_dirname, &
                                                 varname, &
                                                 decomp_2d_read_mode, &
                                                 opt_family=opt_family)

         if (io%family%type == DECOMP_2D_IO_MPI) then

            ! MPI-IO
            call mpi_read_var(ipencil, io, decomp, freal, dreal, fcplx, dcplx)

         else if (io%family%type == DECOMP_2D_IO_ADIOS2) then

            ! ADIOS2
            call adios2_read_var(io, varname, freal, dreal, fcplx, dcplx)

         end if

         call io%end_close()

      end if

   end subroutine read_one

   !
   ! Low-level. Using MPI-IO / ADIOS2 to write planes to a file
   !
   ! Inputs
   !   - ipencil : pencil orientation of the variable
   !   - varname : name of the variable
   !   - decomp : decomp_info for the variable
   !   - nplanes : number of planes in the variable
   !   - opt_mode : writing mode
   !   - opt_family : IO family can be provided to avoid the default one
   !   - opt_io : IO reader / writer can be provided
   !   - opt_dirname : This is mandatory if no IO reader / writer was provided
   !   - freal / dreal / fcplx / dcplx : array
   !   - opt_nb_req : id of the request for non-blocking MPI IO
   !
   ! If opt_io is provided and open
   !    - MPI IO will write to the file handle opt_io%fh
   !    - ADIOS2 IO will write to the engine opt_io%engine
   !
   ! Otherwise
   !    - opt_dirname is mandatory
   !    - MPI IO will write to the file opt_dirname/varname
   !    - ADIOS2 IO will write to the folder opt_dirname
   !
   subroutine write_plane(ipencil, varname, decomp, nplanes, opt_mode, &
                          opt_family, opt_io, opt_dirname, &
                          freal, dreal, fcplx, dcplx, opt_nb_req)

      implicit none

      integer, intent(IN) :: ipencil
      character(len=*), intent(in) :: varname
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(in) :: nplanes
      integer, intent(in), optional :: opt_mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      type(d2d_io), intent(inout), optional :: opt_io
      character(len=*), intent(in), optional :: opt_dirname
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx
      integer, intent(out), optional :: opt_nb_req

      logical :: use_opt_io, opt_nb
      integer :: mode
      type(d2d_io) :: io

      ! Use opt_io only if present and open
      if (present(opt_io)) then
         if (opt_io%is_open) then
            use_opt_io = .true.
         else
            use_opt_io = .false.
         end if
      else
         use_opt_io = .false.
      end if

      ! Safety check
      if (nplanes < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, nplanes, "Error invalid value of nplanes ")
      end if
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if
      if ((.not. use_opt_io) .and. (.not. present(opt_dirname))) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, "Invalid arguments")
      end if

      ! Default write mode is deferred
      if (present(opt_mode)) then
         mode = opt_mode
      else
         mode = decomp_2d_io_deferred
      end if

      ! Default MPI IO for writing is blocking
      if (use_opt_io .and. mode == decomp_2d_io_deferred .and. present(opt_nb_req)) then
         opt_nb = .true.
      else
         opt_nb = .false.
      end if

      ! Default value
      if (present(opt_nb_req)) opt_nb_req = MPI_REQUEST_NULL

      if (use_opt_io) then

         if (opt_io%family%type == DECOMP_2D_IO_MPI) then

            ! MPI-IO
            call mpi_write_plane(ipencil, opt_io, decomp, nplanes, freal, dreal, fcplx, dcplx, opt_nb, opt_nb_req)

         else if (opt_io%family%type == DECOMP_2D_IO_ADIOS2) then

            ! ADIOS2
            call adios2_write_var(opt_io, varname, mode, freal, dreal, fcplx, dcplx)

         end if

      else

         call decomp_2d_io_object_open_and_start(io, &
                                                 opt_dirname, &
                                                 varname, &
                                                 decomp_2d_write_mode, &
                                                 opt_family=opt_family)

         if (io%family%type == DECOMP_2D_IO_MPI) then

            ! MPI-IO
            call mpi_write_plane(ipencil, io, decomp, nplanes, freal, dreal, fcplx, dcplx)

         else if (io%family%type == DECOMP_2D_IO_ADIOS2) then

            ! ADIOS2
            call adios2_write_var(io, varname, mode, freal, dreal, fcplx, dcplx)

         end if

         call io%end_close()

      end if

   end subroutine write_plane

   !
   !
   ! Low-level MPI. This could be moved to a dedicated module.
   !
   !
   ! Write a 3D array to a file
   !
   subroutine mpi_write_var(ipencil, io, decomp, freal, dreal, fcplx, dcplx, opt_nb, opt_nb_req)

      implicit none

      integer, intent(IN) :: ipencil
      type(d2d_io), intent(inout) :: io
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx
      logical, intent(in), optional :: opt_nb
      integer, intent(out), optional :: opt_nb_req

      integer, dimension(3) :: sizes, subsizes, starts

      ! Process the size
      call io_get_size(ipencil, decomp, sizes, starts, subsizes)

      ! Do the MPI IO
      call mpi_read_or_write(.false., io, sizes, subsizes, starts, &
                             freal=freal, &
                             dreal=dreal, &
                             fcplx=fcplx, &
                             dcplx=dcplx, &
                             opt_nb=opt_nb, opt_nb_req=opt_nb_req)

   end subroutine mpi_write_var
   !
   ! Read a 3D array from a file
   !
   subroutine mpi_read_var(ipencil, io, decomp, freal, dreal, fcplx, dcplx, opt_nb, opt_nb_req)

      implicit none

      integer, intent(IN) :: ipencil
      type(d2d_io), intent(inout) :: io
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      real(real32), contiguous, dimension(:, :, :), intent(out), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(out), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(out), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(out), optional :: dcplx
      logical, intent(in), optional :: opt_nb
      integer, intent(out), optional :: opt_nb_req

      integer, dimension(3) :: sizes, subsizes, starts

      ! Process the size
      call io_get_size(ipencil, decomp, sizes, starts, subsizes)

      ! Do the MPI IO
      call mpi_read_or_write(.true., io, sizes, subsizes, starts, &
                             freal=freal, &
                             dreal=dreal, &
                             fcplx=fcplx, &
                             dcplx=dcplx, &
                             opt_nb=opt_nb, opt_nb_req=opt_nb_req)

   end subroutine mpi_read_var
   !
   ! Write planes to a file
   !
   subroutine mpi_write_plane(ipencil, io, decomp, nplanes, freal, dreal, fcplx, dcplx, opt_nb, opt_nb_req)

      implicit none

      integer, intent(IN) :: ipencil
      type(d2d_io), intent(inout) :: io
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(in) :: nplanes
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx
      logical, intent(in), optional :: opt_nb
      integer, intent(out), optional :: opt_nb_req

      integer, dimension(3) :: sizes, subsizes, starts

      ! Process the size
      call io_get_size(ipencil, decomp, sizes, starts, subsizes, nplanes)

      ! Do the MPI IO
      call mpi_read_or_write(.false., io, sizes, subsizes, starts, &
                             freal=freal, &
                             dreal=dreal, &
                             fcplx=fcplx, &
                             dcplx=dcplx, &
                             opt_nb=opt_nb, opt_nb_req=opt_nb_req)

   end subroutine mpi_write_plane
   !
   subroutine mpi_read_or_write(flag_read, io, sizes, subsizes, starts, &
                                freal, dreal, fcplx, dcplx, ints, logs, &
                                opt_nb, opt_nb_req)

      implicit none

      logical, intent(in) :: flag_read
      type(d2d_io), intent(inout) :: io
      integer, dimension(3), intent(in) :: sizes, subsizes, starts
      real(real32), contiguous, dimension(:, :, :), optional :: freal
      real(real64), contiguous, dimension(:, :, :), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), optional :: dcplx
      integer, contiguous, dimension(:, :, :), optional :: ints
      logical, contiguous, dimension(:, :, :), optional :: logs
      logical, intent(in), optional :: opt_nb
      integer, intent(out), optional :: opt_nb_req

      logical :: non_blocking
      integer :: my_mpi_info, ierror, data_type, newtype, type_bytes

      ! Safety check
      if (.not. io%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, &
                              "IO reader / writer was not opened "//io%label)
      end if

      ! Allow non-blocking MPI IO
      if (present(opt_nb)) then
         non_blocking = opt_nb
      else
         non_blocking = .false.
      end if
      if (non_blocking) then
         if (.not. present(opt_nb_req)) then
            call decomp_2d_abort(__FILE__, __LINE__, 0, "Missing argument")
         end if
      end if

      ! Hints for MPI IO
      if (associated(io%mpi_file_set_view_info)) then
         my_mpi_info = io%mpi_file_set_view_info
      else
         my_mpi_info = MPI_INFO_NULL
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
         call decomp_2d_abort(__FILE__, __LINE__, 0, "Invalid inputs for "//io%label)
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
      call MPI_FILE_SET_VIEW(io%fh, io%disp, data_type, newtype, 'native', my_mpi_info, ierror)
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

   end subroutine mpi_read_or_write

   !
   !
   ! Low-level ADIOS2. This could be moved to a dedicated module.
   !
   !
   ! Write a 3D array
   !
   subroutine adios2_write_var(io, varname, mode, &
                               freal, dreal, fcplx, dcplx)

      implicit none

      type(d2d_io), intent(inout) :: io
      character(len=*), intent(in) :: varname
      integer, intent(in) :: mode
      real(real32), contiguous, dimension(:, :, :), intent(IN), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(IN), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(IN), optional :: dcplx

#ifdef ADIOS2
      integer :: write_mode, ierror
      type(adios2_variable) :: var_handle
#endif

      ! Safety checks
      if (.not. io%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, &
                              "IO reader / writer was not opened "//io%label)
      end if
#ifdef ADIOS2
      if (.not. io%is_active) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, &
                              "IO reader / writer has not started "//io%label)
      end if
#endif
      if (mode < decomp_2d_io_deferred .or. &
          mode > decomp_2d_io_sync) then
         call decomp_2d_abort(__FILE__, __LINE__, mode, "Invalid value")
      end if

#ifdef ADIOS2
      ! Get the variable handle
      call adios2_inquire_variable(var_handle, io%family%io, varname, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_inquire_variable "//trim(varname))
      if (.not. var_handle%valid) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "ERROR: trying to write variable before registering! "//trim(varname))
      end if

      if (mode == decomp_2d_io_deferred) then
         write_mode = adios2_mode_deferred
      else if (mode == decomp_2d_io_sync) then
         write_mode = adios2_mode_sync
      end if

      if (io%engine%valid) then
         if (present(freal)) then
            call adios2_put(io%engine, var_handle, freal, write_mode, ierror)
         else if (present(dreal)) then
            call adios2_put(io%engine, var_handle, dreal, write_mode, ierror)
         else if (present(fcplx)) then
            call adios2_put(io%engine, var_handle, fcplx, write_mode, ierror)
         else if (present(dcplx)) then
            call adios2_put(io%engine, var_handle, dcplx, write_mode, ierror)
         end if
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_put")
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "ERROR: decomp2d thinks engine is live, but adios2 engine object is not valid")
      end if
#else
      associate (o => io, m => mode, v => varname, fr => freal, dr => dreal, fc => fcplx, dc => dcplx)
      end associate
#endif

   end subroutine adios2_write_var
   !
   ! Read a 3D array
   !
   subroutine adios2_read_var(io, varname, &
                              freal, dreal, fcplx, dcplx)

      implicit none

      type(d2d_io), intent(inout) :: io
      character(len=*), intent(in) :: varname
      real(real32), contiguous, dimension(:, :, :), intent(out), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(out), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(out), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(out), optional :: dcplx

#ifdef ADIOS2
      integer :: ierror
      integer(kind=8) :: nsteps, curstep
      type(adios2_variable) :: var_handle
#endif

      ! Safety checks
      if (.not. io%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, &
                              "IO reader / writer was not opened "//io%label)
      end if
#ifdef ADIOS2
      if (.not. io%is_active) then
         call decomp_2d_abort(__FILE__, __LINE__, 0, &
                              "IO reader / writer has not started "//io%label)
      end if
#endif

#ifdef ADIOS2
      ! Get the variable handle
      call adios2_inquire_variable(var_handle, io%family%io, varname, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_inquire_variable "//trim(varname))
      if (.not. var_handle%valid) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "ERROR: trying to read variable without registering first! "//trim(varname))
      end if

      call adios2_variable_steps(nsteps, var_handle, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_variable_steps")

      if (present(freal)) then
         call adios2_get(io%engine, var_handle, freal, adios2_mode_deferred, ierror)
      else if (present(dreal)) then
         call adios2_get(io%engine, var_handle, dreal, adios2_mode_deferred, ierror)
      else if (present(fcplx)) then
         call adios2_get(io%engine, var_handle, fcplx, adios2_mode_deferred, ierror)
      else if (present(dcplx)) then
         call adios2_get(io%engine, var_handle, dcplx, adios2_mode_deferred, ierror)
      end if
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_get")

      call adios2_current_step(curstep, io%engine, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_current_step")

      call adios2_inquire_variable(var_handle, io%family%io, varname, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_inquire_variable "//trim(varname))
      if (.not. var_handle%valid) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "ERROR: trying to write variable before registering! "//trim(varname))
      end if
#else
      associate (o => io, v => varname, fr => freal, dr => dreal, fc => fcplx, dc => dcplx)
      end associate
#endif

   end subroutine adios2_read_var

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write a 3D array as part of a big MPI-IO file, starting from
   !  displacement 'disp'; 'disp' will be updated after the writing
   !  operation to prepare the writing of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_var_real(fh, disp, ipencil, var, opt_decomp)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: ipencil
      real(mytype), contiguous, dimension(:, :, :), intent(IN) :: var
      TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

      TYPE(DECOMP_INFO) :: decomp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: ierror, newtype, data_type

      data_type = real_type

#include "io_write_var.inc"

      return
   end subroutine write_var_real

   subroutine write_var_complex(fh, disp, ipencil, var, opt_decomp)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: ipencil
      complex(mytype), contiguous, dimension(:, :, :), intent(IN) :: var
      TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

      TYPE(DECOMP_INFO) :: decomp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: ierror, newtype, data_type

      data_type = complex_type

#include "io_write_var.inc"

      return
   end subroutine write_var_complex

   subroutine write_outflow(dirname, varname, ntimesteps, var, io_name, opt_decomp)

      implicit none

      character(len=*), intent(in) :: dirname, varname, io_name
      integer, intent(IN) :: ntimesteps
      real(mytype), contiguous, dimension(:, :, :), intent(IN) :: var
      TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

      TYPE(DECOMP_INFO) :: decomp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: ierror, data_type
      integer :: idx
#ifdef ADIOS2
      type(adios2_io) :: io_handle
      type(adios2_variable) :: var_handle
#else
      integer :: newtype
#endif

      if (decomp_profiler_io) call decomp_profiler_start("io_write_outflow")

      data_type = real_type

#include "io_write_outflow.f90"

      if (decomp_profiler_io) call decomp_profiler_end("io_write_outflow")

   end subroutine write_outflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Read a 3D array as part of a big MPI-IO file, starting from
   !  displacement 'disp'; 'disp' will be updated after the reading
   !  operation to prepare the reading of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_var_real(fh, disp, ipencil, var, opt_decomp)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: ipencil
      real(mytype), contiguous, dimension(:, :, :), intent(INOUT) :: var
      TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

      TYPE(DECOMP_INFO) :: decomp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: ierror, newtype, data_type

      data_type = real_type

#include "io_read_var.inc"

      return
   end subroutine read_var_real

   subroutine read_var_complex(fh, disp, ipencil, var, opt_decomp)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: ipencil
      complex(mytype), contiguous, dimension(:, :, :), intent(INOUT) :: var
      TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

      TYPE(DECOMP_INFO) :: decomp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: ierror, newtype, data_type

      data_type = complex_type

#include "io_read_var.inc"

      return
   end subroutine read_var_complex

   subroutine read_inflow(dirname, varname, ntimesteps, var, io_name, opt_decomp)

      implicit none

      character(len=*), intent(in) :: dirname, varname, io_name
      integer, intent(IN) :: ntimesteps
      real(mytype), contiguous, dimension(:, :, :), intent(INOUT) :: var
      TYPE(DECOMP_INFO), intent(IN), optional :: opt_decomp

      TYPE(DECOMP_INFO) :: decomp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: ierror, data_type
      integer :: idx
#ifdef ADIOS2
      type(adios2_io) :: io_handle
      type(adios2_variable) :: var_handle
#else
      integer :: newtype
#endif

      if (decomp_profiler_io) call decomp_profiler_start("io_read_inflow")

      data_type = real_type

#include "io_read_inflow.f90"

      if (decomp_profiler_io) call decomp_profiler_end("io_read_inflow")

   end subroutine read_inflow

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write scalar variables as part of a big MPI-IO file, starting from
   !  displacement 'disp'; 'disp' will be updated after the reading
   !  operation to prepare the reading of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_scalar_real(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh             ! file handle
      integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
      integer, intent(IN) :: n              ! number of scalars
      real(mytype), dimension(n), &
         intent(IN) :: var                ! array of scalars

      integer :: m, ierror

      call MPI_FILE_SET_VIEW(fh, disp, real_type, &
                             real_type, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      if (nrank == 0) then
         m = n ! only one rank needs to write
      else
         m = 0
      end if
      call MPI_FILE_WRITE_ALL(fh, var, m, real_type, &
                              MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(mytype_bytes, kind=MPI_OFFSET_KIND)

      return
   end subroutine write_scalar_real

   subroutine write_scalar_complex(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: n
      complex(mytype), dimension(n), intent(IN) :: var

      integer :: m, ierror

      call MPI_FILE_SET_VIEW(fh, disp, complex_type, &
                             complex_type, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      if (nrank == 0) then
         m = n
      else
         m = 0
      end if
      call MPI_FILE_WRITE_ALL(fh, var, m, complex_type, &
                              MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(mytype_bytes, kind=MPI_OFFSET_KIND) &
             * 2_MPI_OFFSET_KIND

      return
   end subroutine write_scalar_complex

   subroutine write_scalar_integer(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: n
      integer, dimension(n), intent(IN) :: var

      integer :: m, ierror

      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, &
                             MPI_INTEGER, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      if (nrank == 0) then
         m = n
      else
         m = 0
      end if
      call MPI_FILE_WRITE_ALL(fh, var, m, MPI_INTEGER, &
                              MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
      call MPI_TYPE_SIZE(MPI_INTEGER, m, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(m, kind=MPI_OFFSET_KIND)

      return
   end subroutine write_scalar_integer

   subroutine write_scalar_logical(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: n
      logical, dimension(n), intent(IN) :: var

      integer :: m, ierror

      call MPI_FILE_SET_VIEW(fh, disp, MPI_LOGICAL, &
                             MPI_LOGICAL, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      if (nrank == 0) then
         m = n
      else
         m = 0
      end if
      call MPI_FILE_WRITE_ALL(fh, var, m, MPI_LOGICAL, &
                              MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
      call MPI_TYPE_SIZE(MPI_LOGICAL, m, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(m, kind=MPI_OFFSET_KIND)

      return
   end subroutine write_scalar_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Read scalar variables as part of a big MPI-IO file, starting from
   !  displacement 'disp'; 'disp' will be updated after the reading
   !  operation to prepare the reading of next chunk of data.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine read_scalar_real(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh             ! file handle
      integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement
      integer, intent(IN) :: n              ! number of scalars
      real(mytype), dimension(n), &
         intent(INOUT) :: var             ! array of scalars

      integer :: ierror

      call MPI_FILE_SET_VIEW(fh, disp, real_type, &
                             real_type, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      call MPI_FILE_READ_ALL(fh, var, n, real_type, &
                             MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_READ_ALL")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(mytype_bytes, kind=MPI_OFFSET_KIND)

      return
   end subroutine read_scalar_real

   subroutine read_scalar_complex(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: n
      complex(mytype), dimension(n), intent(INOUT) :: var

      integer :: ierror

      call MPI_FILE_SET_VIEW(fh, disp, complex_type, &
                             complex_type, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      call MPI_FILE_READ_ALL(fh, var, n, complex_type, &
                             MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_READ_ALL")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(mytype_bytes, kind=MPI_OFFSET_KIND) &
             * 2_MPI_OFFSET_KIND

      return
   end subroutine read_scalar_complex

   subroutine read_scalar_integer(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: n
      integer, dimension(n), intent(INOUT) :: var

      integer :: m, ierror

      call MPI_FILE_SET_VIEW(fh, disp, MPI_INTEGER, &
                             MPI_INTEGER, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      call MPI_FILE_READ_ALL(fh, var, n, MPI_INTEGER, &
                             MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_READ_ALL")
      call MPI_TYPE_SIZE(MPI_INTEGER, m, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(m, kind=MPI_OFFSET_KIND)

      return
   end subroutine read_scalar_integer

   subroutine read_scalar_logical(fh, disp, n, var)

      implicit none

      integer, intent(IN) :: fh
      integer(KIND=MPI_OFFSET_KIND), intent(INOUT) :: disp
      integer, intent(IN) :: n
      logical, dimension(n), intent(INOUT) :: var

      integer :: m, ierror

      call MPI_FILE_SET_VIEW(fh, disp, MPI_LOGICAL, &
                             MPI_LOGICAL, 'native', MPI_INFO_NULL, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
      call MPI_FILE_READ_ALL(fh, var, n, MPI_LOGICAL, &
                             MPI_STATUS_IGNORE, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_READ_ALL")
      call MPI_TYPE_SIZE(MPI_LOGICAL, m, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_SIZE")
      disp = disp + int(n, kind=MPI_OFFSET_KIND) &
             * int(m, kind=MPI_OFFSET_KIND)

      return
   end subroutine read_scalar_logical

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write 3D array data for every specified mesh point
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_every_real(ipencil, var, iskip, jskip, kskip, &
                               filename, from1)

      implicit none

      integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
      real(mytype), contiguous, dimension(:, :, :), intent(IN) :: var
      integer, intent(IN) :: iskip, jskip, kskip
      character(len=*), intent(IN) :: filename
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer(kind=MPI_OFFSET_KIND) :: filesize, disp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: i, j, k, ierror, newtype, fh, key, color, newcomm, data_type
      integer, dimension(3) :: xsz, ysz, zsz, xst, yst, zst, xen, yen, zen, skip

      if (decomp_profiler_io) call decomp_profiler_start("io_write_every_real")

      data_type = real_type

#include "io_write_every.inc"

      if (decomp_profiler_io) call decomp_profiler_end("io_write_every_real")

   end subroutine write_every_real

   subroutine write_every_complex(ipencil, var, iskip, jskip, kskip, &
                                  filename, from1)

      implicit none

      integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
      complex(mytype), contiguous, dimension(:, :, :), intent(IN) :: var
      integer, intent(IN) :: iskip, jskip, kskip
      character(len=*), intent(IN) :: filename
      logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
      ! .false. - save n,2n,3n...

      complex(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer(kind=MPI_OFFSET_KIND) :: filesize, disp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: i, j, k, ierror, newtype, fh, key, color, newcomm, data_type
      integer, dimension(3) :: xsz, ysz, zsz, xst, yst, zst, xen, yen, zen, skip

      if (decomp_profiler_io) call decomp_profiler_start("io_write_every_cplx")

      data_type = complex_type

#include "io_write_every.inc"

      if (decomp_profiler_io) call decomp_profiler_end("io_write_every_cplx")

   end subroutine write_every_complex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Write a 3D data set covering a smaller sub-domain only
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine write_subdomain(ipencil, var, is, ie, js, je, ks, ke, filename)

      implicit none

      integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
      real(mytype), contiguous, dimension(:, :, :), intent(IN) :: var
      integer, intent(IN) :: is, ie, js, je, ks, ke
      character(len=*), intent(IN) :: filename

      real(mytype), allocatable, dimension(:, :, :) :: wk, wk2
      integer(kind=MPI_OFFSET_KIND) :: filesize, disp
      integer, dimension(3) :: sizes, subsizes, starts
      integer :: color, key, errorcode, newcomm, ierror
      integer :: newtype, fh, data_type, i, j, k
      integer :: i1, i2, j1, j2, k1, k2

      if (decomp_profiler_io) call decomp_profiler_start("io_write_subdomain")

      data_type = real_type

      ! validate the input paramters
      if (is < 1 .OR. ie > nx_global .OR. js < 1 .OR. je > ny_global .OR. &
          ks < 1 .OR. ke > nz_global) then
         errorcode = 10
         call decomp_2d_abort(errorcode, &
                              'Invalid subdomain specified in I/O')
      end if

      ! create a communicator for all those MPI ranks containing the subdomain
      color = 1
      key = 1
      if (ipencil == 1) then
         if (xstart(1) > ie .OR. xend(1) < is .OR. xstart(2) > je .OR. xend(2) < js &
             .OR. xstart(3) > ke .OR. xend(3) < ks) then
            color = 2
         end if
      else if (ipencil == 2) then
         if (ystart(1) > ie .OR. yend(1) < is .OR. ystart(2) > je .OR. yend(2) < js &
             .OR. ystart(3) > ke .OR. yend(3) < ks) then
            color = 2
         end if
      else if (ipencil == 3) then
         if (zstart(1) > ie .OR. zend(1) < is .OR. zstart(2) > je .OR. zend(2) < js &
             .OR. zstart(3) > ke .OR. zend(3) < ks) then
            color = 2
         end if
      end if
      call MPI_COMM_SPLIT(decomp_2d_comm, color, key, newcomm, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_COMM_SPLIT")

      if (color == 1) then ! only ranks in this group do IO collectively

         ! generate MPI-IO subarray information

         ! global size of the sub-domain to write
         sizes(1) = ie - is + 1
         sizes(2) = je - js + 1
         sizes(3) = ke - ks + 1

         ! 'subsizes' & 'starts' as required by MPI_TYPE_CREATE_SUBARRAY
         ! note the special code whe subdomain only occupy part of the pencil
         if (ipencil == 1) then

            subsizes(1) = xsize(1)
            starts(1) = xstart(1) - is
            if (xend(1) > ie .AND. xstart(1) < is) then
               subsizes(1) = ie - is + 1
               starts(1) = 0
            else if (xstart(1) < is) then
               subsizes(1) = xend(1) - is + 1
               starts(1) = 0
            else if (xend(1) > ie) then
               subsizes(1) = ie - xstart(1) + 1
            end if
            subsizes(2) = xsize(2)
            starts(2) = xstart(2) - js
            if (xend(2) > je .AND. xstart(2) < js) then
               subsizes(2) = je - js + 1
               starts(2) = 0
            else if (xstart(2) < js) then
               subsizes(2) = xend(2) - js + 1
               starts(2) = 0
            else if (xend(2) > je) then
               subsizes(2) = je - xstart(2) + 1
            end if
            subsizes(3) = xsize(3)
            starts(3) = xstart(3) - ks
            if (xend(3) > ke .AND. xstart(3) < ks) then
               subsizes(3) = ke - ks + 1
               starts(3) = 0
            else if (xstart(3) < ks) then
               subsizes(3) = xend(3) - ks + 1
               starts(3) = 0
            else if (xend(3) > ke) then
               subsizes(3) = ke - xstart(3) + 1
            end if

         else if (ipencil == 2) then

            ! TODO

         else if (ipencil == 3) then

            ! TODO

         end if

         ! copy data from orginal to a temp array
         ! pay attention to blocks only partially cover the sub-domain
         if (ipencil == 1) then

            if (xend(1) > ie .AND. xstart(1) < is) then
               i1 = is
               i2 = ie
            else if (xend(1) > ie) then
               i1 = xstart(1)
               i2 = ie
            else if (xstart(1) < is) then
               i1 = is
               i2 = xend(1)
            else
               i1 = xstart(1)
               i2 = xend(1)
            end if

            if (xend(2) > je .AND. xstart(2) < js) then
               j1 = js
               j2 = je
            else if (xend(2) > je) then
               j1 = xstart(2)
               j2 = je
            else if (xstart(2) < js) then
               j1 = js
               j2 = xend(2)
            else
               j1 = xstart(2)
               j2 = xend(2)
            end if

            if (xend(3) > ke .AND. xstart(3) < ks) then
               k1 = ks
               k2 = ke
            else if (xend(3) > ke) then
               k1 = xstart(3)
               k2 = ke
            else if (xstart(3) < ks) then
               k1 = ks
               k2 = xend(3)
            else
               k1 = xstart(3)
               k2 = xend(3)
            end if

            allocate (wk(i1:i2, j1:j2, k1:k2))
            call alloc_x(wk2, opt_global=.true.)
            wk2 = var
            do k = k1, k2
               do j = j1, j2
                  do i = i1, i2
                     wk(i, j, k) = wk2(i, j, k)
                  end do
               end do
            end do

         else if (ipencil == 2) then

            ! TODO

         else if (ipencil == 3) then

            ! TODO

         end if

         deallocate (wk2)

         ! MPI-IO
         call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts, &
                                       MPI_ORDER_FORTRAN, data_type, newtype, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_CREATE_SUBARRAY")
         call MPI_TYPE_COMMIT(newtype, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_COMMIT")
         call MPI_FILE_OPEN(newcomm, filename, &
                            MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, &
                            fh, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_OPEN")
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh, filesize, ierror)  ! guarantee overwriting
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_SIZE")
         disp = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_VIEW(fh, disp, data_type, &
                                newtype, 'native', MPI_INFO_NULL, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_VIEW")
         call MPI_FILE_WRITE_ALL(fh, wk, &
                                 subsizes(1) * subsizes(2) * subsizes(3), &
                                 data_type, MPI_STATUS_IGNORE, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_WRITE_ALL")
         call MPI_FILE_CLOSE(fh, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_CLOSE")
         call MPI_TYPE_FREE(newtype, ierror)
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_TYPE_FREE")

         deallocate (wk)

      end if

      call decomp_2d_mpi_comm_free(newcomm)

      if (decomp_profiler_io) call decomp_profiler_end("io_write_subdomain")

   end subroutine write_subdomain

   integer function get_io_idx(io_name, engine_name)

      implicit none

      character(len=*), intent(in) :: io_name
      character(len=*), intent(in) :: engine_name

      character(len=(len(io_name) + len(io_sep) + len(engine_name))) :: full_name
      integer :: idx
      logical :: found

      character(len=1024), dimension(:), pointer :: names_ptr

#ifndef ADIOS2
      names_ptr => fh_names
#else
      names_ptr => engine_names
#endif

      full_name = io_name//io_sep//engine_name

      found = .false.
      do idx = 1, MAX_IOH
         if (names_ptr(idx) == full_name) then
            found = .true.
            exit
         end if
      end do

      if (.not. found) then
         idx = -1
      end if

      get_io_idx = idx

   end function get_io_idx

   function gen_iodir_name(io_dir, io_name)

      character(len=*), intent(in) :: io_dir, io_name
      character(len=(len(io_dir) + 5)) :: gen_iodir_name
#ifdef ADIOS2
      integer :: ierror
      type(adios2_io) :: io
      character(len=5) :: ext
#endif

#ifndef ADIOS2
      associate (nm => io_name) ! Silence unused dummy argument
      end associate
      write (gen_iodir_name, "(A)") io_dir
#else
      call adios2_at_io(io, adios, io_name, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_at_io "//trim(io_name))
      if (io%engine_type == "BP4") then
         ext = ".bp4"
      else if (io%engine_type == "HDF5") then
         ext = ".hdf5"
      else if (io%engine_type == "SST") then
         ext = ""
      else
         print *, "ERROR: Unkown engine type! ", io%engine_type
         print *, "-  IO: ", io_name
         print *, "- DIR:", io_dir
         stop
      end if
      write (gen_iodir_name, "(A,A)") io_dir, trim(ext)
#endif

   end function gen_iodir_name

end module decomp_2d_io
