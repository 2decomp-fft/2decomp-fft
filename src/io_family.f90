!! SPDX-License-Identifier: BSD-3-Clause

!
! This is the IO module for families of readers / writers
!    => Initialize / finalize a family of readers / writers
!    => Register a variable in a family of readers / writers
!

module decomp_2d_io_family

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_io_utilities
   use decomp_2d_profiler
   use MPI
#ifdef ADIOS2
   use adios2
#endif
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

#ifdef ADIOS2
   type(adios2_adios), save, pointer :: default_adios
#endif

   !
   ! derived type to store info for a family of readers / writers
   !
   type, public :: d2d_io_family
      integer :: type = decomp_2d_io_none                     ! Type of the family of writer
      character(:), allocatable :: label                      ! Label of the family writer
      ! Hints for MPI IO
      !    The external code can allocate / associate the pointers
      !    If allocated, the external code should deallocate
      integer, public, pointer :: mpi_file_open_info => null()
      integer, public, pointer :: mpi_file_set_view_info => null()
#ifdef ADIOS2
      type(adios2_adios), pointer :: adios                    ! adios2 only
      type(adios2_io) :: io                                   ! adios2 only
#endif
   contains
      procedure :: init => d2d_io_family_init                     ! ADIOS2 writer if possible, MPI otherwise
      procedure :: mpi_init => d2d_io_family_mpi_init             ! Force MPI writer
      procedure :: adios2_init => d2d_io_family_adios2_init       ! Force ADIOS2 writer
      procedure :: fin => d2d_io_family_fin                       ! Clear the writer
      procedure :: register_var3d => d2d_io_family_register_var3d ! Register a 3D variable
      procedure :: register_var2d => d2d_io_family_register_var2d ! Register a 2D variable
      procedure :: log => d2d_io_family_log                       ! Print some information about the object
   end type d2d_io_family

   private

   public :: decomp_2d_io_family_init, &
             decomp_2d_io_family_fin
#ifdef ADIOS2
   public :: decomp_2d_io_family_get_default_adios, &
             decomp_2d_io_family_set_default_adios
#endif

contains

   !
   ! Initialize a new family of readers / writers (default type)
   !
   subroutine d2d_io_family_init(family, label, opt_backend)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: label
      integer, intent(in), optional :: opt_backend

      if (present(opt_backend)) then
         if (opt_backend == DECOMP_2D_IO_MPI) then
            call family%mpi_init(label)
         else if (opt_backend == DECOMP_2D_IO_ADIOS2) then
            call family%adios2_init(label)
         else
            call decomp_2d_abort(__FILE__, __LINE__, opt_backend, "Invalid value for opt_backend")
         end if
      else
#ifdef ADIOS2
         call family%adios2_init(label)
#else
         call family%mpi_init(label)
#endif
      end if

   end subroutine d2d_io_family_init

   !
   ! Initialize a new family of MPI readers / writers
   !
   subroutine d2d_io_family_mpi_init(family, label)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: label

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_mpi_init")

      ! Safety check
      if (family%type /= decomp_2d_io_none) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Family was not cleared "//label)
      end if

      family%type = decomp_2d_io_mpi
      family%label = label

      if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) call family%log()

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_mpi_init")

   end subroutine d2d_io_family_mpi_init

   !
   ! Initialize a new family of ADIOS2 readers / writers
   !
   subroutine d2d_io_family_adios2_init(family, label)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: label

#ifdef ADIOS2
      integer :: ierror
#else
      call decomp_2d_abort(__FILE__, __LINE__, -1, "ADIOS2 is not available")
#endif

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_adios2_init")

      ! Safety check
      if (family%type /= decomp_2d_io_none) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Family was not cleared "//label)
      end if

#ifdef ADIOS2
      family%type = decomp_2d_io_adios2
      family%label = label

      ! Advanced API
      ! The external code can set its own object of type adios2_adios before calling init
      if (.not. associated(family%adios)) then
         if (.not. associated(default_adios)) then
            call decomp_2d_warning(__FILE__, __LINE__, 0, &
                                   "Warning, default adios2_adios is null()")
         end if
         family%adios => default_adios
      end if

      if (family%adios%valid) then
         call adios2_declare_io(family%io, family%adios, label, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_declare_io "//label)
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, "adios object not valid")
      end if
#endif

      if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) call family%log()

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_adios2_init")

   end subroutine d2d_io_family_adios2_init

   !
   ! Clear the given family of readers / writers
   !
   subroutine d2d_io_family_fin(family)

      implicit none

      class(d2d_io_family), intent(inout) :: family

#ifdef ADIOS2
      integer :: ierror
#endif

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_fin")

      ! Safety check
      if (family%type == decomp_2d_io_none) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, "Family was already cleared.")
         return
      end if

#ifdef ADIOS2
      if (family%type == decomp_2d_io_adios2) then
         call adios2_flush_all_engines(family%io, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "adios2_flush_all_engines "//family%label)
         end if
         call adios2_remove_all_variables(family%io, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "adios2_remove_all_variables "//family%label)
         end if
         call adios2_remove_all_attributes(family%io, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "adios2_remove_all_attributes "//family%label)
         end if
         nullify (family%adios)
      end if
#endif

      family%type = decomp_2d_io_none
      deallocate (family%label)
      if (associated(family%mpi_file_open_info)) nullify (family%mpi_file_open_info)
      if (associated(family%mpi_file_set_view_info)) nullify (family%mpi_file_set_view_info)

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_fin")

   end subroutine d2d_io_family_fin

   !
   ! High-level. Register one 3D variable for the given family of readers / writers.
   !
   !    1 <= ipencil <= 3
   !    type can be MPI_REAL
   !                MPI_DOUBLE_PRECISION
   !                MPI_COMPLEX
   !                MPI_DOUBLE_COMPLEX
   !
   subroutine d2d_io_family_register_var3d(family, &
                                           varname, &
                                           ipencil, &
                                           type, &
                                           opt_reduce_prec, &
                                           opt_decomp)

      implicit none

      class(d2d_io_family), intent(in) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil ! (x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(in) :: type
      logical, intent(in), optional :: opt_reduce_prec
      type(decomp_info), intent(in), optional :: opt_decomp

      logical :: reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_register_var3d")

      if (present(opt_reduce_prec)) then
         reduce_prec = opt_reduce_prec
      else
         reduce_prec = default_opt_reduce_prec
      end if

      if (present(opt_decomp)) then
         call register_var3d(family, varname, ipencil, type, reduce_prec, opt_decomp)
      else
         call register_var3d(family, varname, ipencil, type, reduce_prec, decomp_main)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_register_var3d")

   end subroutine d2d_io_family_register_var3d
   !
   ! Low-level. Register one 3D variable for the given family of readers / writers
   !
   subroutine register_var3d(family, varname, ipencil, type, reduce_prec, decomp)

      implicit none

      type(d2d_io_family), intent(in) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil
      integer, intent(in) :: type
      logical, intent(in) :: reduce_prec
      type(decomp_info), intent(in) :: decomp

#ifdef ADIOS2
      integer, dimension(3) :: sizes, subsizes, starts
#endif

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

#ifdef ADIOS2
      if (family%type == decomp_2d_io_adios2) then

         ! Process the size
         call io_get_size(ipencil, decomp, sizes, starts, subsizes)
         if (reduce_prec .and. type == MPI_DOUBLE_PRECISION) then
            call adios2_register_var(family, varname, sizes, starts, subsizes, MPI_REAL)
         else if (reduce_prec .and. type == MPI_DOUBLE_COMPLEX) then
            call adios2_register_var(family, varname, sizes, starts, subsizes, MPI_COMPLEX)
         else
            call adios2_register_var(family, varname, sizes, starts, subsizes, type)
         end if

      end if
#else
      associate (fm => family, vr => varname, pc => ipencil, &
                 tp => type, rd => reduce_prec, od => decomp)
      end associate
#endif

   end subroutine register_var3d

   !
   ! High-level. Register one 2D variable for the given family of readers / writers.
   !
   !    1 <= ipencil <= 3
   !    type can be MPI_REAL
   !                MPI_DOUBLE_PRECISION
   !                MPI_COMPLEX
   !                MPI_DOUBLE_COMPLEX
   !
   subroutine d2d_io_family_register_var2d(family, &
                                           varname, &
                                           ipencil, &
                                           type, &
                                           opt_reduce_prec, &
                                           opt_decomp, &
                                           opt_nplanes)

      implicit none

      class(d2d_io_family), intent(in) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil ! (x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(in) :: type
      logical, intent(in), optional :: opt_reduce_prec
      type(decomp_info), intent(in), optional :: opt_decomp
      integer, intent(in), optional :: opt_nplanes

      logical :: reduce_prec
      integer :: nplanes

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_register_var2d")

      if (present(opt_reduce_prec)) then
         reduce_prec = opt_reduce_prec
      else
         reduce_prec = default_opt_reduce_prec
      end if

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (present(opt_decomp)) then
         call register_var2d(family, varname, ipencil, type, reduce_prec, opt_decomp, nplanes)
      else
         call register_var2d(family, varname, ipencil, type, reduce_prec, decomp_main, nplanes)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_register_var2d")

   end subroutine d2d_io_family_register_var2d
   !
   ! Low-level. Register one 2D variable for the given family of readers / writers
   !
   subroutine register_var2d(family, varname, ipencil, type, reduce_prec, decomp, nplanes)

      implicit none

      type(d2d_io_family), intent(in) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil
      integer, intent(in) :: type
      logical, intent(in) :: reduce_prec
      type(decomp_info), intent(in) :: decomp
      integer, intent(in) :: nplanes

#ifdef ADIOS2
      integer, dimension(3) :: sizes, subsizes, starts
#endif

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if
      if (nplanes < 1) then
         call decomp_2d_abort(__FILE__, __LINE__, nplanes, "Error invalid value of nplanes")
      end if

#ifdef ADIOS2
      if (family%type == decomp_2d_io_adios2) then

         call io_get_size(ipencil, decomp, sizes, starts, subsizes, nplanes)
         if (reduce_prec .and. type == MPI_DOUBLE_PRECISION) then
            call adios2_register_var(family, varname, sizes, starts, subsizes, MPI_REAL)
         else if (reduce_prec .and. type == MPI_DOUBLE_COMPLEX) then
            call adios2_register_var(family, varname, sizes, starts, subsizes, MPI_COMPLEX)
         else
            call adios2_register_var(family, varname, sizes, starts, subsizes, type)
         end if

      end if
#else
      associate (fm => family, vr => varname, pc => ipencil, &
                 tp => type, rd => reduce_prec, od => decomp)
      end associate
#endif

   end subroutine register_var2d

   !
   ! Print some information about the object
   !
   subroutine d2d_io_family_log(family)

      implicit none

      class(d2d_io_family), intent(in) :: family

      integer :: iounit

      ! Exit if needed
      if (.not. d2d_log_is_active()) return

      ! Get the IO unit
      iounit = d2d_log_get_unit()

      ! Print stuff
      write (iounit, *) ""
      if (allocated(family%label)) then
         write (iounit, *) "  d2d_io_family "//family%label
      else
         write (iounit, *) "  d2d_io_family, label undefined"
      end if
      select case (family%type)
      case (DECOMP_2D_IO_NONE)
         write (iounit, *) '    type : None'
      case (DECOMP_2D_IO_MPI)
         write (iounit, *) '    type : MPI'
      case (DECOMP_2D_IO_ADIOS2)
         write (iounit, *) '    type : ADIOS2'
      case default
         write (iounit, *) '    type : ? ', family%type
      end select

      ! Close the IO unit
      call d2d_log_close_unit(iounit)

   end subroutine d2d_io_family_log

   !
   ! This should be called at the beginning. Not mandatory.
   !
   subroutine decomp_2d_io_family_init()

      implicit none

#ifdef ADIOS2
      if (associated(default_adios)) then
         call decomp_2d_warning(__FILE__, __LINE__, 0, &
                                "Default adios2_adios should be null()")
         nullify (default_adios)
      end if
#endif

   end subroutine decomp_2d_io_family_init

   !
   ! This should be called at the end. Not mandatory.
   !
   subroutine decomp_2d_io_family_fin()

      implicit none

#ifdef ADIOS2
      if (associated(default_adios)) nullify (default_adios)
#endif

   end subroutine decomp_2d_io_family_fin

   !
   !
   ! Stuff below is ADIOS2 only
   !
   !
#ifdef ADIOS2
   !
   ! Register one variable
   !
   subroutine adios2_register_var(family, varname, sizes, starts, subsizes, type)

      implicit none

      type(d2d_io_family), intent(in) :: family
      character(len=*), intent(in) :: varname
      integer, dimension(3), intent(in) :: sizes, starts, subsizes
      integer, intent(in) :: type

      type(adios2_variable) :: var_handle
      integer, parameter :: ndims = 3
      logical, parameter :: adios2_constant_dims = .true.
      integer :: data_type
      integer :: ierror, iounit

      ! Register
      if (family%io%valid) then
         call adios2_inquire_variable(var_handle, family%io, varname, ierror)
         ! Variable found : ierror=0. Variable not found : ierror=adios2_error_exception
         if (.not. (ierror == 0 .or. ierror == adios2_error_exception)) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "adios2_inquire_variable "//varname)
         end if
         if (.not. var_handle%valid) then

            ! New variable
            if (d2d_log_is_active()) then
               iounit = d2d_log_get_unit()
               write (iounit, *) "Registering variable for IO: ", varname
               call d2d_log_close_unit(iounit)
            end if

            ! Need to set the ADIOS2 data type
            if (type == MPI_REAL) then
               data_type = adios2_type_real
            else if (type == MPI_DOUBLE_PRECISION) then
               data_type = adios2_type_dp
            else if (type == MPI_COMPLEX) then
               data_type = adios2_type_complex
            else if (type == MPI_DOUBLE_COMPLEX) then
               data_type = adios2_type_complex_dp
            else
               call decomp_2d_abort(__FILE__, __LINE__, type, &
                                    "Trying to write unknown data type!")
            end if

            call adios2_define_variable(var_handle, family%io, varname, data_type, ndims, &
                                        int(sizes, kind=8), int(starts, kind=8), &
                                        int(subsizes, kind=8), adios2_constant_dims, ierror)
            if (ierror /= 0) then
               call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                    "adios2_define_variable "//varname)
            end if
         else
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "Variable already registered "//varname)
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, 0, &
                              "Invalid adios2_io object for registering a variable")
      end if

   end subroutine adios2_register_var

   !
   ! Return a pointer to the default adios2_adios object
   !
   function decomp_2d_io_family_get_default_adios()

      implicit none

      type(adios2_adios), pointer :: decomp_2d_io_family_get_default_adios

      if (.not. associated(default_adios)) then
         call decomp_2d_warning(__FILE__, __LINE__, 0, &
                                "Default adios2_adios is null()")
      end if

      decomp_2d_io_family_get_default_adios => default_adios

   end function decomp_2d_io_family_get_default_adios

   !
   ! Set the default adios2_adios object to the provided target
   !
   subroutine decomp_2d_io_family_set_default_adios(adios)

      implicit none

      type(adios2_adios), target, intent(in), optional :: adios

      if (present(adios)) then
         default_adios => adios
      else
         if (associated(default_adios)) nullify (default_adios)
      end if

   end subroutine decomp_2d_io_family_set_default_adios
#endif

end module decomp_2d_io_family
