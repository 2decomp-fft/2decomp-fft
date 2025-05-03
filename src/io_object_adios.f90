!! SPDX-License-Identifier: BSD-3-Clause

!
! This is the ADIOS2 IO module for readers / writers
!    => Create a family of readers / writers
!    => Create a reader / writer
!    => Open / close a reader / writer
!    => Start / end a reader / writer
!

!
! One can create a family a readers / writers using
!
!    call family%init(label)
!
! Arguments :
!    - label : name of the IO family. This name must be present in the ADIOS2 XML file
!
! Remark : the external code can associate the pointer family%adios to a specific adios2_adios object before calling init
!
! One can then register variables and planes in the family using
!
!    call family%register_var(varname, ipencil, type, opt_reduce_prec, opt_decomp)
!    call family%register_plane(varname, ipencil, type, opt_reduce_prec, opt_decomp, opt_nplanes)
!
! Arguments :
!    - varname : name of the variable
!    - ipencil : pencil orientation of the variable
!    - type : mpi_real / mpi_double_precision / mpi_complex / mpi_double_complex
!    - opt_reduce_prec : flag to convert double precision variables to single precision. default_opt_reduce_prec is used if not provided
!    - opt_decomp : decomp_info of the object. decomp_main is used if not provided
!    - opt_nplanes : number of planes stacked together. One is used if not provided
!

!
! One can open and start a ADIOS2 reader / writer using
!
!    call io%open_start(mode, opt_family)
!
! Arguments :
!    - mode : decomp_2d_write_mode, decomp_2d_read_mode, decomp_2d_append_mode
!    - opt_family : Family of readers / writers. default_family is used if not provided.
!
! To perform IO with the ADIOS2 reader / writer, the external code should use the subroutines defined in the module decomp_2d_io_adios
!

module decomp_2d_io_object_adios

   use adios2
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_io_utilities, only: io_get_size
   use decomp_2d_profiler
   use MPI
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

   ! Default adios2_adios object for IO
   type(adios2_adios), save, target :: default_adios

   !
   ! derived type to store info for a family of readers / writers
   !
   type, public :: d2d_io_family
      character(:), allocatable :: label                      ! Label of the family writer
      type(adios2_adios), pointer :: adios => null()
      type(adios2_io) :: io
   contains
      procedure :: init => d2d_io_family_init                     ! Init the family of readers / writers
      procedure :: fin => d2d_io_family_fin                       ! Clear the family of readers / writers
      procedure :: register_var => d2d_io_family_register_var     ! Register a 3D variable
      procedure :: register_plane => d2d_io_family_register_plane ! Register a 2D variable or stacked 2D variables
      procedure :: log => d2d_io_family_log                       ! Print some information about the object
   end type d2d_io_family

   ! Default family for ADIOS2 IO
   type(d2d_io_family), save, target :: default_family

   ! derived type to store info for a reader / writer
   type, public :: d2d_io_adios
      type(d2d_io_family), pointer :: family        ! Associated family
      character(:), allocatable :: label            ! Label of the writer
      logical :: is_open = .false.                  ! True if the writer is open
      logical :: is_active = .false.                ! True if the writer is active
      type(adios2_engine) :: engine                 ! Only one engine / writer currently
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
             decomp_2d_adios_get_default_adios, &
             decomp_2d_adios_get_default_family

contains

   !
   ! Initialize a new family of readers / writers
   !
   subroutine d2d_io_family_init(family, label)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: label

      integer :: ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_init")

      ! Safety check
      if (allocated(family%label)) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Family was not cleared "//label)
      end if

      family%label = label

      ! Advanced API
      ! The external code can set its own object of type adios2_adios before calling init
      if (.not. associated(family%adios)) then
         if (.not. default_adios%valid) then
            call decomp_2d_warning(__FILE__, __LINE__, 0, &
                                   "Warning, default adios2_adios is not set")
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

      if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) call family%log()

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_init")

   end subroutine d2d_io_family_init

   !
   ! Clear the given family of readers / writers
   !
   subroutine d2d_io_family_fin(family)

      implicit none

      class(d2d_io_family), intent(inout) :: family

      integer :: ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_fin")

      ! Safety check
      if (.not. allocated(family%label)) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, "Family was already cleared.")
         return
      end if

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

      deallocate (family%label)

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_fin")

   end subroutine d2d_io_family_fin

   !
   ! High-level. Register a 3D variable for the given family of readers / writers.
   !
   !    1 <= ipencil <= 3
   !    type can be MPI_REAL
   !                MPI_DOUBLE_PRECISION
   !                MPI_COMPLEX
   !                MPI_DOUBLE_COMPLEX
   !
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

      logical :: reduce_prec

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_register_var")

      if (present(opt_reduce_prec)) then
         reduce_prec = opt_reduce_prec
      else
         reduce_prec = default_opt_reduce_prec
      end if

      if (present(opt_decomp)) then
         call register_var(family, varname, ipencil, type, reduce_prec, opt_decomp, 0)
      else
         call register_var(family, varname, ipencil, type, reduce_prec, decomp_main, 0)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_register_var")

   end subroutine d2d_io_family_register_var

   !
   ! High-level. Register planes for the given family of readers / writers.
   !
   !    1 <= ipencil <= 3
   !    type can be MPI_REAL
   !                MPI_DOUBLE_PRECISION
   !                MPI_COMPLEX
   !                MPI_DOUBLE_COMPLEX
   !
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

      logical :: reduce_prec
      integer :: nplanes

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_register_plane")

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
         call register_var(family, varname, ipencil, type, reduce_prec, opt_decomp, nplanes)
      else
         call register_var(family, varname, ipencil, type, reduce_prec, decomp_main, nplanes)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_register_plane")

   end subroutine d2d_io_family_register_plane

   !
   ! Register one variable for the given family of readers / writers
   !
   subroutine register_var(family, varname, ipencil, typex, reduce_prec, decomp, nplanes)

      implicit none

      type(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil
      integer, intent(in) :: typex
      logical, intent(in) :: reduce_prec
      type(decomp_info), intent(in) :: decomp
      integer, intent(in) :: nplanes

      type(adios2_variable) :: var_handle
      integer, parameter :: ndims = 3
      logical, parameter :: adios2_constant_dims = .true.
      integer :: data_type
      integer :: ierror, iounit
      integer :: type
      integer, dimension(3) :: sizes, subsizes, starts

      ! Safety check
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if
      if (nplanes < 0) then
         call decomp_2d_abort(__FILE__, __LINE__, nplanes, "Error invalid value of nplanes")
      end if

      ! Process the sizes
      if (nplanes == 0) then
         call io_get_size(ipencil, decomp, sizes, starts, subsizes)
      else
         call io_get_size(ipencil, decomp, sizes, starts, subsizes, nplanes)
      end if

      ! Type can be changed depending on reduced precision
      if (reduce_prec .and. typex == MPI_DOUBLE_PRECISION) then
         type = MPI_REAL
      else if (reduce_prec .and. typex == MPI_DOUBLE_COMPLEX) then
         type = MPI_COMPLEX
      else
         type = typex
      end if

      ! Register the variable
      if (family%io%valid) then
         call adios2_inquire_variable(var_handle, family%io, varname, ierror)
         ! Variable found : ierror=0. Variable not found : ierror=adios2_error_exception
         if (.not. (ierror == 0 .or. ierror == adios2_error_exception)) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "adios2_inquire_variable "//varname)
         end if
         if (.not. var_handle%valid) then

            ! New variable
            if (decomp_debug >= D2D_DEBUG_LEVEL_INFO .and. d2d_log_is_active()) then
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
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Invalid adios2_io object for registering a variable")
      end if

   end subroutine register_var

   !
   ! Print some information about the object
   !
   subroutine d2d_io_family_log(family, opt_iounit)

      implicit none

      class(d2d_io_family), intent(in) :: family
      integer, intent(in), optional :: opt_iounit

      integer :: iounit

      ! Exit if needed
      if (.not. d2d_log_is_active()) return

      ! Get the IO unit if needed
      if (present(opt_iounit)) then
         iounit = opt_iounit
      else
         iounit = d2d_log_get_unit()
         write (iounit, *) ""
      end if

      ! Print stuff
      if (allocated(family%label)) then
         write (iounit, *) "    d2d_io_family "//family%label
      else
         write (iounit, *) "    d2d_io_family, label undefined"
      end if
      if (associated(family%adios)) then
         write (iounit, *) '      adios2_adios valid ? ', family%adios%valid
         write (iounit, *) '      adios2_adios f2c ', family%adios%f2c
      else
         write (iounit, *) '      adios2_adios pointer set to NULL'
      end if
      write (iounit, *) '      io is valid ? ', family%io%valid
      write (iounit, *) '      io engine type ', trim(family%io%engine_type)
      write (iounit, *) '      io f2c ', family%io%f2c

      ! Close the IO unit if needed
      if (.not. present(opt_iounit)) call d2d_log_close_unit(iounit)

   end subroutine d2d_io_family_log

   !
   ! Set the default adios2_adios object using a XML file
   !
   subroutine set_default_adios(adios_xml)

      implicit none

      character(len=*), intent(in) :: adios_xml

      integer :: ierror

      ! Safety check
      if (default_adios%valid) then
         if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) then
            call decomp_2d_warning(__FILE__, __LINE__, 0, &
                                   "Warning, default adios2_adios was already set")
         end if
         return
      end if

      call adios2_init(default_adios, trim(adios_xml), decomp_2d_comm, ierror)
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                              "Error initialising ADIOS2 - is "//trim(adios_xml)//" present and valid?")
      end if

   end subroutine set_default_adios

   !
   ! Open the given reader / writer
   !
   subroutine d2d_io_adios_open(writer, mode, family)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer
      integer, intent(in) :: mode
      type(d2d_io_family), target, intent(inout) :: family

      integer :: access_mode, ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_adios_open")

      ! Safety check
      if (writer%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Writer was not closed "//writer%label)
      end if

      ! Prepare the writer
      writer%family => family
      writer%is_open = .true.
      writer%label = family%label

      ! Prepare to open ADIOS2 writer, set access mode
      if (mode == decomp_2d_write_mode) then
         access_mode = adios2_mode_write
      else if (mode == decomp_2d_read_mode) then
         access_mode = adios2_mode_read
      else if (mode == decomp_2d_append_mode) then
         access_mode = adios2_mode_append
      else
         call decomp_2d_abort(__FILE__, __LINE__, mode, "Invalid value for mode")
      end if

      ! Open IO
      if (family%io%valid) then
         if (writer%family%io%engine_type == "BP4") then
            call adios2_open(writer%engine, family%io, &
                             trim(writer%family%label)//".bp4", &
                             access_mode, ierror)
         else if (writer%family%io%engine_type == "BP5") then
            call adios2_open(writer%engine, family%io, &
                             trim(writer%family%label)//".bp5", &
                             access_mode, ierror)
         else if (writer%family%io%engine_type == "HDF5") then
            call adios2_open(writer%engine, family%io, &
                             trim(writer%family%label)//".hdf5", &
                             access_mode, ierror)
         else if (writer%family%io%engine_type == "SST") then
            call adios2_open(writer%engine, family%io, &
                             trim(writer%family%label), &
                             access_mode, ierror)
         else
            call decomp_2d_abort(__FILE__, __LINE__, -1, &
                                 "Engine type "//writer%family%io%engine_type// &
                                 " label "//trim(writer%family%label)//".")
         end if
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "ERROR opening engine "//family%label)
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Couldn't find IO handle "//family%label)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_adios_open")

   end subroutine d2d_io_adios_open

   !
   ! Start IO for the given reader / writer
   !
   subroutine d2d_io_adios_start(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      integer :: ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_adios_start")

      ! Safety check
      if (.not. writer%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Writer was not opened "//writer%label)
      end if

      ! Safety check
      if (writer%is_active) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Writer was not ended "//writer%label)
      end if

      ! Tag the writer
      writer%is_active = .true.

      ! Start the writer if needed
      if (writer%engine%valid) then
         call adios2_begin_step(writer%engine, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "adios2_begin_step "//writer%label)
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "begin step, invalid engine "//writer%label)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_adios_start")

   end subroutine d2d_io_adios_start

   !
   ! Open the given reader / writer and start IO
   !
   subroutine d2d_io_adios_open_start(writer, mode, opt_family)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer
      integer, intent(in) :: mode
      type(d2d_io_family), target, intent(inout), optional :: opt_family

      if (present(opt_family)) then
         call writer%open(mode, opt_family)
      else
         call writer%open(mode, default_family)
      end if
      call writer%start()

   end subroutine d2d_io_adios_open_start

   !
   ! End IO for the given reader / writer
   !
   subroutine d2d_io_adios_end(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      integer :: ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_adios_end")

      ! Safety check
      if (.not. writer%is_active) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                "Writer was not started "//writer%label)
         return
      end if

      ! Tag the writer
      writer%is_active = .false.

      if (writer%engine%valid) then
         call adios2_end_step(writer%engine, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                 "adios2_end_step "//writer%label)
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "trying to end step with invalid engine "//writer%label)
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_adios_end")

   end subroutine d2d_io_adios_end

   !
   ! Close the given reader / writer
   !
   subroutine d2d_io_adios_close(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      integer :: ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_adios_close")

      ! Safety check
      if (.not. writer%is_open) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                "Writer was not opened "//writer%label)
         return
      end if

      ! Close the writer (type was checked at open stage)
      call adios2_close(writer%engine, ierror)
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_close")
      end if

      ! Clear the writer
      nullify (writer%family)
      deallocate (writer%label)
      writer%is_open = .false.

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_adios_close")

   end subroutine d2d_io_adios_close

   !
   ! End IO and close the given reader / writer
   !
   subroutine d2d_io_adios_end_close(writer)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer

      call writer%end()
      call writer%close()

   end subroutine d2d_io_adios_end_close

   !
   ! Write the provided array
   !
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

      integer :: write_mode, ierror
      type(adios2_variable) :: var_handle

      ! Safety checks
      if (.not. writer%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "IO reader / writer was not opened "//writer%label)
      end if
      if (.not. writer%is_active) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "IO reader / writer has not started "//writer%label)
      end if
      if (mode < decomp_2d_io_deferred .or. &
          mode > decomp_2d_io_sync) then
         call decomp_2d_abort(__FILE__, __LINE__, mode, "IO mode : Invalid value")
      end if

      ! Get the variable handle
      call adios2_inquire_variable(var_handle, writer%family%io, varname, ierror)
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

      if (writer%engine%valid) then
         if (present(freal)) then
            call adios2_put(writer%engine, var_handle, freal, write_mode, ierror)
         else if (present(dreal)) then
            call adios2_put(writer%engine, var_handle, dreal, write_mode, ierror)
         else if (present(fcplx)) then
            call adios2_put(writer%engine, var_handle, fcplx, write_mode, ierror)
         else if (present(dcplx)) then
            call adios2_put(writer%engine, var_handle, dcplx, write_mode, ierror)
         end if
         if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_put")
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "ERROR: adios2 engine object is not valid")
      end if

   end subroutine d2d_io_adios_write

   !
   ! Read the provided array
   !
   subroutine d2d_io_adios_read(writer, varname, sel_start, sel_count, freal, dreal, fcplx, dcplx)

      implicit none

      class(d2d_io_adios), intent(inout) :: writer
      character(len=*), intent(in) :: varname
      integer, intent(in), dimension(3) :: sel_start, sel_count
      real(real32), contiguous, dimension(:, :, :), intent(out), optional :: freal
      real(real64), contiguous, dimension(:, :, :), intent(out), optional :: dreal
      complex(real32), contiguous, dimension(:, :, :), intent(out), optional :: fcplx
      complex(real64), contiguous, dimension(:, :, :), intent(out), optional :: dcplx

      integer, parameter :: read_mode = adios2_mode_deferred
      integer :: ierror
      type(adios2_variable) :: var_handle

      ! Safety checks
      if (.not. writer%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "IO reader / writer was not opened "//writer%label)
      end if
      if (.not. writer%is_active) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "IO reader / writer has not started "//writer%label)
      end if

      ! Get the variable handle
      call adios2_inquire_variable(var_handle, writer%family%io, varname, ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_inquire_variable "//trim(varname))
      if (.not. var_handle%valid) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "ERROR: trying to read variable without registering first! "//trim(varname))
      end if

      ! Select a subset of the file
      call adios2_set_selection(var_handle, 3, int(sel_start, kind=8), int(sel_count, kind=8), ierror)
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_set_selection")

      ! Read (deferred)
      if (present(freal)) then
         call adios2_get(writer%engine, var_handle, freal, read_mode, ierror)
      else if (present(dreal)) then
         call adios2_get(writer%engine, var_handle, dreal, read_mode, ierror)
      else if (present(fcplx)) then
         call adios2_get(writer%engine, var_handle, fcplx, read_mode, ierror)
      else if (present(dcplx)) then
         call adios2_get(writer%engine, var_handle, dcplx, read_mode, ierror)
      end if
      if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_get")

   end subroutine d2d_io_adios_read

   !
   ! Print some information about the object
   !
   subroutine d2d_io_adios_log(writer, opt_iounit)

      implicit none

      class(d2d_io_adios), intent(in) :: writer
      integer, intent(in), optional :: opt_iounit

      integer :: iounit

      ! Exit if needed
      if (.not. d2d_log_is_active()) return

      ! Get the IO unit if needed
      if (present(opt_iounit)) then
         iounit = opt_iounit
      else
         iounit = d2d_log_get_unit()
         write (iounit, *) ""
      end if

      ! Print stuff
      if (allocated(writer%label)) then
         write (iounit, *) "  d2d_io_adios "//writer%label
      else
         write (iounit, *) "  d2d_io_adios, label undefined"
      end if
      write (iounit, *) "  is_open ", writer%is_open
      write (iounit, *) "  is_active ", writer%is_active
      if (writer%engine%valid) then
         write (iounit, *) "  engine ", trim(writer%engine%name)//", "//trim(writer%engine%type)
      else
         write (iounit, *) "   engine is not valid"
      end if
      if (associated(writer%family)) then
         if (decomp_debug >= D2D_DEBUG_LEVEL_INFO) then
            call writer%family%log(iounit)
         else
            if (allocated(writer%family%label)) then
               write (iounit, *) "  family ", trim(writer%family%label)
            else
               write (iounit, *) "  family with undefined label"
            end if
         end if
      else
         write (iounit, *) "  writer%family set to NULL"
      end if

      ! Close the IO unit if needed
      if (.not. present(opt_iounit)) call d2d_log_close_unit(iounit)

   end subroutine d2d_io_adios_log

   !
   ! Return a pointer to the default adios2_adios object
   !
   function decomp_2d_adios_get_default_adios()

      implicit none

      type(adios2_adios), pointer :: decomp_2d_adios_get_default_adios

      decomp_2d_adios_get_default_adios => default_adios

   end function decomp_2d_adios_get_default_adios

   !
   ! Return a pointer to the default family used
   !
   function decomp_2d_adios_get_default_family()

      implicit none

      type(d2d_io_family), pointer :: decomp_2d_adios_get_default_family

      decomp_2d_adios_get_default_family => default_family

   end function decomp_2d_adios_get_default_family

   !
   ! This should be called at the beginning
   !
   subroutine decomp_2d_adios_object_init(adios_xml)

      implicit none

      character(len=*), intent(in), optional :: adios_xml

      ! Prepare the default adios2_adios object using a XML file
      if (present(adios_xml)) then
         call set_default_adios(trim(adios_xml))
      else
         call set_default_adios("adios2_config.xml")
      end if

      ! Prepare the default ADIOS2 IO family
      call default_family%init("default")

   end subroutine decomp_2d_adios_object_init

   !
   ! This should be called at the end
   !
   subroutine decomp_2d_adios_object_fin

      implicit none

      integer :: ierror

      ! Clear the default IO family
      call default_family%fin()

      ! Clear the default adios2_adio object
      call adios2_finalize(default_adios, ierror)
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_finalize")
      end if

   end subroutine decomp_2d_adios_object_fin

end module decomp_2d_io_object_adios
