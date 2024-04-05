!! SPDX-License-Identifier: BSD-3-Clause

!
! This is the IO module for readers / writers
!    => Open / close a reader / writer
!    => Start / end a reader / writer
!

!
! TODO
!
! Development
!    Add a log subroutine for the derived type reader / writer ?
!    Have a default reader / writer ? (just like decomp_main for the derived type decomp_info)
!
! Simplification
!    Remove icoarse and the special cases icoarse > 0
!       Restore later
!          Using dedicated decomp_info objects ?
!          Combined with a small interpolation module allowing downsampling ? See interp.f90
!    Remove the reduced precision
!       Restore it later
!          Check adios2_mode_deferred: won't execute until adios2_end_step, adios2_perform_puts / adios2_perform_gets or adios2_close
!                adios2_mode_sync: special case, put data immediately, can be reused after this call
!          Avoid memory allocation inside 2decomp during IO ? => prevent ADIOS2 deferred IO
!          Avoid transpose operations inside 2decomp during IO ? => prevent ADIOS2 deferred IO
!          Reduced precision should be specified at the d2d_io_family level when registering the variable
!    d2d_io_family_register_var
!       The interface can be improved
!       The argument type (kind(0._real32) for instance) could be replaced by the MPI datatype (MPI_REAL for instance)
!    MPI
!       Write a generic subroutine to update the displacement offset
!
! Integration
!    Remove the subroutines decomp_2d_*_io and decomp_2d_register_variable
!    Update and clarify the API (expose the IO API in a markdown file)
!       Rename decomp_2d_write_one => decomp_2d_io_write
!       Rename decomp_2d_read_one => decomp_2d_io_read
!       Remove decomp_2d_write_var
!       Remove decomp_2d_read_var
!       Rename decomp_2d_write_scalar => decomp_2d_io_mpi_write_few (MPI only)
!       Rename decomp_2d_read_scalar => decomp_2d_io_mpi_read_few (MPI only)
!       Rename decomp_2d_write_plane => decomp_2d_io_write_plane
!       Add decomp_2d_io_read_plane
!       Remove decomp_2d_write_every
!       Rename decomp_2d_write_outflow => decomp_2d_io_write_outflow
!       Rename decomp_2d_read_inflow => decomp_2d_io_read_inflow
!       Keep decomp_2d_io_init
!       Rename decomp_2d_io_finalise => decomp_2d_io_fin
!    Adapt the examples accordingly
!

module decomp_2d_io_object

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io_family
   use decomp_2d_mpi
   use decomp_2d_profiler
#ifdef ADIOS2
   use adios2
#endif

   implicit none

   ! Default family
   type(d2d_io_family), save, pointer :: default_family

   ! derived type to store info for a reader / writer
   type, public :: d2d_io
      type(d2d_io_family), pointer :: family        ! Associated family
      character(:), allocatable :: label            ! Label of the writer
      logical :: is_open = .false.                  ! True if the writer is open
#ifdef ADIOS2
      logical :: is_active = .false.                ! True if the writer is active (adios2 only)
      type(adios2_engine) :: engine                 ! adios2 only (only one engine / writer currently)
#endif
      integer :: fh                                 ! File handle (mpi only)
      integer(kind=MPI_OFFSET_KIND) :: disp         ! Displacement offset (mpi only)
      ! Hints for MPI IO
      !    The external code can allocate / associate the pointers
      !    If allocated, the external code should deallocate
      integer, public, pointer :: mpi_file_open_info => null()
      integer, public, pointer :: mpi_file_set_view_info => null()
   contains
      procedure :: open => d2d_io_open              ! Open the IO
      procedure :: start => d2d_io_start            ! Start the IO
      procedure :: open_start => d2d_io_open_start  ! Open and start the IO
      procedure :: end => d2d_io_end                ! End the IO
      procedure :: close => d2d_io_close            ! Close the IO
      procedure :: end_close => d2d_io_end_close    ! End and close the IO
   end type d2d_io

   private

   public :: decomp_2d_io_object_open_and_start, &
             decomp_2d_io_object_get_default_family, &
             decomp_2d_io_object_set_default_family, &
             decomp_2d_io_object_init, &
             decomp_2d_io_object_fin

contains

   !
   ! Open the given reader / writer
   !
   subroutine d2d_io_open(writer, io_dir, mode, family, opt_mpi_info)

      implicit none

      class(d2d_io), intent(inout) :: writer
      character(len=*), intent(in) :: io_dir
      integer, intent(in) :: mode
      type(d2d_io_family), target, intent(in) :: family
      integer, intent(in), optional :: opt_mpi_info

      integer :: my_mpi_info, access_mode, ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_open")

      ! Safety check
      if (writer%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Writer was not closed "//writer%label)
      end if

      ! The external code can provide hints for MPI IO
      if (present(opt_mpi_info)) then
         ! Hints provided directly with the optional argument
         my_mpi_info = opt_mpi_info
      else if (associated(writer%mpi_file_open_info)) then
         ! Hints provided using the IO object
         my_mpi_info = writer%mpi_file_open_info
      else if (associated(family%mpi_file_open_info)) then
         ! Hints provided using the IO family object
         my_mpi_info = family%mpi_file_open_info
      else
         ! Default value for hints
         my_mpi_info = MPI_INFO_NULL
      end if
      if (associated(family%mpi_file_set_view_info) &
          .and. (.not. associated(writer%mpi_file_set_view_info))) then
         writer%mpi_file_set_view_info => family%mpi_file_set_view_info
      end if

      ! Prepare the writer
      writer%family => family
      writer%is_open = .true.
      writer%label = io_dir

      ! Prepare to open
      if (family%type == decomp_2d_io_mpi) then
         ! MPI writer, set initial displacement and access mode
         writer%disp = 0_MPI_OFFSET_KIND
         if (mode == decomp_2d_write_mode .or. mode == decomp_2d_append_mode) then
            access_mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
         else if (mode == decomp_2d_read_mode) then
            access_mode = MPI_MODE_RDONLY
         else
            call decomp_2d_abort(__FILE__, __LINE__, mode, "Invalid value for mode")
         end if

      else if (family%type == decomp_2d_io_adios2) then
#ifdef ADIOS2
         ! ADIOS2 writer, set access mode
         if (mode == decomp_2d_write_mode) then
            access_mode = adios2_mode_write
         else if (mode == decomp_2d_read_mode) then
            access_mode = adios2_mode_read
         else if (mode == decomp_2d_append_mode) then
            access_mode = adios2_mode_append
         else
            call decomp_2d_abort(__FILE__, __LINE__, mode, "Invalid value for mode")
         end if
#endif
      else
         ! Safety check
         call decomp_2d_abort(__FILE__, __LINE__, family%type, &
                              "Invalid family "//family%label)
      end if

      ! Open IO
      if (family%type == decomp_2d_io_mpi) then
         call MPI_FILE_OPEN(decomp_2d_comm, io_dir, &
                            access_mode, my_mpi_info, &
                            writer%fh, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_OPEN "//io_dir)
         end if
         if (mode == decomp_2d_write_mode) then
            ! Guarantee overwriting
            call MPI_FILE_SET_SIZE(writer%fh, 0_MPI_OFFSET_KIND, ierror)
            if (ierror /= 0) then
               call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_SIZE")
            end if
         end if

      else if (family%type == decomp_2d_io_adios2) then
#ifdef ADIOS2
         if (family%io%valid) then
            if (writer%family%io%engine_type == "BP4") then
               call adios2_open(writer%engine, family%io, &
                                trim(io_dir)//".bp4", &
                                access_mode, ierror)
            else if (writer%family%io%engine_type == "HDF5") then
               call adios2_open(writer%engine, family%io, &
                                trim(io_dir)//".hdf5", &
                                access_mode, ierror)
            else if (writer%family%io%engine_type == "SST") then
               call adios2_open(writer%engine, family%io, &
                                trim(io_dir), &
                                access_mode, ierror)
            else
               call decomp_2d_abort(__FILE__, __LINE__, -1, "Engine "//writer%family%io%engine_type//" Dir "//trim(io_dir)//".")
            end if
            if (ierror /= 0) then
               call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                    "ERROR opening engine "//io_dir//" "//family%label)
            end if
         else
            call decomp_2d_abort(__FILE__, __LINE__, -1, &
                                 "Couldn't find IO handle "//io_dir//" "//family%label)
         end if
#endif
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_open")

   end subroutine d2d_io_open

   !
   ! Start IO for the given reader / writer
   !
   subroutine d2d_io_start(writer)

      implicit none

      class(d2d_io), intent(inout) :: writer

#ifdef ADIOS2
      integer :: ierror
#endif

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_start")

      ! Safety check
      if (.not. writer%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Writer was not opened "//writer%label)
      end if

#ifdef ADIOS2
      ! Safety check
      if (writer%is_active) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Writer was not ended "//writer%label)
      end if

      ! Tag the writer
      writer%is_active = .true.

      ! Start the writer if needed
      if (writer%family%type == decomp_2d_io_adios2) then
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
      end if
#endif

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_start")

   end subroutine d2d_io_start

   !
   ! Open the given reader / writer and start IO
   !
   subroutine d2d_io_open_start(writer, io_dir, mode, opt_family, opt_mpi_info)

      implicit none

      class(d2d_io), intent(inout) :: writer
      character(len=*), intent(in) :: io_dir
      integer, intent(in) :: mode
      type(d2d_io_family), target, intent(in), optional :: opt_family
      integer, intent(in), optional :: opt_mpi_info

      if (present(opt_family)) then
         call writer%open(io_dir, mode, opt_family, opt_mpi_info)
      else
         call writer%open(io_dir, mode, default_family, opt_mpi_info)
      end if
      call writer%start()

   end subroutine d2d_io_open_start

   !
   ! End IO for the given reader / writer
   !
   subroutine d2d_io_end(writer)

      implicit none

      class(d2d_io), intent(inout) :: writer

#ifdef ADIOS2
      integer :: ierror
#endif

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_end")

#ifdef ADIOS2
      ! Safety check
      if (.not. writer%is_active) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                "Writer was not started "//writer%label)
         return
      end if

      ! Tag the writer
      writer%is_active = .false.

      if (writer%family%type == decomp_2d_io_adios2) then
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
      end if
#else
      associate (wrt => writer)
      end associate
#endif

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_end")

   end subroutine d2d_io_end

   !
   ! Close the given reader / writer
   !
   subroutine d2d_io_close(writer)

      implicit none

      class(d2d_io), intent(inout) :: writer

      integer :: ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_close")

      ! Safety check
      if (.not. writer%is_open) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                "Writer was not opened "//writer%label)
         return
      end if

      ! Close the writer (type was checked at open stage)
      if (writer%family%type == decomp_2d_io_mpi) then
         call MPI_FILE_CLOSE(writer%fh, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_CLOSE")
         end if
      else if (writer%family%type == decomp_2d_io_adios2) then
#ifdef ADIOS2
         call adios2_close(writer%engine, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_close")
         end if
#endif
      end if

      ! Clear the writer
      nullify (writer%family)
      deallocate (writer%label)
      writer%is_open = .false.
      if (associated(writer%mpi_file_open_info)) nullify (writer%mpi_file_open_info)
      if (associated(writer%mpi_file_set_view_info)) nullify (writer%mpi_file_set_view_info)

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_close")

   end subroutine d2d_io_close

   !
   ! End IO and close the given reader / writer
   !
   subroutine d2d_io_end_close(writer)

      implicit none

      class(d2d_io), intent(inout) :: writer

      call writer%end()
      call writer%close()

   end subroutine d2d_io_end_close

   !
   ! Open and start in write_one, read_one or write_plane
   !
   recursive subroutine decomp_2d_io_object_open_and_start(writer, &
                                                           dirname, &
                                                           varname, &
                                                           mode, &
                                                           opt_family)

      implicit none

      type(d2d_io), intent(out) :: writer
      character(len=*), intent(in) :: dirname, varname
      integer, intent(in) :: mode
      type(d2d_io_family), target, intent(in), optional :: opt_family

      character(:), allocatable :: fullname

      if (present(opt_family)) then
         if (opt_family%type == DECOMP_2D_IO_MPI) then
            fullname = trim(dirname)//"/"//trim(varname)
         else if (opt_family%type == DECOMP_2D_IO_ADIOS2) then
            fullname = trim(dirname)
         else
            fullname = ""
            call decomp_2d_abort(__FILE__, __LINE__, opt_family%type, "Invalid value")
         end if
         call writer%open_start(fullname, mode, opt_family)
         deallocate (fullname)
      else
         call decomp_2d_io_object_open_and_start(writer, &
                                                 dirname, &
                                                 varname, &
                                                 mode, &
                                                 opt_family=default_family)
      end if

   end subroutine decomp_2d_io_object_open_and_start

   !
   ! Provide a pointer to the default family used
   !
   function decomp_2d_io_object_get_default_family()

      implicit none

      type(d2d_io_family), pointer :: decomp_2d_io_object_get_default_family

      decomp_2d_io_object_get_default_family => default_family

   end function decomp_2d_io_object_get_default_family

   !
   ! Set the default family to the provided target
   !
   subroutine decomp_2d_io_object_set_default_family(family)

      implicit none

      type(d2d_io_family), target, intent(in), optional :: family

      if (present(family)) then
         default_family => family
      else
         if (associated(default_family)) nullify (default_family)
      end if

   end subroutine decomp_2d_io_object_set_default_family

   !
   ! This should be called at the beginning. Not mandatory.
   !
   subroutine decomp_2d_io_object_init

      implicit none

      if (associated(default_family)) nullify (default_family)

   end subroutine decomp_2d_io_object_init

   !
   ! This should be called at the end. Not mandatory.
   !
   subroutine decomp_2d_io_object_fin

      implicit none

      if (associated(default_family)) nullify (default_family)

   end subroutine decomp_2d_io_object_fin

end module decomp_2d_io_object
