!! SPDX-License-Identifier: BSD-3-Clause

! This is the MPI-IO module for readers / writers
!    => Open / close a reader / writer

!
! The external code can open a file for MPI IO with the
! type-bound procedure of the objects d2d_io_mpi
!
!    call io_obj%open(filename, mode, opt_dirname, opt_info1, opt_info2)
!
! Arguments :
!   - filename : name of the file
!   - mode : decomp_2d_write_mode, decomp_2d_read_mode or decomp_2d_append_mode
!   - opt_dirname : optional, folder where the file is located
!   - opt_info1 : optional, MPI hints provided to MPI_FILE_OPEN
!   - opt_info2 : optional, MPI hints provided to MPI_FILE_SET_VIEW
!
! The file can be closed with
!
!    call io%close()
!
! Default values for MPI hints can be set using the subroutines
!
!    decomp_2d_io_object_set_default_mpi_file_open_info
!    decomp_2d_io_object_set_default_mpi_file_set_view_info
!

module decomp_2d_io_object_mpi

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_profiler

   implicit none

   integer, save, private :: default_mpi_file_open_info = MPI_INFO_NULL
   integer, save, private :: default_mpi_file_set_view_info = MPI_INFO_NULL

   ! derived type to store info for a MPI-IO reader / writer
   type, public :: d2d_io_mpi
      character(:), allocatable :: label                    ! Label of the writer
      logical :: is_open = .false.                          ! True if the writer is open
      integer, public :: fh                                 ! File handle
      integer(kind=MPI_OFFSET_KIND), public :: disp         ! Displacement offset
      ! Hints for MPI IO
      !    The external code can allocate / associate the pointers
      !    If allocated, the external code should deallocate
      integer, public :: mpi_file_open_info = MPI_INFO_NULL
      integer, public :: mpi_file_set_view_info = MPI_INFO_NULL
   contains
      procedure :: open => d2d_io_mpi_open              ! Open the IO
      procedure :: close => d2d_io_mpi_close            ! Close the IO
      procedure :: log => d2d_io_mpi_log                ! Print some information about the object
   end type d2d_io_mpi

   private

   public :: decomp_2d_io_object_get_default_mpi_file_open_info, &
             decomp_2d_io_object_set_default_mpi_file_open_info, &
             decomp_2d_io_object_get_default_mpi_file_set_view_info, &
             decomp_2d_io_object_set_default_mpi_file_set_view_info, &
             decomp_2d_io_object_mpi_init, &
             decomp_2d_io_object_mpi_fin

contains

   !
   ! Open the given MPI-IO reader / writer
   !
   subroutine d2d_io_mpi_open(writer, filename, mode, &
                              opt_dirname, &
                              opt_mpi_file_open_info, &
                              opt_mpi_file_set_view_info)

      implicit none

      class(d2d_io_mpi), intent(inout) :: writer
      character(len=*), intent(in) :: filename
      integer, intent(in) :: mode
      character(len=*), intent(in), optional :: opt_dirname
      integer, intent(in), optional :: opt_mpi_file_open_info, &
                                       opt_mpi_file_set_view_info

      integer :: access_mode, ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_mpi_open")

      ! Safety check
      if (writer%is_open) then
         if (allocated(writer%label)) then
            call decomp_2d_abort(__FILE__, __LINE__, -1, &
                                 "Writer was not closed "//writer%label)
         else
            call decomp_2d_abort(__FILE__, __LINE__, -1, &
                                 "Writer was not closed (undefined label)")
         end if
      end if

      ! The external code can provide hints for MPI IO
      if (present(opt_mpi_file_open_info)) then
         ! Hints provided with the optional argument
         writer%mpi_file_open_info = opt_mpi_file_open_info
      else
         ! Default hints
         writer%mpi_file_open_info = default_mpi_file_open_info
      end if
      if (present(opt_mpi_file_set_view_info)) then
         ! Hints provided with the optional argument
         writer%mpi_file_set_view_info = opt_mpi_file_set_view_info
      else
         ! Default hints
         writer%mpi_file_set_view_info = default_mpi_file_set_view_info
      end if

      ! Prepare the writer
      writer%is_open = .true.
      writer%label = filename

      ! Set initial displacement and access mode
      writer%disp = 0_MPI_OFFSET_KIND
      if (mode == decomp_2d_write_mode .or. mode == decomp_2d_append_mode) then
         access_mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
      else if (mode == decomp_2d_read_mode) then
         access_mode = MPI_MODE_RDONLY
      else
         call decomp_2d_abort(__FILE__, __LINE__, mode, "Invalid value for mode")
      end if

      ! Open IO
      if (present(opt_dirname)) then
         call MPI_FILE_OPEN(decomp_2d_comm, &
                            opt_dirname//"/"//filename, &
                            access_mode, &
                            writer%mpi_file_open_info, &
                            writer%fh, &
                            ierror)
      else
         call MPI_FILE_OPEN(decomp_2d_comm, &
                            filename, &
                            access_mode, &
                            writer%mpi_file_open_info, &
                            writer%fh, &
                            ierror)
      end if
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_OPEN "//filename)
      end if
      if (mode == decomp_2d_write_mode) then
         ! Guarantee overwriting
         call MPI_FILE_SET_SIZE(writer%fh, 0_MPI_OFFSET_KIND, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_SET_SIZE")
         end if
      end if

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_mpi_open")

   end subroutine d2d_io_mpi_open

   !
   ! Close the given reader / writer
   !
   subroutine d2d_io_mpi_close(writer)

      implicit none

      class(d2d_io_mpi), intent(inout) :: writer

      integer :: ierror

      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_mpi_close")

      ! Safety check
      if (.not. writer%is_open) then
         if (allocated(writer%label)) then
            call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                   "Writer was not opened "//writer%label)
         else
            call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                   "Writer was not opened (undefined label)")
         end if
         return
      end if

      ! Close the MPI-IO reader / writer
      call MPI_FILE_CLOSE(writer%fh, ierror)
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_CLOSE")
      end if

      ! Clear the writer
      deallocate (writer%label)
      writer%is_open = .false.
      writer%disp = 0_MPI_OFFSET_KIND
      writer%mpi_file_open_info = MPI_INFO_NULL
      writer%mpi_file_set_view_info = MPI_INFO_NULL

      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_mpi_close")

   end subroutine d2d_io_mpi_close

   !
   ! Print some information about the object
   !
   subroutine d2d_io_mpi_log(writer, opt_iounit)

      implicit none

      class(d2d_io_mpi), intent(in) :: writer
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

      ! Print tuff
      if (allocated(writer%label)) then
         write (iounit, *) "  d2d_mpi_io "//writer%label
      else
         write (iounit, *) "  d2d_mpi_io, label undefined"
      end if
      write (iounit, *) "  is_open ", writer%is_open
      write (iounit, *) "  fh ", writer%fh
      write (iounit, *) "  disp ", writer%disp
      write (iounit, *) "  file_open_info ", writer%mpi_file_open_info
      write (iounit, *) "  file_set_view_info ", writer%mpi_file_set_view_info

      ! Close the IO unit if needed
      if (.not. present(opt_iounit)) call d2d_log_close_unit(iounit)

   end subroutine d2d_io_mpi_log

   !
   ! Return the default MPI_INFO used when MPI is opening the file
   !
   function decomp_2d_io_object_get_default_mpi_file_open_info() &
      result(file_open_info)

      implicit none

      integer :: file_open_info

      file_open_info = default_mpi_file_open_info

   end function decomp_2d_io_object_get_default_mpi_file_open_info

   !
   ! Set the default MPI_INFO used when MPI is opening the file
   !
   subroutine decomp_2d_io_object_set_default_mpi_file_open_info(file_open_info)

      implicit none

      integer, intent(in), optional :: file_open_info

      if (present(file_open_info)) then
         default_mpi_file_open_info = file_open_info
      else
         default_mpi_file_open_info = MPI_INFO_NULL
      end if

   end subroutine decomp_2d_io_object_set_default_mpi_file_open_info

   !
   ! Return the default MPI_INFO used when calling MPI_FILE_SET_VIEW
   !
   function decomp_2d_io_object_get_default_mpi_file_set_view_info() &
      result(file_set_view_info)

      implicit none

      integer :: file_set_view_info

      file_set_view_info = default_mpi_file_set_view_info

   end function decomp_2d_io_object_get_default_mpi_file_set_view_info

   !
   ! Set the default MPI_INFO used when calling MPI_FILE_SET_VIEW
   !
   subroutine decomp_2d_io_object_set_default_mpi_file_set_view_info(file_set_view_info)

      implicit none

      integer, intent(in), optional :: file_set_view_info

      if (present(file_set_view_info)) then
         default_mpi_file_set_view_info = file_set_view_info
      else
         default_mpi_file_set_view_info = MPI_INFO_NULL
      end if

   end subroutine decomp_2d_io_object_set_default_mpi_file_set_view_info

   !
   ! This should be called at the beginning. Not mandatory.
   !
   !    The external code can specify default MPI_INFO values for
   !    calls to MPI_FILE_OPEN and MPI_FILE_SET_VIEW
   !
   subroutine decomp_2d_io_object_mpi_init(file_open_info, file_set_view_info)

      implicit none

      integer, intent(in), optional :: file_open_info, file_set_view_info

      call decomp_2d_io_object_set_default_mpi_file_open_info(file_open_info)
      call decomp_2d_io_object_set_default_mpi_file_set_view_info(file_set_view_info)

   end subroutine decomp_2d_io_object_mpi_init

   !
   ! This should be called at the end. Not mandatory.
   !
   subroutine decomp_2d_io_object_mpi_fin

      implicit none

      call decomp_2d_io_object_set_default_mpi_file_open_info()
      call decomp_2d_io_object_set_default_mpi_file_set_view_info()

   end subroutine decomp_2d_io_object_mpi_fin

end module decomp_2d_io_object_mpi
