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

!
! This is the IO submodule for writers
!    => Initialize / finalize a writer
!    => Open / close a writer
!    => Start / end a writer
!    => Register a variable in a writer
!

!
! TODO
!
! Development
!    Add a log subroutine for the derived type family of writers ?
!    Have a default family of writers ? (just like decomp_main for the derived type decomp_info)
!    Add a log subroutine for the derived type writer ?
!    Have a default writer ? (just like decomp_main for the derived type decomp_info)
!    Keep type-bound procedures or avoid them ?
!    Move the subroutines gen_iodir_name, coarse_extents and plane_extents inside writer.f90 to hide them
!
! Simplification
!    Remove icoarse and the special cases icoarse > 0
!       Restore later
!          Using dedicated decomp_info objects ?
!          Combined with a small interpolation module allowing downsampling ?
!    Remove the reduced precision
!       Restore it later
!          Avoid memory allocation inside 2decomp during IO
!          Avoid transpose operations inside 2decomp during IO
!    d2d_writer_family_register_var
!       The interface can be improved
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
!       Make gen_iodir_name private
!    Adapt the examples accordingly
!
! Experiment new IO template ? Example for writing one 3D array :
!    4 public write_one subroutines
!       single precision real
!       double precision real
!       single precision complex
!       double precision complex
!       example :
!          write_one_sreal(ipencil, var, varname, writer, opt_decomp)
!             real(0._real32), dimension(:,:,:), intent(in) :: var
!             call write_one_generic(ipencil, varname, writer, opt_decomp, svar = var)
!          write_one_dreal(ipencil, var, varname, writer, opt_decomp)
!             real(0._real64), dimension(:,:,:), intent(in) :: var
!             call write_one_generic(ipencil, varname, writer, opt_decomp, dvar = var)
!    1 private write_one subroutine to do the IO
!       4 optional arguments (SP real, DP real, SP cplx, DP cplx)
!       example :
!          write_one_generic(ipencil, varname, writer, opt_decomp, svar, dvar, cvar, cdvar)
!             real(0._real32), dimension(:,:,:), optional, intent(in) :: svar
!             real(0._real64), dimension(:,:,:), optional, intent(in) :: dvar
!             complex(0._real32), dimension(:,:,:), optional, intent(in) :: cvar
!             complex(0._real64), dimension(:,:,:), optional, intent(in) :: cdvar
!
!    Remarks :
!       no need to include a file, everything happens in one subroutine
!       should allow reduced precision operations easily
!

submodule(decomp_2d_io) decomp_2d_writer

   use decomp_2d_constants
#ifdef ADIOS2
   use adios2
#endif
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

contains

   !
   ! Initialize a new family of writers (default type)
   !
   subroutine d2d_writer_family_init(family, label)

      implicit none

      class(d2d_writer_family), intent(inout) :: family
      character(len=*), intent(in) :: label

#ifdef ADIOS2
      call family%adios2_init(label)
#else
      call family%mpi_init(label)
#endif

   end subroutine d2d_writer_family_init

   !
   ! Initialize a new family of MPI writers
   !
   subroutine d2d_writer_family_mpi_init(family, label)

      implicit none

      class(d2d_writer_family), intent(inout) :: family
      character(len=*), intent(in) :: label

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_family_mpi_init")
#endif

      ! Safety check
      if (family%type /= decomp_2d_writer_none) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Family was not cleared "//label)
      end if

      family%type = decomp_2d_writer_mpi
      family%label = label

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_family_mpi_init")
#endif

   end subroutine d2d_writer_family_mpi_init

   !
   ! Initialize a new family of ADIOS2 writers
   !
   subroutine d2d_writer_family_adios2_init(family, label)

      implicit none

      class(d2d_writer_family), intent(inout) :: family
      character(len=*), intent(in) :: label

#ifdef ADIOS2
      integer :: ierror
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_family_adios2_init")
#endif

      ! Safety check
      if (family%type /= decomp_2d_writer_none) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Family was not cleared "//label)
      end if

#ifndef ADIOS2
      call decomp_2d_abort(__FILE__, __LINE__, -1, "ADIOS2 is not available")
#else
      family%type = decomp_2d_writer_adios2
      family%label = label

      ! Advanced API
      ! The external code can set its own object of type adios2_adios before calling init
      if (.not.associated(family%adios)) then
         family%adios => adios
      endif

      if (family%adios%valid) then
         call adios2_declare_io(family%io, family%adios, label, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "adios2_declare_io "//label)
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, -1, "adios object not valid")
      end if
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_family_adios2_init")
#endif

   end subroutine d2d_writer_family_adios2_init

   !
   ! Clear the given writer
   !
   subroutine d2d_writer_family_fin(family)

      implicit none

      class(d2d_writer_family), intent(inout) :: family

#ifdef ADIOS2
      integer :: ierror
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_family_fin")
#endif

      ! Safety check
      if (family%type == decomp_2d_writer_none) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, "Family was already cleared.")
         return
      end if

#ifdef ADIOS2
      if (family%type == decomp_2d_writer_adios2) then
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
         nullify(family%adios)
      end if
#endif

      family%type = decomp_2d_writer_none
      deallocate (family%label)

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_family_fin")
#endif

   end subroutine d2d_writer_family_fin

   !
   ! Register one variable for the given family of writers
   !    1 <= ipencil <= 3
   !    0 <= iplane <= 3
   !
   !    Warning, iplane > 0 is not compatible with icoarse > 0
   !
   module subroutine d2d_writer_family_register_var(family, varname, ipencil, iplane, &
                                                    type, opt_decomp, opt_nplanes)

      implicit none

      class(d2d_writer_family), intent(inout) :: family
      character(len=*), intent(in) :: varname
      integer, intent(in) :: ipencil ! (x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(in) :: iplane
      integer, intent(in) :: type
      type(decomp_info), intent(in), optional :: opt_decomp
      integer, intent(in), optional :: opt_nplanes

#ifdef ADIOS2
      integer :: nplanes
      type(adios2_variable) :: var_handle
      integer, dimension(3) :: sizes, subsizes, starts
      type(adios2_variable) :: var_handle
      integer, parameter :: ndims = 3
      logical, parameter :: adios2_constant_dims = .true.
      integer :: data_type
      integer :: ierror
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_family_register_var")
#endif

#ifdef ADIOS2
      if (family%type == decomp_2d_writer_adios2) then

         ! Prepare to register
         if (iplane == 0) then
            if (present(opt_decomp)) then
               call coarse_extents(ipencil, 0, sizes, subsizes, starts, opt_decomp)
            else
               call coarse_extents(ipencil, 0, sizes, subsizes, starts)
            end if
         else
            if (present(opt_nplanes)) then
               nplanes = opt_nplanes
            else
               nplanes = 1
            end if
            if (present(opt_decomp)) then
               call plane_extents(sizes, subsizes, starts, iplane, opt_decomp, opt_nplanes=nplanes)
            else
               call plane_extents(sizes, subsizes, starts, iplane, opt_nplanes=nplanes)
            end if
         end if

         ! Register
         if (family%io%valid) then
            call adios2_inquire_variable(var_handle, family%io, varname, ierror)
            if (ierror /= 0) then
               call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                                    "adios2_inquire_variable "//varname)
            end if
            if (.not. var_handle%valid) then

               ! New variable
               if (nrank == 0) then
                  print *, "Registering variable for IO: ", varname
               end if

               ! Need to set the ADIOS2 data type
               if (type == kind(0._real64)) then
                  ! Double
                  data_type = adios2_type_dp
               else if (type == kind(0._real32)) then
                  ! Single
                  data_type = adios2_type_real
               else
                  ! This could be expanded
                  ! adios2_type_complex
                  ! adios2_type_complex_dp
                  ! adios2_type_integer1
                  ! adios2_type_integer2
                  ! adios2_type_integer4
                  ! adios2_type_integer8
                  !
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
            end if
         else
            call decomp_2d_abort(__FILE__, __LINE__, -1, &
                                 "trying to register variable with invalid IO!")
         end if
      end if
#else
      associate (fm => family, vr => varname, pc => ipencil, &
                 pl => iplane, tp => type, od => opt_decomp, on => opt_nplanes)
      end associate
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_family_register_var")
#endif

   end subroutine d2d_writer_family_register_var

   !
   ! Open the given writer
   !
   module subroutine d2d_writer_open(writer, family, io_dir, mode)

      implicit none

      class(d2d_writer), intent(inout) :: writer
      type(d2d_writer_family), target, intent(in) :: family
      character(len=*), intent(in) :: io_dir
      integer, intent(in) :: mode

      integer :: access_mode, ierror

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_open")
#endif

      ! Safety check
      if (writer%is_open) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, &
                              "Writer was not closed "//writer%label)
      end if

      ! Prepare the writer
      writer%family => family
      writer%is_open = .true.
      writer%label = io_dir

      ! Prepare to open
      if (family%type == decomp_2d_writer_mpi) then
         ! MPI writer, set initial displacement and access mode
         writer%disp = 0_MPI_OFFSET_KIND
         if (mode == decomp_2d_write_mode .or. mode == decomp_2d_append_mode) then
            access_mode = MPI_MODE_CREATE + MPI_MODE_WRONLY
         else if (mode == decomp_2d_read_mode) then
            access_mode = MPI_MODE_RDONLY
         else
            call decomp_2d_abort(__FILE__, __LINE__, mode, "Invalid value for mode")
         end if

      else if (family%type == decomp_2d_writer_adios2) then
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
         call decomp_2d_abort(__FILE__, __LINE__, family%type, "Invalid family")
      end if

      ! Open IO
      if (family%type == decomp_2d_writer_mpi) then
         call MPI_FILE_OPEN(decomp_2d_comm, io_dir, &
                            access_mode, MPI_INFO_NULL, &
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

      else if (family%type == decomp_2d_writer_adios2) then
#ifdef ADIOS2
         if (family%io%valid) then
            call adios2_open(writer%engine, family%io, &
                             trim(gen_iodir_name(io_dir, family%label)), &
                             access_mode, ierror)
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

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_open")
#endif

   end subroutine d2d_writer_open

   !
   ! Start IO for the given writer
   !
   module subroutine d2d_writer_start(writer)

      implicit none

      class(d2d_writer), intent(inout) :: writer

#ifdef ADIOS2
      integer :: ierror
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_start")
#endif

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

      ! Start the writer
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
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_start")
#endif

   end subroutine d2d_writer_start

   !
   ! Open the given writer and start IO
   !
   module subroutine d2d_writer_open_start(writer, family, io_dir, mode)

      implicit none

      class(d2d_writer), intent(inout) :: writer
      type(d2d_writer_family), target, intent(in) :: family
      character(len=*), intent(in) :: io_dir
      integer, intent(in) :: mode

      call writer%open(family, io_dir, mode)
      call writer%start()

   end subroutine d2d_writer_open_start

   !
   ! End IO for the given writer
   !
   module subroutine d2d_writer_end(writer)

      implicit none

      class(d2d_writer), intent(inout) :: writer

#ifdef ADIOS2
      integer :: ierror
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_end")
#endif

#ifdef ADIOS2
      ! Safety check
      if (.not. writer%is_active) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                "Writer was not started "//writer%label)
         return
      end if

      ! Tag the writer
      writer%is_active = .false.

      ! FIXME
      ! Currently only one engine per writer
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
#else
      associate (wrt => writer)
      end associate
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_end")
#endif

   end subroutine d2d_writer_end

   !
   ! Close the given writer
   !
   module subroutine d2d_writer_close(writer)

      implicit none

      class(d2d_writer), intent(inout) :: writer

      integer :: ierror

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_writer_close")
#endif

      ! Safety check
      if (.not. writer%is_open) then
         call decomp_2d_warning(__FILE__, __LINE__, -1, &
                                "Writer was not opened "//writer%label)
         return
      end if

      ! Close the writer (type was checked at open stage)
      if (writer%family%type == decomp_2d_writer_mpi) then
         call MPI_FILE_CLOSE(writer%fh, ierror)
         if (ierror /= 0) then
            call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_FILE_CLOSE")
         end if
      else if (writer%family%type == decomp_2d_writer_adios2) then
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

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_writer_close")
#endif

   end subroutine d2d_writer_close

   !
   ! End IO and close the given writer
   !
   module subroutine d2d_writer_end_close(writer)

      implicit none

      class(d2d_writer), intent(inout) :: writer

      call writer%end()
      call writer%close()

   end subroutine d2d_writer_end_close

end submodule decomp_2d_writer
