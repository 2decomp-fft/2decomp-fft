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
! This is the IO module for families of readers / writers
!    => Initialize / finalize a family of readers / writers
!    => Register a variable in a family of readers / writers
!

module decomp_2d_io_family

   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d
#ifdef ADIOS2
   use adios2
#endif
   use, intrinsic :: iso_fortran_env, only: real32, real64

   implicit none

#ifdef ADIOS2
   type(adios2_adios), pointer :: default_adios => null()
#endif

   !
   ! derived type to store info for a family of readers / writers
   !
   type, public :: d2d_io_family
      integer :: type = decomp_2d_io_none                     ! Type of the writer
      character(:), allocatable :: label                      ! Label of the writer
#ifdef ADIOS2
      type(adios2_adios), pointer :: adios                    ! adios2 only
      type(adios2_io) :: io                                   ! adios2 only
#endif
   contains
      procedure :: init => d2d_io_family_init                 ! ADIOS2 writer if possible, MPI otherwise
      procedure :: mpi_init => d2d_io_family_mpi_init         ! Force MPI writer
      procedure :: adios2_init => d2d_io_family_adios2_init   ! Force ADIOS2 writer
      procedure :: fin => d2d_io_family_fin                   ! Clear the writer
      procedure :: register_var => d2d_io_family_register_var ! Register a variable
   end type d2d_io_family

   private

#ifdef ADIOS2
   public :: decomp_2d_io_family_get_default_adios, &
             decomp_2d_io_family_set_default_adios
#endif

contains

   !
   ! Initialize a new family of readers / writers (default type)
   !
   subroutine d2d_io_family_init(family, label)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: label

#ifdef ADIOS2
      call family%adios2_init(label)
#else
      call family%mpi_init(label)
#endif

   end subroutine d2d_io_family_init

   !
   ! Initialize a new family of MPI readers / writers
   !
   subroutine d2d_io_family_mpi_init(family, label)

      implicit none

      class(d2d_io_family), intent(inout) :: family
      character(len=*), intent(in) :: label

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_mpi_init")
#endif

      ! Safety check
      if (family%type /= decomp_2d_io_none) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Family was not cleared "//label)
      end if

      family%type = decomp_2d_io_mpi
      family%label = label

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_mpi_init")
#endif

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
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_adios2_init")
#endif

      ! Safety check
      if (family%type /= decomp_2d_io_none) then
         call decomp_2d_abort(__FILE__, __LINE__, -1, "Family was not cleared "//label)
      end if

#ifndef ADIOS2
      call decomp_2d_abort(__FILE__, __LINE__, -1, "ADIOS2 is not available")
#else
      family%type = decomp_2d_io_adios2
      family%label = label

      ! Advanced API
      ! The external code can set its own object of type adios2_adios before calling init
      if (.not. associated(family%adios)) then
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

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_adios2_init")
#endif

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

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_fin")
#endif

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

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_fin")
#endif

   end subroutine d2d_io_family_fin

   !
   ! Register one variable for the given family of readers / writers
   !    1 <= ipencil <= 3
   !    0 <= iplane <= 3
   !
   module subroutine d2d_io_family_register_var(family, &
                                                varname, &
                                                ipencil, &
                                                iplane, &
                                                type, &
                                                opt_decomp, &
                                                opt_nplanes)

      implicit none

      class(d2d_io_family), intent(inout) :: family
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
      integer, parameter :: ndims = 3
      logical, parameter :: adios2_constant_dims = .true.
      integer :: data_type
      integer :: ierror
#endif

#ifdef PROFILER
      if (decomp_profiler_io) call decomp_profiler_start("d2d_io_family_register_var")
#endif

#ifdef ADIOS2
      if (family%type == decomp_2d_io_adios2) then

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
      if (decomp_profiler_io) call decomp_profiler_end("d2d_io_family_register_var")
#endif

   end subroutine d2d_io_family_register_var

#ifdef ADIOS2
   !
   ! Return a pointer to the default adios2_adios object
   !
   function decomp_2d_io_family_get_default_adios()

      implicit none

      type(adios2_adios), pointer :: decomp_2d_io_family_get_default_adios

      if (.not. associated(default_adios)) then
         call decomp_2d_warning(0, "Warning, default adios2_adios is null()")
      endif

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
         nullify(default_adios)
      endif

   end subroutine decomp_2d_io_family_set_default_adios
#endif

   subroutine coarse_extents(ipencil, icoarse, sizes, subsizes, starts, opt_decomp)

      implicit none

      integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
      integer, intent(IN) :: icoarse !(nstat=1; nvisu=2)
      type(decomp_info), intent(in), optional :: opt_decomp

      integer, dimension(3) :: sizes, subsizes, starts
      type(decomp_info) :: decomp

      if ((icoarse < 0) .or. (icoarse > 2)) then
         call decomp_2d_abort(__FILE__, __LINE__, icoarse, "Error invalid value of icoarse")
      end if
      if ((ipencil < 1) .or. (ipencil > 3)) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Error invalid value of ipencil ")
      end if

      if (icoarse == 0) then
         ! Use full fields

         if (present(opt_decomp)) then
            decomp = opt_decomp
         else
            call get_decomp_info(decomp)
         end if

         sizes(1) = decomp%xsz(1)
         sizes(2) = decomp%ysz(2)
         sizes(3) = decomp%zsz(3)

         if (ipencil == 1) then
            subsizes(1:3) = decomp%xsz(1:3)
            starts(1:3) = decomp%xst(1:3) - 1
         elseif (ipencil == 2) then
            subsizes(1:3) = decomp%ysz(1:3)
            starts(1:3) = decomp%yst(1:3) - 1
         elseif (ipencil == 3) then
            subsizes(1:3) = decomp%zsz(1:3)
            starts(1:3) = decomp%zst(1:3) - 1
         else
            call decomp_2d_abort(-1, "IO/coarse_extents : Wrong value for ipencil")
         end if
      elseif (icoarse == 1) then
         sizes(1) = xszS(1)
         sizes(2) = yszS(2)
         sizes(3) = zszS(3)

         if (ipencil == 1) then
            subsizes(1) = xszS(1)
            subsizes(2) = xszS(2)
            subsizes(3) = xszS(3)
            starts(1) = xstS(1) - 1  ! 0-based index
            starts(2) = xstS(2) - 1
            starts(3) = xstS(3) - 1
         else if (ipencil == 2) then
            subsizes(1) = yszS(1)
            subsizes(2) = yszS(2)
            subsizes(3) = yszS(3)
            starts(1) = ystS(1) - 1
            starts(2) = ystS(2) - 1
            starts(3) = ystS(3) - 1
         else if (ipencil == 3) then
            subsizes(1) = zszS(1)
            subsizes(2) = zszS(2)
            subsizes(3) = zszS(3)
            starts(1) = zstS(1) - 1
            starts(2) = zstS(2) - 1
            starts(3) = zstS(3) - 1
         else
            call decomp_2d_abort(-1, "IO/coarse_extents : Wrong value for ipencil")
         end if
      elseif (icoarse == 2) then
         sizes(1) = xszV(1)
         sizes(2) = yszV(2)
         sizes(3) = zszV(3)

         if (ipencil == 1) then
            subsizes(1) = xszV(1)
            subsizes(2) = xszV(2)
            subsizes(3) = xszV(3)
            starts(1) = xstV(1) - 1  ! 0-based index
            starts(2) = xstV(2) - 1
            starts(3) = xstV(3) - 1
         else if (ipencil == 2) then
            subsizes(1) = yszV(1)
            subsizes(2) = yszV(2)
            subsizes(3) = yszV(3)
            starts(1) = ystV(1) - 1
            starts(2) = ystV(2) - 1
            starts(3) = ystV(3) - 1
         else if (ipencil == 3) then
            subsizes(1) = zszV(1)
            subsizes(2) = zszV(2)
            subsizes(3) = zszV(3)
            starts(1) = zstV(1) - 1
            starts(2) = zstV(2) - 1
            starts(3) = zstV(3) - 1
         else
            call decomp_2d_abort(-1, "IO/coarse_extents : Wrong value for ipencil")
         end if
      end if

   end subroutine coarse_extents

   subroutine plane_extents(sizes, subsizes, starts, iplane, opt_decomp, opt_nplanes)

      integer, intent(in) :: iplane
      type(decomp_info), intent(in), optional :: opt_decomp
      integer, intent(in), optional :: opt_nplanes

      integer, dimension(3), intent(out) :: sizes, subsizes, starts

      integer :: nplanes
      type(decomp_info) :: decomp

      if (present(opt_decomp)) then
         decomp = opt_decomp
      else
         call get_decomp_info(decomp)
      end if

      if (present(opt_nplanes)) then
         nplanes = opt_nplanes
      else
         nplanes = 1
      end if

      if (iplane == 1) then
         sizes(1) = nplanes
         sizes(2) = decomp%ysz(2)
         sizes(3) = decomp%zsz(3)
         subsizes(1) = nplanes
         subsizes(2) = decomp%xsz(2)
         subsizes(3) = decomp%xsz(3)
         starts(1) = 0
         starts(2) = decomp%xst(2) - 1
         starts(3) = decomp%xst(3) - 1
      else if (iplane == 2) then
         sizes(1) = decomp%xsz(1)
         sizes(2) = nplanes
         sizes(3) = decomp%zsz(3)
         subsizes(1) = decomp%ysz(1)
         subsizes(2) = nplanes
         subsizes(3) = decomp%ysz(3)
         starts(1) = decomp%yst(1) - 1
         starts(2) = 0
         starts(3) = decomp%yst(3) - 1
      else if (iplane == 3) then
         sizes(1) = decomp%xsz(1)
         sizes(2) = decomp%ysz(2)
         sizes(3) = nplanes
         subsizes(1) = decomp%zsz(1)
         subsizes(2) = decomp%zsz(2)
         subsizes(3) = nplanes
         starts(1) = decomp%zst(1) - 1
         starts(2) = decomp%zst(2) - 1
         starts(3) = 0
      else
         print *, "Can't work with plane ", iplane
         stop
      end if

   end subroutine plane_extents

end module decomp_2d_io_family
