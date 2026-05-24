!! SPDX-License-Identifier: BSD-3-Clause

! This is the FFTW (version 3.x) implementation of the FFT library

module decomp_2d_fft

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_profiler
   use iso_c_binding
   use cudafor
   use cufft
   use m_decomp_pool

   implicit none

   private        ! Make everything private unless declared public

   ! engine-specific global variables
   ! integer, save :: plan_type = FFTW_MEASURE

   ! FFTW plans
   ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
   ! For c2c transforms:
   !     use plan(-1,j) for  forward transform;
   !     use plan( 1,j) for backward transform;
   ! For r2c/c2r transforms:
   !     use plan(0,j) for r2c transforms;
   !     use plan(2,j) for c2r transforms;
   integer*4, contiguous, pointer, save :: plan(:, :) => null()
   complex*8, device, contiguous, pointer, dimension(:) :: cufft_workspace => null()

   integer, parameter, public :: D2D_FFT_BACKEND = D2D_FFT_BACKEND_CUFFT

   ! Derived type with all the quantities needed to perform FFT
   type decomp_2d_fft_engine
      ! Engine-specific stuff
      integer*4, private :: plan(-1:2, 3)
      complex*8, private, device, allocatable, dimension(:) :: cufft_workspace
      ! All the engines have this
      integer, private :: format
      logical, private :: initialised = .false.
      integer, private :: nx_fft, ny_fft, nz_fft
      type(decomp_info), pointer, public :: ph => null()
      type(decomp_info), private :: ph_target ! ph => ph_target or ph => decomp_main
      type(decomp_info), public :: sp
      complex(mytype), allocatable, private :: wk2_c2c(:, :, :)
      complex(mytype), contiguous, pointer, private :: wk2_r2c(:, :, :) => null()
      complex(mytype), allocatable, private :: wk13(:, :, :)
      logical, private :: inplace
      logical, private :: skip_x_c2c, skip_y_c2c, skip_z_c2c
   contains
      procedure, public :: init => decomp_2d_fft_engine_init
      procedure, public :: fin => decomp_2d_fft_engine_fin
      procedure, public :: use_it => decomp_2d_fft_engine_use_it
      generic, public :: fft => c2c, r2c, c2r
      procedure, private :: c2c => decomp_2d_fft_engine_fft_c2c
      procedure, private :: r2c => decomp_2d_fft_engine_fft_r2c
      procedure, private :: c2r => decomp_2d_fft_engine_fft_c2r
   end type decomp_2d_fft_engine

   ! Workspace to store the intermediate Y-pencil data
   complex(mytype), contiguous, pointer, dimension(:, :, :) :: wk2_r2c => null(), &
                                                               wk2_c2c => null(), &
                                                               wk13 => null()

   ! common code used for all engines, including global variables,
   ! generic interface definitions and several subroutines
#include "fft_common.f90"

   ! Return a cuFFT plan for multiple 1D FFTs in X direction: C2C case
   subroutine c2c_1m_x_plan(plan1, decomp, cufft_type, worksize)

      implicit none

      integer*4, intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: cufft_type

      integer :: istat
      integer(int_ptr_kind()), intent(out) :: worksize
      integer, pointer :: null_fptr
      call c_f_pointer(c_null_ptr, null_fptr)

      istat = cufftCreate(plan1)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftCreate")
      istat = cufftSetAutoAllocation(plan1, 0)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetAutoAllocation")
      istat = cufftMakePlanMany(plan1, 1, decomp%xsz(1), &
                                decomp%xsz(1), 1, decomp%xsz(1), &
                                decomp%xsz(1), 1, decomp%xsz(1), &
                                cufft_type, decomp%xsz(2) * decomp%xsz(3), worksize)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftMakePlanMany")

   end subroutine c2c_1m_x_plan

   ! Return a cuFFT plan for multiple 1D FFTs in Y direction: C2C case
   subroutine c2c_1m_y_plan(plan1, decomp, cufft_type, worksize)

      implicit none

      integer*4, intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: cufft_type

      ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
      ! done one Z-plane at a time. So plan for 2D data sets here.
      integer :: istat
      integer(int_ptr_kind()), intent(out) :: worksize
      integer, pointer :: null_fptr
      call c_f_pointer(c_null_ptr, null_fptr)

      istat = cufftCreate(plan1)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftCreate")
      istat = cufftSetAutoAllocation(plan1, 0)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetAutoAllocation")
      istat = cufftMakePlanMany(plan1, 1, decomp%ysz(2), &
                                decomp%ysz(2), decomp%ysz(1), 1, &
                                decomp%ysz(2), decomp%ysz(1), 1, &
                                cufft_type, decomp%ysz(1), worksize)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftMakePlanMany")

   end subroutine c2c_1m_y_plan

   ! Return a cuFFT plan for multiple 1D FFTs in Z direction: C2C case
   subroutine c2c_1m_z_plan(plan1, decomp, cufft_type, worksize)

      implicit none

      integer*4, intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp
      integer, intent(IN) :: cufft_type

      integer :: istat
      integer(int_ptr_kind()), intent(out) :: worksize
      integer, pointer :: null_fptr
      call c_f_pointer(c_null_ptr, null_fptr)

      istat = cufftCreate(plan1)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftCreate")
      istat = cufftSetAutoAllocation(plan1, 0)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetAutoAllocation")
      istat = cufftMakePlanMany(plan1, 1, decomp%zsz(3), &
                                decomp%zsz(3), decomp%zsz(1) * decomp%zsz(2), 1, &
                                decomp%zsz(3), decomp%zsz(1) * decomp%zsz(2), 1, &
                                cufft_type, decomp%zsz(1) * decomp%zsz(2), worksize)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftMakePlanMany")

   end subroutine c2c_1m_z_plan

   ! Return a cuFFT plan for multiple 1D FFTs in X direction: R2C case
   subroutine r2c_1m_x_plan(plan1, decomp_ph, decomp_sp, cufft_type, worksize)

      implicit none

      integer*4, intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      integer, intent(IN) :: cufft_type

      integer :: istat
      integer(int_ptr_kind()), intent(out) :: worksize
      integer, pointer :: null_fptr
      call c_f_pointer(c_null_ptr, null_fptr)

      istat = cufftCreate(plan1)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftCreate")
      istat = cufftSetAutoAllocation(plan1, 0)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetAutoAllocation")
      istat = cufftMakePlanMany(plan1, 1, decomp_ph%xsz(1), &
                                decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
                                decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
                                cufft_type, decomp_ph%xsz(2) * decomp_ph%xsz(3), worksize)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftMakePlanMany")

   end subroutine r2c_1m_x_plan

   ! Return a cuFFT plan for multiple 1D FFTs in X direction: C2R case
   subroutine c2r_1m_x_plan(plan1, decomp_sp, decomp_ph, cufft_type, worksize)

      implicit none

      integer*4, intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      integer, intent(IN) :: cufft_type

      integer :: istat
      integer(int_ptr_kind()), intent(out) :: worksize
      integer, pointer :: null_fptr
      call c_f_pointer(c_null_ptr, null_fptr)

      istat = cufftCreate(plan1)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftCreate")
      istat = cufftSetAutoAllocation(plan1, 0)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetAutoAllocation")
      istat = cufftMakePlanMany(plan1, 1, decomp_ph%xsz(1), &
                                decomp_sp%xsz(1), 1, decomp_sp%xsz(1), &
                                decomp_ph%xsz(1), 1, decomp_ph%xsz(1), &
                                cufft_type, decomp_ph%xsz(2) * decomp_ph%xsz(3), worksize)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftMakePlanMany")

   end subroutine c2r_1m_x_plan

   ! Return a cuFFT plan for multiple 1D FFTs in X direction: R2C case
   subroutine r2c_1m_z_plan(plan1, decomp_ph, decomp_sp, cufft_type, worksize)

      implicit none

      integer*4, intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      integer, intent(IN) :: cufft_type

      integer :: istat
      integer(int_ptr_kind()), intent(out) :: worksize
      integer, pointer :: null_fptr
      call c_f_pointer(c_null_ptr, null_fptr)

      istat = cufftCreate(plan1)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftCreate")
      istat = cufftSetAutoAllocation(plan1, 0)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetAutoAllocation")
      istat = cufftMakePlanMany(plan1, 1, decomp_ph%zsz(3), &
                                decomp_ph%zsz(3), decomp_ph%zsz(1) * decomp_ph%zsz(2), 1, &
                                decomp_sp%zsz(3), decomp_sp%zsz(1) * decomp_sp%zsz(2), 1, &
                                cufft_type, decomp_ph%zsz(1) * decomp_ph%zsz(2), worksize)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftMakePlanMany")

   end subroutine r2c_1m_z_plan

   ! Return a cuFFT plan for multiple 1D FFTs in X direction: C2R case
   subroutine c2r_1m_z_plan(plan1, decomp_sp, decomp_ph, cufft_type, worksize)

      implicit none

      integer*4, intent(OUT) :: plan1
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph
      TYPE(DECOMP_INFO), intent(IN) :: decomp_sp
      integer, intent(IN) :: cufft_type

      integer :: istat
      integer(int_ptr_kind()), intent(out) :: worksize
      integer, pointer :: null_fptr
      call c_f_pointer(c_null_ptr, null_fptr)

      istat = cufftCreate(plan1)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftCreate")
      istat = cufftSetAutoAllocation(plan1, 0)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetAutoAllocation")
      istat = cufftMakePlanMany(plan1, 1, decomp_ph%zsz(3), &
                                decomp_sp%zsz(3), decomp_sp%zsz(1) * decomp_sp%zsz(2), 1, &
                                decomp_ph%zsz(3), decomp_ph%zsz(1) * decomp_ph%zsz(2), 1, &
                                cufft_type, decomp_ph%zsz(1) * decomp_ph%zsz(2), worksize)
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftMakePlanMany")

   end subroutine c2r_1m_z_plan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time initialisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_fft_engine(engine)

      implicit none

      type(decomp_2d_fft_engine), target, intent(inout) :: engine

      integer(int_ptr_kind()) :: cufft_ws, ws
      integer :: i, j, istat

      !
      ! Allocate the workspace for intermediate y-pencil data
      ! The largest memory block needed is the one for c2c transforms
      !
      call alloc_y(engine%wk2_c2c, engine%ph)
      !
      ! A smaller memory block is needed for r2c and c2r transforms
      ! wk2_c2c and wk2_r2c start at the same memory location
      !
      !    Size of wk2_c2c : ph%ysz(1), ph%ysz(2), ph%ysz(3)
      !    Size of wk2_r2c : sp%ysz(1), sp%ysz(2), sp%ysz(3)
      !
      call c_f_pointer(c_loc(engine%wk2_c2c), engine%wk2_r2c, engine%sp%ysz)
      !
      ! Allocate the workspace for r2c and c2r transforms
      !
      ! wk13 can not be easily fused with wk2_*2c due to statements such as
      ! transpose_y_to_x(wk2_r2c, wk13, sp)
      ! transpose_y_to_z(wk2_r2c, wk13, sp)
      !
      if (engine%format == PHYSICAL_IN_X) then
         call alloc_x(engine%wk13, engine%sp)
      else if (engine%format == PHYSICAL_IN_Z) then
         call alloc_z(engine%wk13, engine%sp)
      end if

      call decomp_2d_fft_log("cuFFT")

      cufft_ws = 0
#ifdef DOUBLE_PREC
      if (format == PHYSICAL_IN_X) then
         ! For C2C transforms
         call c2c_1m_x_plan(plan(-1, 1), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(-1, 2), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(-1, 3), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(1, 3), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(1, 2), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(1, 1), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         ! For R2C/C2R tranforms
         call r2c_1m_x_plan(plan(0, 1), ph, sp, CUFFT_D2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(0, 2), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(0, 3), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(2, 3), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(2, 2), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2r_1m_x_plan(plan(2, 1), sp, ph, CUFFT_Z2D, ws)
         cufft_ws = max(cufft_ws, ws)

      else if (format == PHYSICAL_IN_Z) then

         ! For C2C transforms
         call c2c_1m_z_plan(plan(-1, 3), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(-1, 2), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(-1, 1), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(1, 1), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(1, 2), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(1, 3), ph, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)

         ! For R2C/C2R tranforms
         call r2c_1m_z_plan(plan(0, 3), ph, sp, CUFFT_D2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(0, 2), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(0, 1), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(2, 1), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(2, 2), sp, CUFFT_Z2Z, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2r_1m_z_plan(plan(2, 3), sp, ph, CUFFT_Z2D, ws)
         cufft_ws = max(cufft_ws, ws)

      end if
#else
      if (format == PHYSICAL_IN_X) then
         ! For C2C transforms
         call c2c_1m_x_plan(plan(-1, 1), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(-1, 2), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(-1, 3), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(1, 3), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(1, 2), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(1, 1), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         ! For R2C/C2R tranforms
         call r2c_1m_x_plan(plan(0, 1), ph, sp, CUFFT_R2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(0, 2), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(0, 3), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(2, 3), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(2, 2), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2r_1m_x_plan(plan(2, 1), sp, ph, CUFFT_C2R, ws)
         cufft_ws = max(cufft_ws, ws)

      else if (format == PHYSICAL_IN_Z) then

         ! For C2C transforms
         call c2c_1m_z_plan(plan(-1, 3), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(-1, 2), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(-1, 1), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(1, 1), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(1, 2), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_z_plan(plan(1, 3), ph, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)

         ! For R2C/C2R tranforms
         call r2c_1m_z_plan(plan(0, 3), ph, sp, CUFFT_R2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(0, 2), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(0, 1), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_x_plan(plan(2, 1), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2c_1m_y_plan(plan(2, 2), sp, CUFFT_C2C, ws)
         cufft_ws = max(cufft_ws, ws)
         call c2r_1m_z_plan(plan(2, 3), sp, ph, CUFFT_C2R, ws)
         cufft_ws = max(cufft_ws, ws)

      end if
#endif
      cufft_ws = cufft_ws / sizeof(1._mytype)
      allocate (cufft_workspace(cufft_ws))
      do j = 1, 3
         do i = -1, 2
            istat = cufftSetWorkArea(plan(i, j), cufft_workspace)
            if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftSetWorkArea")
         end do
      end do

   end subroutine init_fft_engine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time finalisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine finalize_fft_engine(engine)

      implicit none

      type(decomp_2d_fft_engine), optional :: engine

      integer :: i, j, istat

      nullify (wk2_c2c)
      nullify (wk2_r2c)
      nullify (wk13)

      if (present(engine)) then

         if (allocated(engine%wk2_c2c)) deallocate (engine%wk2_c2c)
         if (associated(engine%wk2_r2c)) nullify (engine%wk2_r2c)
         if (allocated(engine%wk13)) deallocate (engine%wk13)

         do j = 1, 3
            do i = -1, 2
               istat = cufftDestroy(engine%plan(i, j))
               if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftDestroy")
            end do
         end do

      else

         nullify (plan)
         nullify (cufft_workspace)

      end if

   end subroutine finalize_fft_engine

   ! Use engine-specific stuff
   subroutine use_fft_engine(engine)

      implicit none

      type(decomp_2d_fft_engine), target, intent(in) :: engine

      wk2_c2c => engine%wk2_c2c
      wk2_r2c => engine%wk2_r2c
      wk13 => engine%wk13
      plan => engine%plan
      cufft_workspace => engine%cufft_workspace

   end subroutine use_fft_engine

   ! Following routines calculate multiple one-dimensional FFTs to form
   ! the basis of three-dimensional FFTs.

   ! c2c transform, multiple 1D FFTs in x direction
   subroutine c2c_1m_x(inout, isign, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      integer, intent(IN) :: isign
      integer*4, intent(IN) :: plan1

      integer :: istat

      if (skip_x_c2c) return

#ifdef DOUBLE_PREC
      !$acc host_data use_device(inout)
      istat = cufftExecZ2Z(plan1, inout, inout, isign)
      !$acc end host_data
#else
      !$acc host_data use_device(inout)
      istat = cufftExecC2C(plan1, inout, inout, isign)
      !$acc end host_data
#endif
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftExecC2C/Z2Z")

   end subroutine c2c_1m_x

   ! c2c transform, multiple 1D FFTs in y direction
   subroutine c2c_1m_y(inout, isign, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      integer, intent(IN) :: isign
      integer*4, intent(IN) :: plan1

      integer :: s3, k, istat

      if (skip_y_c2c) return

      ! transform on one Z-plane at a time
      s3 = size(inout, 3)
      do k = 1, s3
#ifdef DOUBLE_PREC
         !$acc host_data use_device(inout)
         istat = cufftExecZ2Z(plan1, inout(:, :, k), inout(:, :, k), isign)
         !$acc end host_data
#else
         !$acc host_data use_device(inout)
         istat = cufftExecC2C(plan1, inout(:, :, k), inout(:, :, k), isign)
         !$acc end host_data
#endif
         if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftExecC2C/Z2Z")
      end do

   end subroutine c2c_1m_y

   ! c2c transform, multiple 1D FFTs in z direction
   subroutine c2c_1m_z(inout, isign, plan1)

      implicit none

      complex(mytype), dimension(:, :, :), intent(INOUT) :: inout
      integer, intent(IN) :: isign
      integer*4, intent(IN) :: plan1

      integer :: istat

      if (skip_z_c2c) return

#ifdef DOUBLE_PREC
      !$acc host_data use_device(inout)
      istat = cufftExecZ2Z(plan1, inout, inout, isign)
      !$acc end host_data
#else
      !$acc host_data use_device(inout)
      istat = cufftExecC2C(plan1, inout, inout, isign)
      !$acc end host_data
#endif
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftExecC2C/Z2Z")

   end subroutine c2c_1m_z

   ! r2c transform, multiple 1D FFTs in x direction
   subroutine r2c_1m_x(input, output)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN)  ::  input
      complex(mytype), dimension(:, :, :), intent(OUT) :: output
      integer :: istat

      if (skip_x_c2c) call decomp_2d_warning(__FILE__, __LINE__, 1, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      !$acc host_data use_device(input,output)
      istat = cufftExecD2Z(plan(0, 1), input, output)
      !$acc end host_data
#else
      !$acc host_data use_device(input,output)
      istat = cufftExecR2C(plan(0, 1), input, output)
      !$acc end host_data
#endif
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftExecR2C/D2Z")

   end subroutine r2c_1m_x

   ! r2c transform, multiple 1D FFTs in z direction
   subroutine r2c_1m_z(input, output)

      implicit none

      real(mytype), dimension(:, :, :), intent(IN)     ::  input
      complex(mytype), dimension(:, :, :), intent(OUT) :: output

      integer :: istat

      if (skip_z_c2c) call decomp_2d_warning(__FILE__, __LINE__, 2, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      !$acc host_data use_device(input,output)
      istat = cufftExecD2Z(plan(0, 3), input, output)
      !$acc end host_data
#else
      !$acc host_data use_device(input,output)
      istat = cufftExecR2C(plan(0, 3), input, output)
      !$acc end host_data
#endif
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftExecR2C/D2Z")

   end subroutine r2c_1m_z

   ! c2r transform, multiple 1D FFTs in x direction
   subroutine c2r_1m_x(input, output)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN)  ::  input
      real(mytype), dimension(:, :, :), intent(OUT)    :: output

      integer :: istat

      if (skip_x_c2c) call decomp_2d_warning(__FILE__, __LINE__, 3, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      !$acc host_data use_device(input,output)
      istat = cufftExecZ2D(plan(2, 1), input, output)
      !$acc end host_data
#else
      !$acc host_data use_device(input,output)
      istat = cufftExecC2R(plan(2, 1), input, output)
      !$acc end host_data
#endif
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftExecC2R/Z2D")

   end subroutine c2r_1m_x

   ! c2r transform, multiple 1D FFTs in z direction
   subroutine c2r_1m_z(input, output)

      implicit none

      complex(mytype), dimension(:, :, :), intent(IN) :: input
      real(mytype), dimension(:, :, :), intent(OUT) :: output

      integer :: istat

      if (skip_z_c2c) call decomp_2d_warning(__FILE__, __LINE__, 4, &
                                             "r2c / c2r transform can not be skipped")

#ifdef DOUBLE_PREC
      !$acc host_data use_device(input,output)
      istat = cufftExecZ2D(plan(2, 3), input, output)
      !$acc end host_data
#else
      !$acc host_data use_device(input,output)
      istat = cufftExecC2R(plan(2, 3), input, output)
      !$acc end host_data
#endif
      if (istat /= 0) call decomp_2d_abort(__FILE__, __LINE__, istat, "cufftExecC2R/Z2D")

   end subroutine c2r_1m_z

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D FFT - complex to complex
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2c(in, out, isign)

      implicit none

      ! Arguments
      complex(mytype), dimension(:, :, :), intent(INOUT) :: in
      complex(mytype), dimension(:, :, :), intent(OUT) :: out
      integer, intent(IN) :: isign

      ! Local variables
      complex(mytype), allocatable, dimension(:, :, :) :: wk1

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2c")

      !$acc data create(wk2_c2c) present(in,out)

      if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_FORWARD .OR. &
          format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_BACKWARD) then

         ! ===== 1D FFTs in X =====
         if (inplace) then
            call c2c_1m_x(in, isign, plan(isign, 1))
         else
            call alloc_x(wk1, ph)
            !$acc enter data create(wk1) async
            !$acc wait
            !$acc kernels default(present)
            wk1(:, :, :) = in(:, :, :)
            !$acc end kernels
            call c2c_1m_x(wk1, isign, plan(isign, 1))
         end if

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1)==1) then
            ! Single rank or slab : input available in X and Y
            if (inplace) then
               call c2c_1m_y(in, isign, plan(isign, 2))
            else
               call c2c_1m_y(wk1, isign, plan(isign, 2))
            end if
         else
            if (dims(2)==1) then
               ! Slab : output available in Y and Z
               if (inplace) then
                  call transpose_x_to_y(in, out, ph)
               else
                  call transpose_x_to_y(wk1, out, ph)
               end if
               call c2c_1m_y(out, isign, plan(isign, 2))
            else
               ! Pencil : use Y buffer
               if (inplace) then
                  call transpose_x_to_y(in, wk2_c2c, ph)
               else
                  call transpose_x_to_y(wk1, wk2_c2c, ph)
               end if
               call c2c_1m_y(wk2_c2c, isign, plan(isign, 2))
            end if
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1)==1) then
            ! Single rank : use transpose to copy input inside output
            ! Slab : input available in X and Y
            if (inplace) then
               call transpose_y_to_z(in, out, ph)
            else
               call transpose_y_to_z(wk1, out, ph)
            end if
         else if (dims(2) > 1) then
            ! Slab or pencil : transpose if needed
            call transpose_y_to_z(wk2_c2c, out, ph)
         end if
         call c2c_1m_z(out, isign, plan(isign, 3))

      else if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_BACKWARD &
               .OR. &
               format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_FORWARD) then

         ! ===== 1D FFTs in Z =====
         if (inplace) then
            call c2c_1m_z(in, isign, plan(isign, 3))
         else
            call alloc_z(wk1, ph)
            !$acc enter data create(wk1) async
            !$acc wait
            !$acc kernels default(present)
            wk1(:, :, :) = in(:, :, :)
            !$acc end kernels
            call c2c_1m_z(wk1, isign, plan(isign, 3))
         end if

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(2)==1) then
            ! Single rank or slab : input available in Y and Z
            if (inplace) then
               call c2c_1m_y(in, isign, plan(isign, 2))
            else
               call c2c_1m_y(wk1, isign, plan(isign, 2))
            end if
         else
            if (dims(1)==1) then
               ! Slab : output available in X and Y
               if (inplace) then
                  call transpose_z_to_y(in, out, ph)
               else
                  call transpose_z_to_y(wk1, out, ph)
               end if
               call c2c_1m_y(out, isign, plan(isign, 2))
            else
               ! Pencil : use Y buffer
               if (inplace) then
                  call transpose_z_to_y(in, wk2_c2c, ph)
               else
                  call transpose_z_to_y(wk1, wk2_c2c, ph) 
               end if 
               call c2c_1m_y(wk2_c2c, isign, plan(isign, 2))
            end if
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(2)==1) then
            ! Single rank : use transpose to copy input inside output
            ! Slab : input available in Y and Z
            if (inplace) then
               call transpose_y_to_x(in, out, ph)
            else
               call transpose_y_to_x(wk1, out, ph)
            end if
         else if (dims(1) > 1) then
            ! Slab or pencil : transpose if needed
            call transpose_y_to_x(wk2_c2c, out, ph)
         end if
         call c2c_1m_x(out, isign, plan(isign, 1))

      end if

      ! Free memory
      if (allocated(wk1)) then
         !$acc exit data delete(wk1) async
         !$acc wait
         deallocate (wk1)
      end if

      !$acc end data

      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2c")

   end subroutine fft_3d_c2c

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D forward FFT - real to complex
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_r2c(in_r, out_c)

      implicit none

      ! Arguments
      real(mytype), dimension(:, :, :), intent(IN) :: in_r
      complex(mytype), dimension(:, :, :), intent(OUT) :: out_c

      if (decomp_profiler_fft) call decomp_profiler_start("fft_r2c")

      !$acc data create(wk13,wk2_r2c) present(in_r,out_c)

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in X =====
         if (dims(1)==1 .and. dims(2)==1) then
            ! Single rank : output available in X, Y and Z
            call r2c_1m_x(in_r, out_c)
         else if (dims(1)==1) then
            ! Slab : Y buffer available in X
            call r2c_1m_x(in_r, wk2_r2c)
         else
            ! Default : use X buffer
            call r2c_1m_x(in_r, wk13)
         end if

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1)==1 .and. dims(2)==1) then
            ! Single rank : output available in X, Y and Z
            call c2c_1m_y(out_c, -1, plan(0, 2))
         else if (dims(2)==1) then
            ! Slab : output available in Y and Z
            call transpose_x_to_y(wk13, out_c, sp)
            call c2c_1m_y(out_c, -1, plan(0, 2))
         else
            ! Default : transpose if needed
            if (dims(1) > 1) call transpose_x_to_y(wk13, wk2_r2c, sp)
            call c2c_1m_y(wk2_r2c, -1, plan(0, 2))
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(2) > 1) call transpose_y_to_z(wk2_r2c, out_c, sp)
         call c2c_1m_z(out_c, -1, plan(0, 3))

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in Z =====
         if (dims(1)==1 .and. dims(2)==1) then
            ! Single rank : output available in X, Y and Z
            call r2c_1m_z(in_r, out_c)
         else if (dims(2)==1) then
            ! Slab : Y buffer available in Z
            call r2c_1m_z(in_r, wk2_r2c)
         else
            ! Default : use Z buffer
            call r2c_1m_z(in_r, wk13)
         end if

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(1)==1 .and. dims(2)==1) then
            ! Single rank : output available in X, Y and Z
            call c2c_1m_y(out_c, -1, plan(0, 2))
         else if (dims(1)==1) then
            ! Slab : output available in X and Y
            call transpose_z_to_y(wk13, out_c, sp)
            call c2c_1m_y(out_c, -1, plan(0, 2))
         else
            ! Default : transpose if needed
            if (dims(2) > 1) call transpose_z_to_y(wk13, wk2_r2c, sp)
            call c2c_1m_y(wk2_r2c, -1, plan(0, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) call transpose_y_to_x(wk2_r2c, out_c, sp)
         call c2c_1m_x(out_c, -1, plan(0, 1))

      end if

      !$acc end data

      if (decomp_profiler_fft) call decomp_profiler_end("fft_r2c")

   end subroutine fft_3d_r2c

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D inverse FFT - complex to real
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2r(in_c, out_r)

      implicit none

      ! Arguments
      complex(mytype), dimension(:, :, :), intent(INOUT) :: in_c
      real(mytype), dimension(:, :, :), intent(OUT) :: out_r

      ! Local variables
      complex(mytype), allocatable, dimension(:, :, :) :: wk1

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2r")

      !$acc data create(wk2_r2c,wk13) present(in_c,out_r)

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in Z =====
         if (inplace) then
            call c2c_1m_z(in_c, 1, plan(2, 3))
         else
            call alloc_z(wk1, sp)
            !$acc enter data create(wk1) async
            !$acc wait
            !$acc kernels default(present)
            wk1(:, :, :) = in_c(:, :, :)
            !$acc end kernels
            call c2c_1m_z(wk1, 1, plan(2, 3))
         end if

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(2)==1) then
            ! Slab : input available in Y and Z
            if (inplace) then
               call c2c_1m_y(in_c, 1, plan(2, 2))
            else
               call c2c_1m_y(wk1, 1, plan(2, 2))
            end if
         else if (dims(1)==1) then
            ! Slab : final buffer available in X and Y
            if (inplace) then
               call transpose_z_to_y(in_c, wk13, sp)
            else
               call transpose_z_to_y(wk1, wk13, sp)
            end if
            call c2c_1m_y(wk13, 1, plan(2, 2))
         else
            ! Pencil : use Y buffer
            if (inplace) then
               call transpose_z_to_y(in_c, wk2_r2c, sp)
            else
               call transpose_z_to_y(wk1, wk2_r2c, sp)
            end if
            call c2c_1m_y(wk2_r2c, 1, plan(2, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(2)==1 .and. dims(1)==1) then
            ! Single rank : input available in X, Y and Z
            if (inplace) then
               call c2r_1m_x(in_c, out_r)
            else
               call c2r_1m_x(wk1, out_r)
            end if
         else
            ! Default : transpose if needed, use final buffer
            if (dims(2)==1) then
               if (inplace) then
                  call transpose_y_to_x(in_c, wk13, sp)
               else
                  call transpose_y_to_x(wk1, wk13, sp)
               end if
            else if (dims(1) > 1) then
               call transpose_y_to_x(wk2_r2c, wk13, sp)
            end if
            call c2r_1m_x(wk13, out_r)
         end if

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in X =====
         if (inplace) then
            call c2c_1m_x(in_c, 1, plan(2, 1))
         else
            call alloc_x(wk1, sp)
            !$acc enter data create(wk1) async
            !$acc wait
            !$acc kernels default(present)
            wk1(:, :, :) = in_c(:, :, :)
            !$acc end kernels
            call c2c_1m_x(wk1, 1, plan(2, 1))
         end if

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1)==1) then
            ! Slab : input available in X and Y
            if (inplace) then
               call c2c_1m_y(in_c, 1, plan(2, 2))
            else
               call c2c_1m_y(wk1, 1, plan(2, 2))
            end if
         else
            if (dims(2)==1) then
               ! Slab : final buffer available in Y and Z
               if (dims(1) > 1) then
                  if (inplace) then
                     call transpose_x_to_y(in_c, wk13, sp)
                  else
                     call transpose_x_to_y(wk1, wk13, sp)
                  end if
               end if
               call c2c_1m_y(wk13, 1, plan(2, 2))
            else
               ! Pencil : use Y buffer
               if (inplace) then
                  call transpose_x_to_y(in_c, wk2_r2c, sp)
               else
                  call transpose_x_to_y(wk1, wk2_r2c, sp)
               end if
               call c2c_1m_y(wk2_r2c, 1, plan(2, 2))
            end if
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1)==1 .and. dims(2)==1) then
            ! Single rank : input available in X, Y and Z
            if (inplace) then
               call c2r_1m_z(in_c, out_r)
            else
               call c2r_1m_z(wk1, out_r)
            end if
         else
            ! Default : transpose if needed, use final buffer
            if (dims(1)==1) then
               if (inplace) then
                  call transpose_y_to_z(in_c, wk13, sp)
               else
                  call transpose_y_to_z(wk1, wk13, sp)
               end if
            else if (dims(2) > 1) then
               call transpose_y_to_z(wk2_r2c, wk13, sp)
            end if
            call c2r_1m_z(wk13, out_r)
         end if

      end if

      ! Free memory
      if (allocated(wk1)) then
         !$acc exit data delete(wk1) async
         !$acc wait
         deallocate (wk1)
      end if

      !$acc end data

      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2r")

   end subroutine fft_3d_c2r

end module decomp_2d_fft
