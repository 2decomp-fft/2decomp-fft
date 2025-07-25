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

         if (dims(1) > 1) then
            if (inplace) then
               call transpose_x_to_y(in, wk2_c2c, ph)
            else
               call transpose_x_to_y(wk1, wk2_c2c, ph)
            end if
            call c2c_1m_y(wk2_c2c, isign, plan(isign, 2))
         else
            if (inplace) then
               call c2c_1m_y(in, isign, plan(isign, 2))
            else
               call c2c_1m_y(wk1, isign, plan(isign, 2))
            end if
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_c2c, out, ph)
         else
            if (inplace) then
               call transpose_y_to_z(in, out, ph)
            else
               call transpose_y_to_z(wk1, out, ph)
            end if
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
         if (dims(1) > 1) then
            if (inplace) then
               call transpose_z_to_y(in, wk2_c2c, ph)
            else
               call transpose_z_to_y(wk1, wk2_c2c, ph)
            end if
            call c2c_1m_y(wk2_c2c, isign, plan(isign, 2))
         else  ! out==wk2_c2c if 1D decomposition
            if (inplace) then
               call transpose_z_to_y(in, out, ph)
            else
               call transpose_z_to_y(wk1, out, ph)
            end if
            call c2c_1m_y(out, isign, plan(isign, 2))
         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
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

      !use nvtx

      implicit none

      ! Arguments
      real(mytype), dimension(:, :, :), intent(IN) :: in_r
      complex(mytype), dimension(:, :, :), intent(OUT) :: out_c

      ! Local variables
#ifdef DEBUG
      integer :: i, j, k
      integer, dimension(3) :: dim3d
#endif

      if (decomp_profiler_fft) call decomp_profiler_start("fft_r2c")

      !$acc data create(wk13,wk2_r2c) present(in_r,out_c)

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in X =====
         call r2c_1m_x(in_r, wk13)

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            call transpose_x_to_y(wk13, wk2_r2c, sp)
            call c2c_1m_y(wk2_r2c, -1, plan(0, 2))
         else
            call c2c_1m_y(wk13, -1, plan(0, 2))
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_r2c, out_c, sp)
         else
            call transpose_y_to_z(wk13, out_c, sp)
         end if
         call c2c_1m_z(out_c, -1, plan(0, 3))

      else if (format == PHYSICAL_IN_Z) then

#ifdef DEBUG
         dim3d = shape(in_r)
         do k = 1, dim3d(3), dim3d(3) / 8
            do j = 1, dim3d(2), dim3d(2) / 8
               do i = 1, dim3d(1), dim3d(1) / 8
                  print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(in_r(i, j, k))
               end do
            end do
         end do
#endif

         ! ===== 1D FFTs in Z =====
         !call nvtxStartRange("Z r2c_1m_z")
         call r2c_1m_z(in_r, wk13)
         !call nvtxEndRange("Z r2c_1m_z")

         ! ===== Swap Z --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            !call nvtxStartRange("Z1 transpose_z_to_y")
            call transpose_z_to_y(wk13, wk2_r2c, sp)
            !call nvtxEndRange
            !call nvtxStartRange("Z1 c2c_1m_y")
            call c2c_1m_y(wk2_r2c, -1, plan(0, 2))
            !call nvtxEndRange

#ifdef DEBUG
            write (*, *) 'c2c_1m_y'
            dim3d = shape(wk2_r2c)
            do k = 1, dim3d(3), dim3d(3) / 8
               do j = 1, dim3d(2), dim3d(2) / 8
                  do i = 1, dim3d(1), dim3d(1) / 8
                     print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(wk2_r2c(i, j, k)), &
                        aimag(wk2_r2c(i, j, k))
                  end do
               end do
            end do
            write (*, *)
            write (*, *)
#endif

         else  ! out_c==wk2_r2c if 1D decomposition
            !call nvtxStartRange("Z transpose_z_to_y")
            call transpose_z_to_y(wk13, out_c, sp)
            !call nvtxEndRange
            !call nvtxStartRange("Z c2c_1m_y")
            call c2c_1m_y(out_c, -1, plan(0, 2))
            !call nvtxEndRange

#ifdef DEBUG
            write (*, *) 'c2c_1m_y2'
            dim3d = shape(out_c)
            do k = 1, dim3d(3), dim3d(3) / 8
               do j = 1, dim3d(2), dim3d(2) / 8
                  do i = 1, dim3d(1), dim3d(1) / 8
                     print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(out_c(i, j, k)), &
                        aimag(out_c(i, j, k))
                  end do
               end do
            end do
            write (*, *)
            write (*, *)
#endif

         end if

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            !call nvtxStartRange("Z1 transpose_y_to_x")
            call transpose_y_to_x(wk2_r2c, out_c, sp)
            !call nvtxEndRange
         end if
         !call nvtxStartRange("c2c_1m_x")
         call c2c_1m_x(out_c, -1, plan(0, 1))
         !call nvtxEndRange

#ifdef DEBUG
         write (*, *) 'c2c_1m_x'
         dim3d = shape(out_c)
         do k = 1, dim3d(3), dim3d(3) / 8
            do j = 1, dim3d(2), dim3d(2) / 8
               do i = 1, dim3d(1), dim3d(1) / 8
                  print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(out_c(i, j, k)), &
                     aimag(out_c(i, j, k))
               end do
            end do
         end do
         write (*, *)
         write (*, *)
#endif

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
#ifdef DEBUG
      integer :: i, j, k
      integer, dimension(3) :: dim3d
#endif

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
         if (inplace) then
            call transpose_z_to_y(in_c, wk2_r2c, sp)
         else
            call transpose_z_to_y(wk1, wk2_r2c, sp)
         end if
         call c2c_1m_y(wk2_r2c, 1, plan(2, 2))

         ! ===== Swap Y --> X; 1D FFTs in X =====
         if (dims(1) > 1) then
            call transpose_y_to_x(wk2_r2c, wk13, sp)
            call c2r_1m_x(wk13, out_r)
         else
            call c2r_1m_x(wk2_r2c, out_r)
         end if

      else if (format == PHYSICAL_IN_Z) then

#ifdef DEBUG
         write (*, *) 'Back Init c2c_1m_x line 788'
         dim3d = shape(in_c)
         do k = 1, dim3d(3), dim3d(3) / 8
            do j = 1, dim3d(2), dim3d(2) / 8
               do i = 1, dim3d(1), dim3d(1) / 8
                  print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(in_c(i, j, k)), &
                     aimag(in_c(i, j, k))
               end do
            end do
         end do
         write (*, *)
         write (*, *)
#endif

         ! ===== 1D FFTs in X =====
         if (inplace) then
            call c2c_1m_x(in_c, 1, plan(2, 1))

#ifdef DEBUG
            write (*, *) 'Back c2c_1m_x overwrite line 804'
            dim3d = shape(in_c)
            do k = 1, dim3d(3), dim3d(3) / 8
               do j = 1, dim3d(2), dim3d(2) / 8
                  do i = 1, dim3d(1), dim3d(1) / 8
                     print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(in_c(i, j, k)), &
                        aimag(in_c(i, j, k))
                  end do
               end do
            end do
            write (*, *)
            write (*, *)
#endif

         else
            call alloc_x(wk1, sp)
            !$acc enter data create(wk1) async
            !$acc wait
            !$acc kernels default(present)
            wk1(:, :, :) = in_c(:, :, :)
            !$acc end kernels
            call c2c_1m_x(wk1, 1, plan(2, 1))

#ifdef DEBUG
            write (*, *) 'Back2 c2c_1m_x line 821'
            dim3d = shape(wk1)
            do k = 1, dim3d(3), dim3d(1) / 8
               do j = 1, dim3d(2), dim3d(2) / 8
                  do i = 1, dim3d(1), dim3d(1) / 8
                     print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(wk1(i, j, k)), &
                        aimag(wk1(i, j, k))
                  end do
               end do
            end do
            write (*, *)
            write (*, *)
#endif

         end if

         ! ===== Swap X --> Y; 1D FFTs in Y =====
         if (dims(1) > 1) then
            if (inplace) then
               call transpose_x_to_y(in_c, wk2_r2c, sp)
            else
               call transpose_x_to_y(wk1, wk2_r2c, sp)
            end if
            call c2c_1m_y(wk2_r2c, 1, plan(2, 2))

#ifdef DEBUG
            write (*, *) 'Back c2c_1m_y line 844'
            dim3d = shape(wk2_r2c)
            do k = 1, dim3d(3), dim3d(3) / 8
               do j = 1, dim3d(2), dim3d(2) / 8
                  do i = 1, dim3d(1), dim3d(1) / 8
                     print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(wk2_r2c(i, j, k)), &
                        aimag(wk2_r2c(i, j, k))
                  end do
               end do
            end do
            write (*, *)
            write (*, *)
#endif

         else  ! in_c==wk2_r2c if 1D decomposition
            if (inplace) then
               call c2c_1m_y(in_c, 1, plan(2, 2))

#ifdef DEBUG
               write (*, *) 'Back2 c2c_1m_y line 860'
               dim3d = shape(in_c)
               do k = 1, dim3d(3), dim3d(3) / 8
                  do j = 1, dim3d(2), dim3d(2) / 8
                     do i = 1, dim3d(1), dim3d(1) / 8
                        print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(in_c(i, j, k)), &
                           aimag(in_c(i, j, k))
                     end do
                  end do
               end do
               write (*, *)
               write (*, *)
#endif

            else
               call c2c_1m_y(wk1, 1, plan(2, 2))

#ifdef DEBUG
               write (*, *) 'Back3 c2c_1m_y line 875'
               dim3d = shape(wk1)
               do k = 1, dim3d(3), dim3d(3) / 8
                  do j = 1, dim3d(2), dim3d(2) / 8
                     do i = 1, dim3d(1), dim3d(1) / 8
                        print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(wk1(i, j, k)), &
                           aimag(wk1(i, j, k))
                     end do
                  end do
               end do
               write (*, *)
               write (*, *)
#endif

            end if
         end if

         ! ===== Swap Y --> Z; 1D FFTs in Z =====
         if (dims(1) > 1) then
            call transpose_y_to_z(wk2_r2c, wk13, sp)
         else
            if (inplace) then
               call transpose_y_to_z(in_c, wk13, sp)
            else
               call transpose_y_to_z(wk1, wk13, sp)
            end if
         end if

#ifdef DEBUG
         write (*, *) 'Back2 after tr_y2z'
         dim3d = shape(wk13)
         do k = 1, dim3d(3), dim3d(3) / 8
            do j = 1, dim3d(2), dim3d(2) / 8
               do i = 1, dim3d(1), dim3d(1) / 8
                  print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(wk13(i, j, k)), &
                     aimag(wk13(i, j, k))
               end do
            end do
         end do
         write (*, *)
         write (*, *)
#endif

         call c2r_1m_z(wk13, out_r)

#ifdef DEBUG
         write (*, *) 'Back2 c2r_1m_z out_r line 902'
         dim3d = shape(out_r)
         do k = 1, dim3d(3), dim3d(3) / 8
            do j = 1, dim3d(2), dim3d(2) / 8
               do i = 1, dim3d(1), dim3d(1) / 8
                  print "(i3,1x,i3,1x,i3,1x,e12.5,1x,e12.5)", i, j, k, real(out_r(i, j, k))
               end do
            end do
         end do
         write (*, *)
         write (*, *)
#endif

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
