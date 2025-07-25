!! SPDX-License-Identifier: BSD-3-Clause

! This is the Intel MKL implementation of the FFT library

module decomp_2d_fft

   use iso_c_binding, only: c_f_pointer, c_loc
   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use decomp_2d_profiler
   use m_decomp_pool
   use MKL_DFTI   ! MKL FFT module

   implicit none

   private        ! Make everything private unless declared public

   ! engine-specific global variables

   ! Descriptors for MKL FFT, one for each set of 1D FFTs
   type(DFTI_DESCRIPTOR), pointer, save :: c2c_x => null(), & ! c2c transforms
                                           c2c_y => null(), &
                                           c2c_z => null(), &
                                           r2c_x => null(), & ! r2c/c2r, physical in x
                                           c2c_y2 => null(), &
                                           c2c_z2 => null(), &
                                           c2r_x => null(), &
                                           r2c_z => null(), & ! r2c/c2r, physical in z
                                           c2c_x2 => null(), &
                                           c2r_z => null()

   integer, parameter, public :: D2D_FFT_BACKEND = D2D_FFT_BACKEND_MKL

   ! Derived type with all the quantities needed to perform FFT
   type decomp_2d_fft_engine
      ! Engine-specific stuff
      type(DFTI_DESCRIPTOR), private, pointer :: c2c_x => null(), & ! c2c transforms
                                                 c2c_y => null(), &
                                                 c2c_z => null(), &
                                                 r2c_x => null(), & ! r2c/c2r, physical in x
                                                 c2c_y2 => null(), &
                                                 c2c_z2 => null(), &
                                                 c2r_x => null(), &
                                                 r2c_z => null(), & ! r2c/c2r, physical in z
                                                 c2c_x2 => null(), &
                                                 c2r_z => null()
      ! All the engines have this
      integer, private :: format
      logical, private :: initialised = .false.
      integer, private :: nx_fft, ny_fft, nz_fft
      type(decomp_info), pointer, public :: ph => null()
      type(decomp_info), private :: ph_target ! ph => ph_target or ph => decomp_main
      type(decomp_info), public :: sp
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

   ! common code used for all engines, including global variables,
   ! generic interface definitions and several subroutines
#include "fft_common.f90"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time initialisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine init_fft_engine(engine)

      implicit none

      type(decomp_2d_fft_engine), target, intent(inout) :: engine

      call decomp_2d_fft_log("MKL")

      ! For C2C transforms
      call c2c_1m_x_plan(engine%c2c_x, ph)
      call c2c_1m_y_plan(engine%c2c_y, ph)
      call c2c_1m_z_plan(engine%c2c_z, ph)

      ! For R2C/C2R tranfroms with physical space in X-pencil
      if (format == PHYSICAL_IN_X) then
         call r2c_1m_x_plan(engine%r2c_x, ph, sp, -1)
         call c2c_1m_y_plan(engine%c2c_y2, sp)
         call c2c_1m_z_plan(engine%c2c_z2, sp)
         call r2c_1m_x_plan(engine%c2r_x, ph, sp, 1)

         ! For R2C/C2R tranfroms with physical space in Z-pencil
      else if (format == PHYSICAL_IN_Z) then
         call r2c_1m_z_plan(engine%r2c_z, ph, sp, -1)
         call c2c_1m_y_plan(engine%c2c_y2, sp)
         call c2c_1m_x_plan(engine%c2c_x2, sp)
         call r2c_1m_z_plan(engine%c2r_z, ph, sp, 1)
      end if

      return
   end subroutine init_fft_engine

   ! Return an MKL plan for multiple 1D c2c FFTs in X direction
   subroutine c2c_1m_x_plan(desc, decomp)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: status

#ifdef DOUBLE_PREC
      status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
                                    DFTI_COMPLEX, 1, decomp%xsz(1))
#else
      status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
                                    DFTI_COMPLEX, 1, decomp%xsz(1))
#endif
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCreateDescriptor")
      status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
                            decomp%xsz(2) * decomp%xsz(3))
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      if (inplace) then
         status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_INPLACE)
      else
         status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      end if
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, decomp%xsz(1))
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, decomp%xsz(1))
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiCommitDescriptor(desc)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCommitDescriptor")

      return
   end subroutine c2c_1m_x_plan

   ! Return an MKL plan for multiple 1D c2c FFTs in Y direction
   subroutine c2c_1m_y_plan(desc, decomp)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: status, strides(2)

      ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
      ! done one Z-plane at a time. So plan for 2D data sets here.

#ifdef DOUBLE_PREC
      status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
                                    DFTI_COMPLEX, 1, decomp%ysz(2))
#else
      status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
                                    DFTI_COMPLEX, 1, decomp%ysz(2))
#endif
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCreateDescriptor")
      if (inplace) then
         status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_INPLACE)
      else
         status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      end if
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, decomp%ysz(1))
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, 1)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, 1)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      strides(1) = 0
      strides(2) = decomp%ysz(1)
      status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiCommitDescriptor(desc)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCommitDescriptor")

      return
   end subroutine c2c_1m_y_plan

   ! Return an MKL plan for multiple 1D c2c FFTs in Z direction
   subroutine c2c_1m_z_plan(desc, decomp)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      TYPE(DECOMP_INFO), intent(IN) :: decomp

      integer :: status, strides(2)

#ifdef DOUBLE_PREC
      status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
                                    DFTI_COMPLEX, 1, decomp%zsz(3))
#else
      status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
                                    DFTI_COMPLEX, 1, decomp%zsz(3))
#endif
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCreateDescriptor")
      if (inplace) then
         status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_INPLACE)
      else
         status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      end if
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
                            decomp%zsz(1) * decomp%zsz(2))
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, 1)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, 1)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      strides(1) = 0
      strides(2) = decomp%zsz(1) * decomp%zsz(2)
      status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiCommitDescriptor(desc)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCommitDescriptor")

      return
   end subroutine c2c_1m_z_plan

   ! Return an MKL plan for multiple 1D r2c FFTs in X direction
   subroutine r2c_1m_x_plan(desc, decomp_ph, decomp_sp, direction)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph, decomp_sp
      integer, intent(IN) :: direction ! (-1=r2c; 1=c2r)

      integer :: status

      ! c2r and r2c plans are almost the same, just swap input/output

#ifdef DOUBLE_PREC
      status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
                                    DFTI_REAL, 1, decomp_ph%xsz(1))
#else
      status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
                                    DFTI_REAL, 1, decomp_ph%xsz(1))
#endif
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCreateDescriptor")
      status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
                            decomp_ph%xsz(2) * decomp_ph%xsz(3))
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
                            DFTI_COMPLEX_COMPLEX)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      if (direction == -1) then  ! r2c
         status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, &
                               decomp_ph%xsz(1))
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
         status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, &
                               decomp_sp%xsz(1))
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      else if (direction == 1) then  ! c2r
         status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, &
                               decomp_sp%xsz(1))
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
         status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, &
                               decomp_ph%xsz(1))
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      end if
      status = DftiCommitDescriptor(desc)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCommitDescriptor")

      return
   end subroutine r2c_1m_x_plan

   ! Return an MKL plan for multiple 1D r2c FFTs in Z direction
   subroutine r2c_1m_z_plan(desc, decomp_ph, decomp_sp, direction)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      TYPE(DECOMP_INFO), intent(IN) :: decomp_ph, decomp_sp
      integer, intent(IN) :: direction ! (-1=r2c; 1=c2r)

      integer :: status, strides(2)

      ! c2r and r2c plans are almost the same, just swap input/output

#ifdef DOUBLE_PREC
      status = DftiCreateDescriptor(desc, DFTI_DOUBLE, &
                                    DFTI_REAL, 1, decomp_ph%zsz(3))
#else
      status = DftiCreateDescriptor(desc, DFTI_SINGLE, &
                                    DFTI_REAL, 1, decomp_ph%zsz(3))
#endif
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCreateDescriptor")
      status = DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, &
                            decomp_ph%zsz(1) * decomp_ph%zsz(2))
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_CONJUGATE_EVEN_STORAGE, &
                            DFTI_COMPLEX_COMPLEX)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_INPUT_DISTANCE, 1)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, 1)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      strides(1) = 0
      strides(2) = decomp_ph%zsz(1) * decomp_ph%zsz(2)
      if (direction == -1) then
         status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
      else if (direction == 1) then
         status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
      end if
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      strides(2) = decomp_sp%zsz(1) * decomp_sp%zsz(2)
      if (direction == -1) then
         status = DftiSetValue(desc, DFTI_OUTPUT_STRIDES, strides)
      else if (direction == 1) then
         status = DftiSetValue(desc, DFTI_INPUT_STRIDES, strides)
      end if
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiSetValue")
      status = DftiCommitDescriptor(desc)
      if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiCommitDescriptor")

      return
   end subroutine r2c_1m_z_plan

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This routine performs one-time finalisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine finalize_fft_engine(engine)

      implicit none

      type(decomp_2d_fft_engine), optional :: engine

      integer :: status

      if (present(engine)) then

         status = DftiFreeDescriptor(engine%c2c_x)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
         status = DftiFreeDescriptor(engine%c2c_y)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
         status = DftiFreeDescriptor(engine%c2c_z)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
         if (engine%format == PHYSICAL_IN_X) then
            status = DftiFreeDescriptor(engine%r2c_x)
            if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
            status = DftiFreeDescriptor(engine%c2c_z2)
            if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
            status = DftiFreeDescriptor(engine%c2r_x)
            if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
         else if (engine%format == PHYSICAL_IN_Z) then
            status = DftiFreeDescriptor(engine%r2c_z)
            if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
            status = DftiFreeDescriptor(engine%c2c_x2)
            if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
            status = DftiFreeDescriptor(engine%c2r_z)
            if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")
         end if
         status = DftiFreeDescriptor(engine%c2c_y2)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "DftiFreeDescriptor")

      else

         if (associated(c2c_x)) nullify (c2c_x)
         if (associated(c2c_y)) nullify (c2c_y)
         if (associated(c2c_z)) nullify (c2c_z)
         if (associated(r2c_x)) nullify (r2c_x)
         if (associated(c2c_z2)) nullify (c2c_z2)
         if (associated(c2r_x)) nullify (c2r_x)
         if (associated(r2c_z)) nullify (r2c_z)
         if (associated(c2c_x2)) nullify (c2c_x2)
         if (associated(c2r_z)) nullify (c2r_z)
         if (associated(c2c_y2)) nullify (c2c_y2)

      end if

   end subroutine finalize_fft_engine

   ! Use engine-specific stuff
   subroutine use_fft_engine(engine)

      implicit none

      type(decomp_2d_fft_engine), target, intent(in) :: engine

      c2c_x => engine%c2c_x
      c2c_y => engine%c2c_y
      c2c_z => engine%c2c_z
      if (format == PHYSICAL_IN_X) then
         r2c_x => engine%r2c_x
         c2c_z2 => engine%c2c_z2
         c2r_x => engine%c2r_x
      else
         r2c_z => engine%r2c_z
         c2c_x2 => engine%c2c_x2
         c2r_z => engine%c2r_z
      end if
      c2c_y2 => engine%c2c_y2

   end subroutine use_fft_engine

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D FFT - complex to complex
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2c(in, out, isign)

      implicit none

      ! Arguments
      complex(mytype), dimension(:, :, :), intent(IN) :: in
      complex(mytype), dimension(:, :, :), intent(OUT) :: out
      integer, intent(IN) :: isign

      ! Local variables
      complex(mytype), contiguous, pointer, dimension(:, :, :) :: wk2_c2c, wk1, wk2b, wk3
      integer :: k, status

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2c")

      ! Init pointers
      call decomp_pool_get(wk2_c2c, ph%ysz)
      nullify (wk1)
      nullify (wk2b)
      nullify (wk3)

      if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_FORWARD .OR. &
          format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_BACKWARD) then

         ! ===== 1D FFTs in X =====
         !       if (isign==DECOMP_2D_FFT_FORWARD) then
         !          status = DftiComputeForward(c2c_x, in(:,1,1), wk1(:,1,1))
         !       else if (isign==DECOMP_2D_FFT_BACKWARD) then
         !          status = DftiComputeBackward(c2c_x, in(:,1,1), wk1(:,1,1))
         !       end if
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_x, in, isign)
         else
            call decomp_pool_get(wk1, ph%xsz)
            status = wrapper_c2c(c2c_x, in, wk1, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

         ! ===== Swap X --> Y =====
         if (inplace) then
            call transpose_x_to_y(in, wk2_c2c, ph)
         else
            call transpose_x_to_y(wk1, wk2_c2c, ph)
         end if

         ! ===== 1D FFTs in Y =====
         if (inplace) then
            do k = 1, ph%ysz(3)
               status = wrapper_c2c_inplace(c2c_y, wk2_c2c(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         else
            call decomp_pool_get(wk2b, ph%ysz)
            do k = 1, ph%ysz(3) ! one Z-plane at a time
               !          if (isign==DECOMP_2D_FFT_FORWARD) then
               !             status = DftiComputeForward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
               !          else if (isign==DECOMP_2D_FFT_BACKWARD) then
               !             status = DftiComputeBackward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
               !          end if
               status = wrapper_c2c(c2c_y, wk2_c2c(:, :, k), wk2b(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         end if

         ! ===== Swap Y --> Z =====
         if (inplace) then
            call transpose_y_to_z(wk2_c2c, out, ph)
         else
            call decomp_pool_get(wk3, ph%zsz)
            call transpose_y_to_z(wk2b, wk3, ph)
         end if

         ! ===== 1D FFTs in Z =====
         !       if (isign==DECOMP_2D_FFT_FORWARD) then
         !          status = DftiComputeForward(c2c_z, wk3(:,1,1), out(:,1,1))
         !       else if (isign==DECOMP_2D_FFT_BACKWARD) then
         !          status = DftiComputeBackward(c2c_z, wk3(:,1,1), out(:,1,1))
         !       end if
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_z, out, isign)
         else
            status = wrapper_c2c(c2c_z, wk3, out, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

      else if (format == PHYSICAL_IN_X .AND. isign == DECOMP_2D_FFT_BACKWARD &
               .OR. &
               format == PHYSICAL_IN_Z .AND. isign == DECOMP_2D_FFT_FORWARD) then

         ! ===== 1D FFTs in Z =====
         !       if (isign==DECOMP_2D_FFT_FORWARD) then
         !          status = DftiComputeForward(c2c_z, in(:,1,1), wk1(:,1,1))
         !       else if (isign==DECOMP_2D_FFT_BACKWARD) then
         !          status = DftiComputeBackward(c2c_z, in(:,1,1), wk1(:,1,1))
         !       end if
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_z, in, isign)
         else
            call decomp_pool_get(wk1, ph%zsz)
            status = wrapper_c2c(c2c_z, in, wk1, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

         ! ===== Swap Z --> Y =====
         if (inplace) then
            call transpose_z_to_y(in, wk2_c2c, ph)
         else
            call transpose_z_to_y(wk1, wk2_c2c, ph)
         end if

         ! ===== 1D FFTs in Y =====
         if (inplace) then
            do k = 1, ph%ysz(3)
               status = wrapper_c2c_inplace(c2c_y, wk2_c2c(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         else
            call decomp_pool_get(wk2b, ph%ysz)
            do k = 1, ph%ysz(3) ! one Z-plane at a time
               !          if (isign==DECOMP_2D_FFT_FORWARD) then
               !             status = DftiComputeForward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
               !          else if (isign==DECOMP_2D_FFT_BACKWARD) then
               !             status = DftiComputeBackward(c2c_y, wk2(:,1,k), wk2b(:,1,k))
               !          end if
               status = wrapper_c2c(c2c_y, wk2_c2c(:, :, k), wk2b(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         end if

         ! ===== Swap Y --> X =====
         if (inplace) then
            call transpose_y_to_x(wk2_c2c, out, ph)
         else
            call decomp_pool_get(wk3, ph%xsz)
            call transpose_y_to_x(wk2b, wk3, ph)
         end if

         ! ===== 1D FFTs in X =====
         !       if (isign==DECOMP_2D_FFT_FORWARD) then
         !          status = DftiComputeForward(c2c_x, wk3(:,1,1), out(:,1,1))
         !       else if (isign==DECOMP_2D_FFT_BACKWARD) then
         !          status = DftiComputeBackward(c2c_x, wk3(:,1,1), out(:,1,1))
         !       end if
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_x, out, isign)
         else
            status = wrapper_c2c(c2c_x, wk3, out, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

      end if

      ! Free memory
      call decomp_pool_free(wk2_c2c)
      if (associated(wk1)) call decomp_pool_free(wk1)
      if (associated(wk2b)) call decomp_pool_free(wk2b)
      if (associated(wk3)) call decomp_pool_free(wk3)

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

      ! Local variables
      complex(mytype), pointer, contiguous, dimension(:, :, :) :: wk2_r2c, wk13, wk2b, wk3
      integer :: k, status, isign

      if (decomp_profiler_fft) call decomp_profiler_start("fft_r2c")

      ! Init pointers to null
      call decomp_pool_get(wk2_r2c, sp%ysz)
      if (format == PHYSICAL_IN_X) then
         call decomp_pool_get(wk13, sp%xsz)
      else
         call decomp_pool_get(wk13, sp%zsz)
      end if
      nullify (wk2b)
      nullify (wk3)

      isign = DECOMP_2D_FFT_FORWARD

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in X =====
         !       status = DftiComputeForward(r2c_x, in_r(:,1,1), wk1(:,1,1))
         status = wrapper_r2c(r2c_x, in_r, wk13)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_r2c")

         ! ===== Swap X --> Y =====
         call transpose_x_to_y(wk13, wk2_r2c, sp)

         ! ===== 1D FFTs in Y =====
         if (inplace) then
            do k = 1, sp%ysz(3)
               status = wrapper_c2c_inplace(c2c_y2, wk2_r2c(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         else
            call decomp_pool_get(wk2b, sp%ysz)
            do k = 1, sp%ysz(3)
               !          status = DftiComputeForward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
               status = wrapper_c2c(c2c_y2, wk2_r2c(:, :, k), wk2b(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         end if

         ! ===== Swap Y --> Z =====
         if (inplace) then
            call transpose_y_to_z(wk2_r2c, out_c, sp)
         else
            call decomp_pool_get(wk3, sp%zsz)
            call transpose_y_to_z(wk2b, wk3, sp)
         end if

         ! ===== 1D FFTs in Z =====
         !       status = DftiComputeForward(c2c_z2, wk3(:,1,1), out_c(:,1,1))
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_z2, out_c, isign)
         else
            status = wrapper_c2c(c2c_z2, wk3, out_c, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in Z =====
         !       status = DftiComputeForward(r2c_z, in_r(:,1,1), wk1(:,1,1))
         status = wrapper_r2c(r2c_z, in_r, wk13)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_r2c")

         ! ===== Swap Z --> Y =====
         call transpose_z_to_y(wk13, wk2_r2c, sp)

         ! ===== 1D FFTs in Y =====
         if (inplace) then
            do k = 1, sp%ysz(3)
               status = wrapper_c2c_inplace(c2c_y2, wk2_r2c(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         else
            call decomp_pool_get(wk2b, sp%ysz)
            do k = 1, sp%ysz(3)
               !          status = DftiComputeForward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
               status = wrapper_c2c(c2c_y2, wk2_r2c(:, :, k), wk2b(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         end if

         ! ===== Swap Y --> X =====
         if (inplace) then
            call transpose_y_to_x(wk2_r2c, out_c, sp)
         else
            call decomp_pool_get(wk3, sp%xsz)
            call transpose_y_to_x(wk2b, wk3, sp)
         end if

         ! ===== 1D FFTs in X =====
         !       status = DftiComputeForward(c2c_x2, wk3(:,1,1), out_c(:,1,1))
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_x2, out_c, isign)
         else
            status = wrapper_c2c(c2c_x2, wk3, out_c, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

      end if

      ! Free memory
      call decomp_pool_free(wk2_r2c)
      call decomp_pool_free(wk13)
      if (associated(wk2b)) call decomp_pool_free(wk2b)
      if (associated(wk3)) call decomp_pool_free(wk3)

      if (decomp_profiler_fft) call decomp_profiler_end("fft_r2c")

      return
   end subroutine fft_3d_r2c

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine fft_3d_c2r(in_c, out_r)

      implicit none

      ! Arguments
      complex(mytype), dimension(:, :, :), intent(IN) :: in_c
      real(mytype), dimension(:, :, :), intent(OUT) :: out_r

      ! Local variables
      complex(mytype), pointer, contiguous, dimension(:, :, :) :: wk2_r2c, wk13, wk1, wk2b
      integer :: k, status, isign

      if (decomp_profiler_fft) call decomp_profiler_start("fft_c2r")

      ! Init pointers to null
      call decomp_pool_get(wk2_r2c, sp%ysz)
      if (format == PHYSICAL_IN_X) then
         call decomp_pool_get(wk13, sp%xsz)
      else
         call decomp_pool_get(wk13, sp%zsz)
      end if
      nullify (wk1)
      nullify (wk2b)

      isign = DECOMP_2D_FFT_BACKWARD

      if (format == PHYSICAL_IN_X) then

         ! ===== 1D FFTs in Z =====
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_z2, in_c, isign)
         else
            call decomp_pool_get(wk1, sp%zsz)
            !       status = DftiComputeBackward(c2c_z2, in_c(:,1,1), wk1(:,1,1))
            status = wrapper_c2c(c2c_z2, in_c, wk1, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

         ! ===== Swap Z --> Y =====
         if (inplace) then
            call transpose_z_to_y(in_c, wk2_r2c, sp)
         else
            call transpose_z_to_y(wk1, wk2_r2c, sp)
         end if

         ! ===== 1D FFTs in Y =====
         if (inplace) then
            do k = 1, sp%ysz(3)
               status = wrapper_c2c_inplace(c2c_y2, wk2_r2c(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         else
            call decomp_pool_get(wk2b, sp%ysz)
            do k = 1, sp%ysz(3)
               !          status = DftiComputeBackward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
               status = wrapper_c2c(c2c_y2, wk2_r2c(:, :, k), wk2b(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         end if

         ! ===== Swap Y --> X =====
         if (inplace) then
            call transpose_y_to_x(wk2_r2c, wk13, sp)
         else
            call transpose_y_to_x(wk2b, wk13, sp)
         end if

         ! ===== 1D FFTs in X =====
         !       status = DftiComputeBackward(c2r_x, wk3(:,1,1), out_r(:,1,1))
         status = wrapper_c2r(c2r_x, wk13, out_r)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2r")

      else if (format == PHYSICAL_IN_Z) then

         ! ===== 1D FFTs in X =====
         if (inplace) then
            status = wrapper_c2c_inplace(c2c_x2, in_c, isign)
         else
            call decomp_pool_get(wk1, sp%xsz)
            !       status = DftiComputeBackward(c2c_x2, in_c(:,1,1), wk1(:,1,1))
            status = wrapper_c2c(c2c_x2, in_c, wk1, isign)
         end if
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")

         ! ===== Swap X --> Y =====
         if (inplace) then
            call transpose_x_to_y(in_c, wk2_r2c, sp)
         else
            call transpose_x_to_y(wk1, wk2_r2c, sp)
         end if

         ! ===== 1D FFTs in Y =====
         if (inplace) then
            do k = 1, sp%ysz(3)
               status = wrapper_c2c_inplace(c2c_y2, wk2_r2c(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         else
            call decomp_pool_get(wk2b, sp%ysz)
            do k = 1, sp%ysz(3)
               !          status = DftiComputeBackward(c2c_y2, wk2(:,1,k), wk2b(:,1,k))
               status = wrapper_c2c(c2c_y2, wk2_r2c(:, :, k), wk2b(:, :, k), isign)
               if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2c")
            end do
         end if

         ! ===== Swap Y --> Z =====
         if (inplace) then
            call transpose_y_to_z(wk2_r2c, wk13, sp)
         else
            call transpose_y_to_z(wk2b, wk13, sp)
         end if

         ! ===== 1D FFTs in Z =====
         !       status = DftiComputeBackward(c2r_z, wk3(:,1,1), out_r(:,1,1))
         status = wrapper_c2r(c2r_z, wk13, out_r)
         if (status /= 0) call decomp_2d_abort(__FILE__, __LINE__, status, "wrapper_c2r")

      end if

      ! Free memory
      call decomp_pool_free(wk2_r2c)
      call decomp_pool_free(wk13)
      if (associated(wk1)) call decomp_pool_free(wk1)
      if (associated(wk2b)) call decomp_pool_free(wk2b)

      if (decomp_profiler_fft) call decomp_profiler_end("fft_c2r")

   end subroutine fft_3d_c2r

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Wrapper functions so that one can pass 3D arrays to DftiCompute
   !  -- MKL accepts only 1D arrays as input/output for its multi-
   !     dimensional FFTs.
   !  -- Using EQUIVALENCE as suggested by MKL documents is impossible
   !     for allocated arrays, not to mention bad coding style
   !  -- All code commented out above may well work but not safe. There
   !     is no guarantee that compiler wouldn't make copies of 1D arrays
   !     (which would contain only one slice of the original 3D data)
   !     rather than referring to the same memory address, i.e. 3D array
   !     A and 1D array A(:,1,1) may refer to different memory location.
   !  -- Using the following wrappers is safe and standard conforming.
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   integer function wrapper_c2c(desc, in, out, isign)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      complex(mytype), dimension(*) :: in, out
      integer :: isign, status

      if ((associated(desc, c2c_x) .and. skip_x_c2c) .or. &
          (associated(desc, c2c_x2) .and. skip_x_c2c) .or. &
          (associated(desc, c2c_y) .and. skip_y_c2c) .or. &
          (associated(desc, c2c_y2) .and. skip_y_c2c) .or. &
          (associated(desc, c2c_z) .and. skip_z_c2c) .or. &
          (associated(desc, c2c_z2) .and. skip_z_c2c)) then
         if (associated(desc, c2c_x)) out(1:product(ph%xsz)) = in(1:product(ph%xsz))
         if (associated(desc, c2c_y)) out(1:ph%ysz(1) * ph%ysz(2)) = in(1:ph%ysz(1) * ph%ysz(2))
         if (associated(desc, c2c_z)) out(1:product(ph%zsz)) = in(1:product(ph%zsz))
         if (associated(desc, c2c_x2)) out(1:product(sp%xsz)) = in(1:product(sp%xsz))
         if (associated(desc, c2c_y2)) out(1:sp%ysz(1) * sp%ysz(2)) = in(1:sp%ysz(1) * sp%ysz(2))
         if (associated(desc, c2c_z2)) out(1:product(sp%zsz)) = in(1:product(sp%zsz))
         wrapper_c2c = 0
         return
      end if

      if (isign == DECOMP_2D_FFT_FORWARD) then
         status = DftiComputeForward(desc, in, out)
      else if (isign == DECOMP_2D_FFT_BACKWARD) then
         status = DftiComputeBackward(desc, in, out)
      end if

      wrapper_c2c = status

   end function wrapper_c2c

   integer function wrapper_c2c_inplace(desc, inout, isign)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      complex(mytype), dimension(*) :: inout
      integer :: isign, status

      if ((associated(desc, c2c_x) .and. skip_x_c2c) .or. &
          (associated(desc, c2c_x2) .and. skip_x_c2c) .or. &
          (associated(desc, c2c_y) .and. skip_y_c2c) .or. &
          (associated(desc, c2c_y2) .and. skip_y_c2c) .or. &
          (associated(desc, c2c_z) .and. skip_z_c2c) .or. &
          (associated(desc, c2c_z2) .and. skip_z_c2c)) then
         wrapper_c2c_inplace = 0
         return
      end if

      if (isign == DECOMP_2D_FFT_FORWARD) then
         status = DftiComputeForward(desc, inout)
      else if (isign == DECOMP_2D_FFT_BACKWARD) then
         status = DftiComputeBackward(desc, inout)
      end if

      wrapper_c2c_inplace = status

   end function wrapper_c2c_inplace

   integer function wrapper_r2c(desc, in, out)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      real(mytype), dimension(*) :: in
      complex(mytype), dimension(*) :: out

      if ((associated(desc, r2c_x) .and. skip_x_c2c) .or. &
          (associated(desc, r2c_z) .and. skip_z_c2c)) &
         call decomp_2d_warning(__FILE__, __LINE__, 1, &
                                "r2c / c2r transform can not be skipped")

      wrapper_r2c = DftiComputeForward(desc, in, out)

   end function wrapper_r2c

   integer function wrapper_c2r(desc, in, out)

      implicit none

      type(DFTI_DESCRIPTOR), pointer :: desc
      complex(mytype), dimension(*) :: in
      real(mytype), dimension(*) :: out

      if ((associated(desc, c2r_x) .and. skip_x_c2c) .or. &
          (associated(desc, c2r_z) .and. skip_z_c2c)) &
         call decomp_2d_warning(__FILE__, __LINE__, 2, &
                                "r2c / c2r transform can not be skipped")

      wrapper_c2r = DftiComputeBackward(desc, in, out)

   end function wrapper_c2r

end module decomp_2d_fft
