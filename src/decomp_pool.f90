!! SPDX-License-Identifier: BSD-3-Clause

!
! Provide a memory pool without exposing the underlying object
!
! Simplified interface for the external code
!
! The advanced interface is available when using directly the variable decomp_pool
!
module m_decomp_pool

   use iso_fortran_env, only : real32
   use iso_c_binding, only : c_size_t, c_loc, c_f_pointer
   use decomp_2d_constants, only : mytype, real_type
   use decomp_2d_mpi, only : decomp_2d_abort
   use m_info
   use m_mem_pool

   implicit none

   ! Main memory pool for the library 2decomp
   type(mem_pool), target, save :: decomp_pool

   ! true when the memory pool is ready
   logical, protected, save :: decomp_pool_ready = .false.

   ! Default type in the memory pool module
   integer, save :: decomp_pool_default_type = real_type

   ! Default : private
   private
   public :: decomp_pool, &
             decomp_pool_init, &
             decomp_pool_fin, &
             decomp_pool_get, &
             decomp_pool_free, &
             decomp_pool_ready, &
             decomp_pool_default_type, &
             decomp_pool_set_default_type

   interface decomp_pool_get
      module procedure decomp_pool_get_real
      module procedure decomp_pool_get_real1D
      module procedure decomp_pool_get_cplx
      module procedure decomp_pool_get_cplx1D
   end interface decomp_pool_get

   interface decomp_pool_free
      module procedure decomp_pool_free_real                                                 
      module procedure decomp_pool_free_real1D
      module procedure decomp_pool_free_cplx                                                 
      module procedure decomp_pool_free_cplx1D
   end interface decomp_pool_free

contains

   !
   ! Initialize the decomp pool
   !
   !   - decomp provides the size of the 3D array
   !   - type provides the type of the 3D arrays
   !   - blk_n provides the number of blocks to allocate
   !
   subroutine decomp_pool_init(decomp, type, blk_n)

      implicit none

      ! Arguments
      class(info), intent(in) :: decomp
      integer, intent(in), optional :: type
      integer, intent(in), optional :: blk_n

      if (present(type)) then
         call decomp_pool%init(shape = (/type, decomp%xsz/), blk_n = blk_n)
         call decomp_pool%new_shape(type, decomp)
      else
         call decomp_pool%init(shape = (/decomp_pool_default_type, decomp%xsz/), blk_n = blk_n)
         call decomp_pool%new_shape(decomp_pool_default_type, decomp)
      end if

      decomp_pool_ready = .true.

   end subroutine decomp_pool_init

   !
   ! Finalize the decomp pool
   !
   subroutine decomp_pool_fin()

      implicit none

      call decomp_pool%fin()

      decomp_pool_ready = .false.

   end subroutine decomp_pool_fin

   !
   ! Add a free memory block
   !
   subroutine decomp_pool_new()

      implicit none

      call decomp_pool%new(.false.)

   end subroutine decomp_pool_new

   !
   ! Get a block from the free memory list
   !
   ! If no block is available, a new one is created
   !
   subroutine decomp_pool_get_real(ptr, shape)

      implicit none

      ! Arguments
      real(mytype), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in), optional :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      call decomp_pool%get(ptr, shape)

   end subroutine decomp_pool_get_real

   subroutine decomp_pool_get_real1D(ptr)

      implicit none

      ! Arguments
      real(mytype), intent(out), dimension(:), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Local variable
      real(mytype), dimension(:, :, :), contiguous, pointer :: ptr3D
#ifdef _GPU
      attributes(device) :: ptr3D
#endif

      call decomp_pool%get(ptr3D)
      if (mytype == KIND(0._real32)) then
         call c_f_pointer(c_loc(ptr3D), ptr, (/decomp_pool%get_size()/))
      else
         call c_f_pointer(c_loc(ptr3D), ptr, (/decomp_pool%get_size()/2_c_size_t/))
      end if
      nullify(ptr3D)

   end subroutine decomp_pool_get_real1D

   subroutine decomp_pool_get_cplx(ptr, shape)

      implicit none

      ! Arguments
      complex(mytype), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in), optional :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      call decomp_pool%get(ptr, shape)

   end subroutine decomp_pool_get_cplx

   subroutine decomp_pool_get_cplx1D(ptr)                                                  
      
      implicit none                                                                        
      
      ! Arguments
      complex(mytype), intent(out), dimension(:), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Local variable
      complex(mytype), dimension(:, :, :), contiguous, pointer :: ptr3D
#ifdef _GPU
      attributes(device) :: ptr3D
#endif
      
      call decomp_pool%get(ptr3D)
      if (mytype == KIND(0._real32)) then
         call c_f_pointer(c_loc(ptr3D), ptr, (/decomp_pool%get_size()/2_c_size_t/))                         
      else
         call c_f_pointer(c_loc(ptr3D), ptr, (/decomp_pool%get_size()/4_c_size_t/))              
      end if
      nullify(ptr3D)                                                                       
   
   end subroutine decomp_pool_get_cplx1D

   !
   ! Return the provided array to the free memory pool
   !
   ! The pointer provided will be set to null
   !
   subroutine decomp_pool_free_real(ptr)

      implicit none

      ! Arguments
      real(mytype), intent(inout), dimension(:, :, :), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      call decomp_pool%free(ptr)

   end subroutine decomp_pool_free_real

   subroutine decomp_pool_free_real1D(ptr)

      implicit none

      ! Arguments
      real(mytype), intent(inout), dimension(:), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      call decomp_pool%free(c_loc(ptr))
      nullify(ptr)

   end subroutine decomp_pool_free_real1D

   subroutine decomp_pool_free_cplx(ptr)

      implicit none

      ! Arguments
      complex(mytype), intent(inout), dimension(:, :, :), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      call decomp_pool%free(ptr)

   end subroutine decomp_pool_free_cplx

   subroutine decomp_pool_free_cplx1D(ptr)

      implicit none

      ! Arguments
      complex(mytype), intent(inout), dimension(:), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      call decomp_pool%free(c_loc(ptr))
      nullify(ptr)

   end subroutine decomp_pool_free_cplx1D

   !
   ! Update the default type in the module
   !
   subroutine decomp_pool_set_default_type(type)

      implicit none

      integer, intent(in) :: type

      ! Safety check
      if (type /= mem_pool_freal .and. &
          type /= mem_pool_dreal .and. &
          type /= mem_pool_fcplx .and. &
          type /= mem_pool_dcplx .and. &
          type /= mem_pool_int64 .and. &
          type /= mem_pool_int32 .and. &
          type /= mem_pool_int16 .and. &
          type /= mem_pool_int8 .and. &
          type /= mem_pool_bool) then
         call decomp_2d_abort(__FILE__, __LINE__, type, "Invalid value.")
      end if

      decomp_pool_default_type = type

   end subroutine decomp_pool_set_default_type

end module m_decomp_pool
