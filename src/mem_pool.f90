!! SPDX-License-Identifier: BSD-3-Clause

!
! Provide a memory pool
!
module m_mem_pool

   use iso_fortran_env, only: real32, real64, int64, int32, int16, int8, &
                              output_unit, error_unit
   use iso_c_binding, only: c_size_t, c_loc, c_associated, c_f_pointer, c_ptr, c_null_ptr
   use decomp_2d_constants
   use decomp_2d_mpi, only : nrank, decomp_2d_abort
   use m_blk
   use m_info
   use mpi
#ifdef _GPU
   use cudafor
#endif

   implicit none

   !
   ! Maintain a list of free memory blocks with uniform size
   !
   type :: mem_pool
      ! Size of the blocks in floats
      integer(c_size_t), private :: size
      ! Head of the lists
      type(blk), private, pointer :: free_head, busy_head
      type(blk), dimension(:), allocatable, private :: blk_target
      ! True when the list is active
      logical, private :: available = .false.
      ! Default shapes for the 3D blocks
      integer, dimension(3, 4), public :: shapes
   contains
      ! Initialize the memory pool
      procedure :: init => mem_pool_init
      ! Finalize the memory pool
      procedure :: fin => mem_pool_fin
      ! Create a block inside the free memory pool
      procedure :: new => mem_pool_new
      ! Deallocate free blocks to reduce memory usage
      procedure :: purge => mem_pool_purge
      ! Get a block from the free memory pool
      procedure :: get_raw => mem_pool_get_raw
      generic :: get => get_freal, get_dreal, get_fcplx, get_dcplx, &
                        get_int64, get_int32, get_int16, get_int8, get_bool
      procedure, private :: get_freal => mem_pool_get_freal
      procedure, private :: get_dreal => mem_pool_get_dreal
      procedure, private :: get_fcplx => mem_pool_get_fcplx
      procedure, private :: get_dcplx => mem_pool_get_dcplx
      procedure, private :: get_int64 => mem_pool_get_int64
      procedure, private :: get_int32 => mem_pool_get_int32
      procedure, private :: get_int16 => mem_pool_get_int16
      procedure, private :: get_int8 => mem_pool_get_int8
      procedure, private :: get_bool => mem_pool_get_bool
      ! Return a block to the free memory pool
      generic :: free => free_raw, free_freal, free_dreal, free_fcplx, free_dcplx, &
                         free_int64, free_int32, free_int16, free_int8, free_bool
      procedure, private :: free_raw => mem_pool_free_raw
      procedure, private :: free_freal => mem_pool_free_freal
      procedure, private :: free_dreal => mem_pool_free_dreal
      procedure, private :: free_fcplx => mem_pool_free_fcplx
      procedure, private :: free_dcplx => mem_pool_free_dcplx
      procedure, private :: free_int64 => mem_pool_free_int64
      procedure, private :: free_int32 => mem_pool_free_int32
      procedure, private :: free_int16 => mem_pool_free_int16
      procedure, private :: free_int8 => mem_pool_free_int8
      procedure, private :: free_bool => mem_pool_free_bool
      ! Add a new shape to the memory pool
      procedure :: new_shape => mem_pool_new_shape
      ! Return the size of the memory pool
      procedure :: get_size => mem_pool_get_size
      ! Resize the memory pool
      procedure :: resize => mem_pool_resize
      ! Print the setup of the memory pool
      procedure :: print => mem_pool_print
   end type mem_pool

   integer, parameter :: mem_pool_none = 0, &
                         mem_pool_freal = MPI_REAL, &             ! compatible with 2decomp real_type
                         mem_pool_dreal = MPI_DOUBLE_PRECISION, & ! compatible with 2decomp real_type
                         mem_pool_fcplx = MPI_COMPLEX, &          ! compatible with 2decomp complex_type
                         mem_pool_dcplx = MPI_DOUBLE_COMPLEX, &   ! compatible with 2decomp complex_type
                         mem_pool_int64 = MPI_INT64_T, &
                         mem_pool_int32 = MPI_INT32_T, &
                         mem_pool_int16 = MPI_INT16_T, &
                         mem_pool_int8 = MPI_INT8_T, &
                         mem_pool_bool = MPI_C_BOOL

   ! Default type in the memory pool module : mem_pool_freal or mem_pool_dreal
   !
   ! The external code can change it using mem_pool_set_default_type
   integer, protected, save :: mem_pool_default_type = real_type

   ! Default : private
   private
   public :: mem_pool, &
             mem_pool_freal, &
             mem_pool_dreal, &
             mem_pool_fcplx, &
             mem_pool_dcplx, &
             mem_pool_int64, &
             mem_pool_int32, &
             mem_pool_int16, &
             mem_pool_int8, &
             mem_pool_bool, &
             mem_pool_default_type, &
             mem_pool_set_default_type

contains

   !
   ! Initialize the memory pool
   !
   !   - decomp or shape provides the type and the size of the 3D array
   !   - blk_n provides the number of blocks to allocate
   !   - blk_init flag to initialize the block
   !
   subroutine mem_pool_init(self, decomp, shape, blk_n, blk_init)

      implicit none

      ! Arguments
      class(mem_pool), target, intent(inout) :: self
      class(info), intent(in), optional :: decomp
      integer, intent(in), optional :: shape(:)
      integer, intent(in), optional :: blk_n
      logical, intent(in), optional :: blk_init

      ! Local variables
      integer :: i

      ! Safety check
      if (self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool already initialised")
      if (present(shape)) then
         if (size(shape) < 3 .or. size(shape) > 4) &
            call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument shape")
      end if

      ! Tag the memory pool
      self%available = .true.

      ! One list for the free blocks, one for the busy ones
      allocate (self%blk_target(2))
      self%free_head => self%blk_target(1)
      self%busy_head => self%blk_target(2)

      ! Compute the size and the default shape
      self%size = 0_c_size_t
      self%shapes = mem_pool_none
      if (present(decomp)) then
         call update_size_shapes(self%size, self%shapes, mem_pool_default_type, decomp%xsz)
         call update_size_shapes(self%size, self%shapes, mem_pool_default_type, decomp%ysz)
         call update_size_shapes(self%size, self%shapes, mem_pool_default_type, decomp%zsz)
      else if (present(shape)) then
         if (size(shape) == 3) then
            call update_size_shapes(self%size, self%shapes, mem_pool_default_type, shape(:))
         else
            call update_size_shapes(self%size, self%shapes, shape(1), shape(2:))
         end if
      else
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid arguments")
      end if

      ! Create free blocks if needed
      if (present(blk_n)) then
         ! Safety check
         if (blk_n < 0) &
            call decomp_2d_abort(__FILE__, __LINE__, blk_n, "Invalid value in argument blk_n")
         if (present(blk_init)) then
            do i = 1, blk_n
               call self%new(init=blk_init)
            end do
         else
            do i = 1, blk_n
               call self%new(init=.false.)
            end do
         end if
      end if

   end subroutine mem_pool_init

   !
   ! Finalize the memory pool
   !
   subroutine mem_pool_fin(self)

      implicit none

      ! Arguments
      class(mem_pool), target, intent(inout) :: self

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Error if some blocks are busy
      if (associated(self%busy_head%next)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "All blocks are not released")

      ! Clear free blocks
      call self%purge()

      ! Tag the memory pool
      self%available = .false.

      ! Restore the size and the default shapes
      self%size = 0_c_size_t
      self%shapes = mem_pool_none

      ! Free the local blocks
      nullify (self%free_head)
      nullify (self%busy_head)
      deallocate (self%blk_target)

   end subroutine mem_pool_fin

   !
   ! Add a free memory block
   !
   subroutine mem_pool_new(self, init)

      implicit none

      ! Arguments
      class(mem_pool), target, intent(inout) :: self
      logical, intent(in) :: init

      ! Local variable
#ifdef _GPU
      type(blk), pointer, device :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Add a block at the beginning of the free list
#ifdef _GPU
      ptr => self%free_head
      call ptr%new(self%size, init)
      nullify(ptr)
#else
      call self%free_head%new(self%size, init)
#endif

   end subroutine mem_pool_new

   !
   ! Deallocate and destroy free blocks to reduce memory usage
   !
   subroutine mem_pool_purge(self)

      implicit none

      ! Argument
      class(mem_pool), intent(inout) :: self

      ! Local variable
      type(blk), pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Return if the memory pool was not initialized
      if (.not. self%available) return

      ! Clear free blocks
      if (associated(self%free_head%next)) then
         ptr => self%free_head%next
         do
            call ptr%fin
            deallocate (ptr)
            if (.not. associated(self%free_head%next)) exit
            ptr => self%free_head%next
         end do
         nullify (ptr)
      end if

   end subroutine mem_pool_purge

   !
   ! Get a block from the free memory list
   !
   ! If no block is available, a new one is created
   !
   function mem_pool_get_raw(self, init) result(ptr)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      logical, intent(in) :: init
#ifdef _GPU
      type(c_devptr) :: ptr
#else
      type(c_ptr) :: ptr
#endif

      ! Local variable
      type(blk), pointer :: blk_ptr
#ifdef _GPU
      attributes(device) :: blk_ptr
      type(blk), pointer, device :: blk_ptr2
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! If no block is available, create a free one
      if (.not. associated(self%free_head%next)) call self%new(init)

      blk_ptr => self%free_head%next
      call blk_ptr%get()
#ifdef _GPU
      blk_ptr2 => self%busy_head
      call blk_ptr%put(blk_ptr2)
      nullify(blk_ptr2)
#else
      call blk_ptr%put(self%busy_head)
#endif
      nullify (blk_ptr)

      ! Done
      ptr = self%busy_head%next%ref

   end function mem_pool_get_raw

   subroutine mem_pool_get_freal(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      real(real32), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in), optional :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Local variables
      integer, dimension(:), allocatable :: shp

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape, or the default one
      if (present(shape)) then
         shp = shape
         if (product(int(shp, c_size_t)) > self%size) &
            call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape")
      else if (self%shapes(1, 1) /= mem_pool_none) then
         shp = self%shapes(:, 1)
      else
         call decomp_2d_abort(__FILE__, __LINE__, 2, "No shape available")
      end if

      ! Get a block and associate it
      call c_f_pointer(self%get_raw(.false.), ptr, shp)

      ! Free memory
      deallocate (shp)

   end subroutine mem_pool_get_freal

   subroutine mem_pool_get_dreal(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      real(real64), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in), optional :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Local variables
      integer, dimension(:), allocatable :: shp

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape, or the default one
      if (present(shape)) then
         shp = shape
         if (2_c_size_t*product(int(shp, c_size_t)) > self%size) &
            call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape")
      else if (self%shapes(1, 2) /= mem_pool_none) then
         shp = self%shapes(:, 2)
      else
         call decomp_2d_abort(__FILE__, __LINE__, 2, "No shape available")
      end if

      ! Get a block and associate it
      call c_f_pointer(self%get_raw(.false.), ptr, shp)

      ! Free memory
      deallocate (shp)

   end subroutine mem_pool_get_dreal

   subroutine mem_pool_get_fcplx(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      complex(real32), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in), optional :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Local variables
      integer, dimension(:), allocatable :: shp

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape, or the default one
      if (present(shape)) then
         shp = shape
         if (2_c_size_t*product(int(shp, c_size_t)) > self%size) &
            call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape")
      else if (self%shapes(1, 3) /= mem_pool_none) then
         shp = self%shapes(:, 3)
      else
         call decomp_2d_abort(__FILE__, __LINE__, 2, "No shape available")
      end if

      ! Get a block and associate it
      call c_f_pointer(self%get_raw(.false.), ptr, shp)

      ! Free memory
      deallocate (shp)

   end subroutine mem_pool_get_fcplx

   subroutine mem_pool_get_dcplx(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      complex(real64), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in), optional :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Local variables
      integer, dimension(:), allocatable :: shp

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape, or the default one
      if (present(shape)) then
         shp = shape
         if (4_c_size_t*product(int(shp, c_size_t)) > self%size) &
            call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape")
      else if (self%shapes(1, 4) /= mem_pool_none) then
         shp = self%shapes(:, 4)
      else
         call decomp_2d_abort(__FILE__, __LINE__, 2, "No shape available")
      end if

      ! Output is ready
      call c_f_pointer(self%get_raw(.false.), ptr, shp)

      ! Free memory
      deallocate (shp)

   end subroutine mem_pool_get_dcplx

   subroutine mem_pool_get_int64(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int64), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in) :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape
      if (2_c_size_t*product(int(shape, c_size_t)) > self%size) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape.")

      ! Output is ready
      call c_f_pointer(self%get_raw(.false.), ptr, shape)

   end subroutine mem_pool_get_int64

   subroutine mem_pool_get_int32(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int32), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in) :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape
      if (product(int(shape, c_size_t)) > self%size) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape.")

      ! Output is ready
      call c_f_pointer(self%get_raw(.false.), ptr, shape)

   end subroutine mem_pool_get_int32

   subroutine mem_pool_get_int16(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int16), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in) :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape
      if (product(int(shape, c_size_t)) > self%size * 2_c_size_t) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape.")

      ! Output is ready
      call c_f_pointer(self%get_raw(.false.), ptr, shape)

   end subroutine mem_pool_get_int16

   subroutine mem_pool_get_int8(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int8), intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in) :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape
      if (product(int(shape, c_size_t)) > self%size * 4_c_size_t) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape.")

      ! Output is ready
      call c_f_pointer(self%get_raw(.false.), ptr, shape)

   end subroutine mem_pool_get_int8

   subroutine mem_pool_get_bool(self, ptr, shape)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      logical, intent(out), dimension(:, :, :), contiguous, pointer :: ptr
      integer, intent(in) :: shape(:)
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")

      ! Use the provided shape
      if (product(int(shape, c_size_t)) > self%size * 32_c_size_t) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid shape.")

      ! Output is ready
      call c_f_pointer(self%get_raw(.false.), ptr, shape)

   end subroutine mem_pool_get_bool

   !
   ! Return the provided array to the free memory pool
   !
   ! If a pointer is provided, it will be set to null
   !
   subroutine mem_pool_free_raw(self, raw)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
#ifdef _GPU
      type(c_devptr), intent(in) :: raw
#else
      type(c_ptr), intent(in) :: raw
#endif

      ! Local variables
      type(blk), pointer :: block
#ifdef _GPU
      attributes(device) :: block
      type(blk), pointer, device :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
#ifndef _GPU
      if (.not. c_associated(raw)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")
#endif

      ! Locate and extract the block in the busy list
#ifdef _GPU
      ptr => self%busy_head
      call ptr%find(raw, block)
      nullify(ptr)
#else
      call self%busy_head%find(raw, block)
#endif
      call block%get()

      ! Put it in the free list
#ifdef _GPU
      ptr => self%free_head
      call block%put(ptr)
      nullify(ptr)
#else
      call block%put(self%free_head)
#endif

      ! Free memory
      nullify (block)

   end subroutine mem_pool_free_raw

   subroutine mem_pool_free_freal(self, freal)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      real(real32), intent(inout), dimension(:, :, :), contiguous, pointer :: freal
#ifdef _GPU
      attributes(device) :: freal
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(freal)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(freal))
#else
      call self%free_raw(c_loc(freal))
#endif

      ! Free memory
      nullify (freal)

   end subroutine mem_pool_free_freal

   subroutine mem_pool_free_dreal(self, dreal)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      real(real64), intent(inout), dimension(:, :, :), contiguous, pointer :: dreal
#ifdef _GPU
      attributes(device) :: dreal
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(dreal)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(dreal))
#else
      call self%free_raw(c_loc(dreal))
#endif

      ! Free memory
      nullify (dreal)

   end subroutine mem_pool_free_dreal

   subroutine mem_pool_free_fcplx(self, fcplx)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      complex(real32), intent(inout), dimension(:, :, :), contiguous, pointer :: fcplx
#ifdef _GPU
      attributes(device) :: fcplx
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(fcplx)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(fcplx))
#else
      call self%free_raw(c_loc(fcplx))
#endif

      ! Free memory
      nullify (fcplx)

   end subroutine mem_pool_free_fcplx

   subroutine mem_pool_free_dcplx(self, dcplx)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      complex(real64), intent(inout), dimension(:, :, :), contiguous, pointer :: dcplx
#ifdef _GPU
      attributes(device) :: dcplx
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(dcplx)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(dcplx))
#else
      call self%free_raw(c_loc(dcplx))
#endif

      ! Free memory
      nullify (dcplx)

   end subroutine mem_pool_free_dcplx

   subroutine mem_pool_free_int64(self, ptr)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int64), intent(inout), dimension(:, :, :), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(ptr)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(ptr))
#else
      call self%free_raw(c_loc(ptr))
#endif

      ! Free memory
      nullify (ptr)

   end subroutine mem_pool_free_int64

   subroutine mem_pool_free_int32(self, ptr)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int32), intent(inout), dimension(:, :, :), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(ptr)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(ptr))
#else
      call self%free_raw(c_loc(ptr))
#endif

      ! Free memory
      nullify (ptr)

   end subroutine mem_pool_free_int32

   subroutine mem_pool_free_int16(self, ptr)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int16), intent(inout), dimension(:, :, :), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(ptr)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(ptr))
#else
      call self%free_raw(c_loc(ptr))
#endif

      ! Free memory
      nullify (ptr)

   end subroutine mem_pool_free_int16

   subroutine mem_pool_free_int8(self, ptr)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer(int8), intent(inout), dimension(:, :, :), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(ptr)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(ptr))
#else
      call self%free_raw(c_loc(ptr))
#endif

      ! Free memory
      nullify (ptr)

   end subroutine mem_pool_free_int8

   subroutine mem_pool_free_bool(self, ptr)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      logical, intent(inout), dimension(:, :, :), contiguous, pointer :: ptr
#ifdef _GPU
      attributes(device) :: ptr
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (.not. associated(ptr)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")

      ! Release the block
#ifdef _GPU
      call self%free_raw(c_devloc(ptr))
#else
      call self%free_raw(c_loc(ptr))
#endif

      ! Free memory
      nullify (ptr)

   end subroutine mem_pool_free_bool

   !
   ! Add a new shape to the memory pool
   !
   subroutine mem_pool_new_shape(self, type, decomp, shp)

      implicit none

      ! Arguments
      class(mem_pool), intent(inout) :: self
      integer, intent(in) :: type
      class(info), intent(in), optional :: decomp
      integer, intent(in), optional :: shp(:)

      ! Local variables
      integer(c_size_t) :: new_size, fact

      ! Compute scaling factor
      if (type == mem_pool_freal) then
         fact = 1_c_size_t
      else if (type == mem_pool_dreal) then
         fact = 2_c_size_t
      else if (type == mem_pool_fcplx) then
         fact = 2_c_size_t
      else if (type == mem_pool_dcplx) then
         fact = 4_c_size_t
      else if (type == mem_pool_int64) then
         fact = 2_c_size_t
      else if (type == mem_pool_int32) then
         fact = 1_c_size_t
      else if (type == mem_pool_int16 .or. &
               type == mem_pool_int8 .or. &
               type == mem_pool_bool) then
         fact = 1_c_size_t ! not optimal. FIXME ?
      else
         fact = 0_c_size_t
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid argument")
      end if

      ! Resize if needed
      if (present(decomp)) then
         new_size = fact*product(int(decomp%xsz, kind=c_size_t))
         new_size = max(new_size, fact*product(int(decomp%ysz, kind=c_size_t)))
         new_size = max(new_size, fact*product(int(decomp%zsz, kind=c_size_t)))
      else if (present(shp)) then
         new_size = fact*product(int(shp, kind=c_size_t))
      else
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Invalid arguments")
      end if
      if (new_size > self%size) call self%resize(new_size)

      ! Update default shape if possible
      if (present(decomp)) then
         if (type == mem_pool_freal) then
            if (any(self%shapes(:, 1) == mem_pool_none)) self%shapes(:, 1) = decomp%xsz
         else if (type == mem_pool_dreal) then
            if (any(self%shapes(:, 2) == mem_pool_none)) self%shapes(:, 2) = decomp%xsz
         else if (type == mem_pool_fcplx) then
            if (any(self%shapes(:, 3) == mem_pool_none)) self%shapes(:, 3) = decomp%xsz
         else if (type == mem_pool_dcplx) then
            if (any(self%shapes(:, 4) == mem_pool_none)) self%shapes(:, 4) = decomp%xsz
         end if
      else if (present(shp)) then
         if (size(shp) /= 3) return
         if (type == mem_pool_freal) then
            if (any(self%shapes(:, 1) == mem_pool_none)) self%shapes(:, 1) = shp
         else if (type == mem_pool_dreal) then
            if (any(self%shapes(:, 2) == mem_pool_none)) self%shapes(:, 2) = shp
         else if (type == mem_pool_fcplx) then
            if (any(self%shapes(:, 3) == mem_pool_none)) self%shapes(:, 3) = shp
         else if (type == mem_pool_dcplx) then
            if (any(self%shapes(:, 4) == mem_pool_none)) self%shapes(:, 4) = shp
         end if
      end if

   end subroutine mem_pool_new_shape

   !
   ! Return the size of the memory blocks
   !
   integer(c_size_t) function mem_pool_get_size(self)

      implicit none

      class(mem_pool), intent(in) :: self

      mem_pool_get_size = self%size

   end function mem_pool_get_size

   !
   ! Resize the memory pool and all the blocks
   !
   subroutine mem_pool_resize(self, size)

      implicit none

      ! Arguments
      class(mem_pool), target, intent(inout) :: self
      integer(c_size_t), intent(in) :: size

      ! Local variable
      type(blk), pointer :: tmp
#ifdef _GPU
      attributes(device) :: tmp
#endif

      ! Safety check
      if (.not. self%available) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Memory pool must be initialized")
      if (associated(self%busy_head%next)) &
         call decomp_2d_abort(__FILE__, __LINE__, 2, "Release all blocks to resize")

      ! Resize all free blocks
      self%size = size
      if (associated(self%free_head%next)) then
         tmp => self%free_head%next
         do
            call tmp%resize(size, .false.)
            if (.not. associated(tmp%next)) exit
            tmp => tmp%next
         end do
         nullify (tmp)
      end if

   end subroutine mem_pool_resize

   !
   ! Print the content of the list
   !
   subroutine mem_pool_print(self, opt_unit, name)

      implicit none

      ! Arguments
      class(mem_pool), intent(in) :: self
      integer, optional, intent(in) :: opt_unit
      character(*), optional, intent(in) :: name

      ! Local variable
      integer :: output
      type(blk), pointer :: tmp
#ifdef _GPU
      attributes(device) :: tmp
#endif

      ! If possible, use the given IO unit
      if (present(opt_unit)) then
         output = opt_unit
      else if (nrank == 0) then
         output = output_unit
      else
         output = error_unit
      end if
      if (output == error_unit) return

      if (present(name)) write (output, *) name

      write (output, *) " Size of the blocks ", self%size

      if (self%available) then
         write (output, *) " Memory pool active"
      else
         write (output, *) " Memory pool inactive"
         return
      end if

      write (output, *) " Default shape :"
      if (self%shapes(1, 1) == mem_pool_none) then
         write (output, *) " - Single prec, real : none"
      else
         write (output, *) " - Single prec, real : ", self%shapes(:, 1)
      end if
      if (self%shapes(1, 2) == mem_pool_none) then
         write (output, *) " - Double prec, real : none"
      else
         write (output, *) " - Double prec, real : ", self%shapes(:, 2)
      end if
      if (self%shapes(1, 3) == mem_pool_none) then
         write (output, *) " - Single prec, cplx : none"
      else
         write (output, *) " - Single prec, cplx : ", self%shapes(:, 3)
      end if
      if (self%shapes(1, 4) == mem_pool_none) then
         write (output, *) " - Double prec, cplx : none"
      else
         write (output, *) " - Double prec, cplx : ", self%shapes(:, 4)
      end if

      if (associated(self%free_head)) then
         if (associated(self%free_head%next)) then
            write (output, *) " List of free blocks :"
            tmp => self%free_head%next
            do
               write (output, *) " - ", transfer(tmp%ref, 0_c_size_t)
               tmp => tmp%next
               if (.not. associated(tmp)) exit
            end do
            nullify (tmp)
         else
            write (output, *) " List of free blocks is empty"
         end if
      end if

      if (associated(self%busy_head)) then
         if (associated(self%busy_head%next)) then
            write (output, *) " List of active blocks :"
            tmp => self%busy_head%next
            do
               write (output, *) " - ", transfer(tmp%ref, 0_c_size_t)
               tmp => tmp%next
               if (.not. associated(tmp)) exit
            end do
            nullify (tmp)
         else
            write (output, *) " List of active blocks is empty"
         end if
      end if

   end subroutine mem_pool_print

   !
   ! Update the default type in the module
   !
   subroutine mem_pool_set_default_type(type)

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

      mem_pool_default_type = type

   end subroutine mem_pool_set_default_type

   !
   ! Update the size and the shape
   !
   subroutine update_size_shapes(size, shapes, type, shape)

      implicit none

      ! Argument
      integer(c_size_t), intent(inout) :: size
      integer, dimension(:, :), intent(inout) :: shapes
      integer, intent(in) :: type, shape(:)

      ! Local variables
      integer(c_size_t) :: new_size, fact

      ! Compute scaling factor and update default shape
      if (type == mem_pool_freal) then
         fact = 1_c_size_t
         if (any(shapes(:, 1) == mem_pool_none)) shapes(:, 1) = shape
      else if (type == mem_pool_dreal) then
         fact = 2_c_size_t
         if (any(shapes(:, 2) == mem_pool_none)) shapes(:, 2) = shape
      else if (type == mem_pool_fcplx) then
         fact = 2_c_size_t
         if (any(shapes(:, 3) == mem_pool_none)) shapes(:, 3) = shape
      else if (type == mem_pool_dcplx) then
         fact = 4_c_size_t
         if (any(shapes(:, 4) == mem_pool_none)) shapes(:, 4) = shape
      else if (type == mem_pool_int64) then
         fact = 2_c_size_t
      else if (type == mem_pool_int32) then
         fact = 1_c_size_t
      else if (type == mem_pool_int16 .or. &
               type == mem_pool_int8 .or. &
               type == mem_pool_bool) then
         fact = 1_c_size_t ! not optimal. FIXME ?
      else
         fact = 0_c_size_t
         call decomp_2d_abort(__FILE__, __LINE__, type, "Invalid argument type")
      end if
      new_size = fact*product(int(shape, kind=c_size_t))

      ! Update size
      if (new_size > size) size = new_size

   end subroutine update_size_shapes

end module m_mem_pool
