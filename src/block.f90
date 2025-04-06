!! SPDX-License-Identifier: BSD-3-Clause

!
! Low-level module for a list of memory blocks
!
! The external code will not use directly this object
!
module m_blk

   use iso_fortran_env, only: real32, output_unit, error_unit
   use iso_c_binding, only: c_ptr, c_null_ptr, c_size_t, c_associated, c_loc
   use decomp_2d_constants
   use decomp_2d_mpi, only : nrank, decomp_2d_abort

   implicit none

   !
   ! Derived type for a memory block inside a list
   !
   ! The external code will not use directly this object
   !
   type :: blk
      ! Memory block
      real(real32), allocatable, dimension(:), private :: dat
#ifdef _GPU
      attributes(device) :: dat
#endif
      ! Address of the memory block
      type(c_ptr), public :: ref = c_null_ptr
      ! True when the block is allocated
      logical, private :: allocated = .false.
      ! Previous block in the list
      type(blk), pointer, public :: prev => null()
      ! Next block in the list
      type(blk), pointer, public :: next => null()
   contains
      ! Initialize and allocate a new block in the list
      procedure :: new => blk_new
      ! Clear and free the block
      procedure :: fin => blk_fin
      ! Get one block (remove it from the list)
      procedure :: get => blk_get
      ! Put one block in the list
      procedure :: put => blk_put
      ! Resize the block
      procedure :: resize => blk_resize
      ! Find a block using the address
      procedure :: find => blk_find
      ! Print the setup of the block
      procedure :: print => blk_print
   end type blk

   ! Default : private
   private
   public :: blk

contains

   !
   ! Initialize and allocate a new block in the list
   !
   !   - size : number of elements (sinfle precision real)
   !   - init : flag to initialize the array to zero
   !
   subroutine blk_new(self, size, init)

      implicit none

      ! Arguments
      class(blk), target, intent(inout) :: self
      integer(c_size_t), intent(in) :: size
      logical, intent(in) :: init

      ! Local variables
      integer :: ierr
      integer(c_size_t) :: i
      type(blk), pointer :: ptr

      ! Save a pointer to the next block if needed
      if (associated(self%next)) then
         ptr => self%next
         nullify (self%next)
      else
         nullify (ptr)
      end if

      ! Create a new block
      allocate (self%next, stat=ierr)
      if (ierr /= 0) &
         call decomp_2d_abort(__FILE__, __LINE__, ierr, "Block creation failed")

      ! Allocate the memory and store the address
      allocate (self%next%dat(size), stat=ierr)
      if (ierr /= 0) &
         call decomp_2d_abort(__FILE__, __LINE__, ierr, "Allocation failed")
      self%next%ref = c_loc(self%next%dat)

      ! Initialize the memory if needed
      if (init) then
         !$acc parallel loop default(present)
         !$omp parallel do
         do i = 1_c_size_t, size
            self%next%dat(i) = 0._real32
         end do
         !$omp end parallel do
         !$acc end loop
      end if

      ! Tag the block
      self%next%allocated = .true.

      ! The new block (self%next) is located after the provided one
      self%next%prev => self

      ! If needed, put the saved block (ptr) after the new one (self%next)
      if (associated(ptr)) then
         self%next%next => ptr
         self%next%next%prev => self%next
         nullify (ptr)
      else
         nullify (self%next%next)
      end if

   end subroutine blk_new

   !
   ! Clear and finalize the provided block
   !
   subroutine blk_fin(self)

      implicit none

      class(blk), target, intent(inout) :: self

      ! Local variable
      integer :: ierr

      ! Safety check
      if (.not. self%allocated) &
         call decomp_2d_abort(__FILE__, __LINE__, 1, "Block must be allocated")

      ! Remove the block from the list
      call self%get()

      ! Tag the block
      self%allocated = .false.

      ! Free memory
      deallocate (self%dat, stat=ierr)
      if (ierr /= 0) &
         call decomp_2d_abort(__FILE__, __LINE__, 1, "Deallocation failed")
      self%ref = c_null_ptr

   end subroutine blk_fin

   !
   ! Get one block (remove it from the list)
   !
   subroutine blk_get(self)

      implicit none

      class(blk), target, intent(inout) :: self

      ! Safety check
      if (.not. self%allocated) &
         call decomp_2d_abort(__FILE__, __LINE__, 1, "Block must be allocated")

      ! Update the list
      if (associated(self%next)) then
         self%prev%next => self%next
         self%next%prev => self%prev
      else
         nullify (self%prev%next)
      end if

      ! Clear the previous and next pointers
      nullify (self%prev)
      nullify (self%next)

   end subroutine blk_get

   !
   ! Put one block (restore it in the list)
   !
   subroutine blk_put(self, head)

      implicit none

      class(blk), target, intent(inout) :: self
      type(blk), pointer, intent(inout) :: head

      ! If needed, move head%next to self%next
      if (associated(head%next)) then
         self%next => head%next
         self%next%prev => self
      end if
      ! Put the block after head
      head%next => self
      self%prev => head

   end subroutine blk_put

   !
   ! Resize the memory block
   !
   !   - size : number of elements (sinfle precision real)
   !   - init : flag to initialize the array to zero
   !
   subroutine blk_resize(self, size, init)

      implicit none

      ! Arguments
      class(blk), target, intent(inout) :: self
      integer(c_size_t), intent(in) :: size
      logical, intent(in) :: init

      ! Local variable
      integer :: ierr
      integer(c_size_t) :: i

      ! Safety check
      if (.not. self%allocated) &
         call decomp_2d_abort(__FILE__, __LINE__, 1, "Block must be allocated")

      ! Free memory
      deallocate (self%dat, stat=ierr)
      if (ierr /= 0) &
         call decomp_2d_abort(__FILE__, __LINE__, ierr, "Deallocation failed")

      ! Allocate again
      allocate (self%dat(size), stat=ierr)
      if (ierr /= 0) &
         call decomp_2d_abort(__FILE__, __LINE__, ierr, "Allocation failed")

      ! Update the address
      self%ref = c_loc(self%dat)

      ! Initialize the memory if needed
      if (init) then
         !$acc parallel loop default(present)
         !$omp parallel do
         do i = 1_c_size_t, size
            self%dat(i) = 0._real32
         end do
         !$omp end parallel do
         !$acc end loop
      end if

   end subroutine blk_resize

   !
   ! Find the block using the address
   !
   !   - self is the beginning of the list
   !   - ref is the address of the memory block
   !   - ptr is a pointer to the block located
   !
   subroutine blk_find(self, ref, ptr)

      implicit none

      ! Arguments
      class(blk), target, intent(in) :: self
      type(c_ptr), intent(in) :: ref
      type(blk), pointer, intent(out) :: ptr

      ! Start at head
      ptr => self

      ! Process the list
      do
         ! It is a match !
         if (c_associated(ptr%ref, ref)) return
         ! Next
         if (associated(ptr%next)) then
            ptr => ptr%next
         else
            exit
         end if
      end do

      call decomp_2d_abort(__FILE__, __LINE__, 1, "Address not available")

   end subroutine blk_find

   !
   ! Print the content of the block
   !
   subroutine blk_print(self, opt_unit, name)

      implicit none

      ! Arguments
      class(blk), intent(in) :: self
      integer, optional, intent(in) :: opt_unit
      character(*), optional, intent(in) :: name

      ! Local variable
      integer :: output

      ! If possible, use the given IO unit
      ! Default : only rank 0 will print
      if (present(opt_unit)) then
         output = opt_unit
      else if (nrank == 0) then
         output = output_unit
      else
         output = error_unit
      end if
      if (output == error_unit) return

      if (present(name)) write (output, *) name
      if (c_associated(self%ref)) write (output, *) "  Data at ", transfer(self%ref, 0_c_size_t)
      if (associated(self%prev)) then
         if (c_associated(self%prev%ref)) then
            write (output, *) "  Prev at ", transfer(self%prev%ref, 0_c_size_t)
         else
            write (output, *) "  Prev without data"
         end if
      else
         write (output, *) "  Prev is null"
      end if
      if (associated(self%next)) then
         if (c_associated(self%next%ref)) then
            write (output, *) "  Next at ", transfer(self%next%ref, 0_c_size_t)
         else
            write (output, *) "  Next without data"
         end if
      else
         write (output, *) "  Next is null"
      end if

   end subroutine blk_print

end module m_blk
