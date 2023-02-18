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
! This submodule allow interpolation between grids
!

submodule(decomp_2d) d2d_interp

   use decomp_2d_constants

   implicit none

contains

   !
   ! Small wrapper for 3D interpolation
   !
   module subroutine d2d_interp_var3d_1(ipencil, varin, varout, decompout, interp)

      implicit none

      ! Arguments
      integer, intent(in) :: ipencil
      type(decomp_info), intent(in) :: decompout
      real(mytype), dimension(:, :, :), intent(in) :: varin
      real(mytype), dimension(:, :, :), intent(out) :: varout
      integer, intent(in), optional :: interp

      call d2d_interp_var3d(ipencil, varin, decomp_main, varout, decompout, interp)

   end subroutine d2d_interp_var3d_1

   !
   ! Small wrapper for 3D interpolation
   !
   module subroutine d2d_interp_var3d_2(ipencil, varin, decompin, varout, interp)

      implicit none

      ! Arguments
      integer, intent(in) :: ipencil
      type(decomp_info), intent(in) :: decompin
      real(mytype), dimension(:, :, :), intent(in) :: varin
      real(mytype), dimension(:, :, :), intent(out) :: varout
      integer, intent(in), optional :: interp

      call d2d_interp_var3d(ipencil, varin, decompin, varout, decomp_main, interp)

   end subroutine d2d_interp_var3d_2

   !
   ! Interpolate the given 3D field
   !
   module subroutine d2d_interp_var3d_0(ipencil, varin, decompin, varout, decompout, interp)

      implicit none

      ! Arguments
      integer, intent(in) :: ipencil
      type(decomp_info), intent(in) :: decompin, decompout
      real(mytype), dimension(:, :, :), intent(in) :: varin
      real(mytype), dimension(:, :, :), intent(out) :: varout
      integer, intent(in), optional :: interp

      ! Local variable
      integer :: method

      if (present(interp)) then
         method = interp
      else
         method = DECOMP_2D_INTERP_BASIC
      endif

      ! Safety check
      if (ipencil < 1 .or. ipencil > 3) then
         call decomp_2d_abort(__FILE__, __LINE__, ipencil, "Invalid value for ipencil")
      endif

      ! Only one interpolation method currently
      if (method == DECOMP_2D_INTERP_BASIC) then
         call d2d_interp_var3d_basic(ipencil, varin, decompin, varout, decompout)
      else
         call decomp_2d_abort(__FILE__, __LINE__, method, "Invalid value for interp")
      endif

   end subroutine d2d_interp_var3d_0

   !
   ! Closest neighbour interpolation without transpose operations
   !
   subroutine d2d_interp_var3d_basic(ipencil, varin, decompin, varout, decompout)

      implicit none

      ! Arguments
      integer, intent(in) :: ipencil
      type(decomp_info), intent(in) :: decompin, decompout
      real(mytype), dimension(:, :, :), intent(in) :: varin
      real(mytype), dimension(:, :, :), intent(out) :: varout

      ! Local variables
      integer, dimension(3) :: totin, totout

      totin(1) = decompin%xsz(1)
      totin(2) = decompin%ysz(2)
      totin(3) = decompin%zsz(3)

      totout(1) = decompout%xsz(1)
      totout(2) = decompout%ysz(2)
      totout(3) = decompout%zsz(3)

      if (ipencil==1) then
         call interp_var3d_basic(varin, decompin%xsz, decompin%xst, totin, &
                                 varout, decompout%xsz, decompout%xst, totout)
      else if (ipencil==2) then
         call interp_var3d_basic(varin, decompin%ysz, decompin%yst, totin, &
                                 varout, decompout%ysz, decompout%yst, totout)
      else
         call interp_var3d_basic(varin, decompin%zsz, decompin%zst, totin, &
                                 varout, decompout%zsz, decompout%zst, totout)
      endif

   end subroutine d2d_interp_var3d_basic
   !
   subroutine interp_var3d_basic(varin, szin, stin, totin, &
                                 varout, szout, stout, totout)

      implicit none

      ! Arguments
      integer, dimension(3), intent(in) :: szin, stin, totin, szout, stout, totout
      real(mytype), dimension(szin(1), szin(2), szin(3)), intent(in) :: varin
      real(mytype), dimension(szout(1), szout(2), szout(3)), intent(out) :: varout

      ! Local variables
      integer :: i, j, k             ! local coordinate in varout
      integer :: iglob, jglob, kglob ! global coordinate in varout
      integer :: itarg, jtarg, ktarg ! local coordinate in varin

      ! If exact downsampling x2
      ! Keeps 2, 4, 6, ...
      ! Skip 1, 3, 5, ...
      do k = 1, szout(3)
         kglob = k + stout(3) - 1
         ktarg = nint(real(kglob*totin(3))/totout(3)) - stin(3) + 1
         ktarg = max(1, min(szin(3), ktarg))
         do j = 1, szout(2)
            jglob = j + stout(2) - 1
            jtarg = nint(real(jglob*totin(2))/totout(2)) - stin(2) + 1
            jtarg = max(1, min(szin(2), jtarg))
            do i = 1, szout(1)
               iglob = i + stout(1) - 1
               itarg = nint(real(iglob*totin(1))/totout(1)) - stin(1) + 1
               itarg = max(1, min(szin(1), itarg))
               varout(i, j, k) = varin(itarg, jtarg, ktarg)
            enddo
         enddo
      enddo

   end subroutine interp_var3d_basic

end submodule d2d_interp
