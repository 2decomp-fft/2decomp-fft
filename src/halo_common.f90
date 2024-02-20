!=======================================================================
! This is part of the 2DECOMP&FFT library
!
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil)
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2021 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This file contain common code to be included by subroutines
! 'update_halo_...' in halo.f90

    if (present(opt_global)) then
       global = opt_global
    else
       global = .false.
    end if

    s1 = size(in, 1)
    s2 = size(in, 2)
    s3 = size(in, 3)

    if (present(opt_pencil)) then
       ipencil = opt_pencil
    else
       ! Historic/default behaviour
       if (s1 == decomp%xsz(1)) then
          ipencil = 1
          if (first_call_x) then
             first_call_x = .false.
             call decomp_2d_warning(__FILE__, __LINE__, &
                                    0, "Deprecated interface - calling halo in X without explicit pencil")
          end if
       else if (s2 == decomp%ysz(2)) then
          ipencil = 2
          if (first_call_y) then
             first_call_y = .false.
             call decomp_2d_warning(__FILE__, __LINE__, &
                                    0, "Deprecated interface - calling halo in Y without explicit pencil")
          end if
       else if (s3 == decomp%zsz(3)) then
          ipencil = 3
          if (first_call_z) then
             first_call_z = .false.
             call decomp_2d_warning(__FILE__, __LINE__, &
                                    0, "Deprecated interface - calling halo in Z without explicit pencil")
          end if
       else
          ipencil = 0
          call decomp_2d_abort(__FILE__, __LINE__, 1, "Invalid decomposition size")
       end if
    end if

    ! Calculate the starting index and ending index of output
    if (ipencil == 1) then  ! X-pencil input
       if (global) then
          xs = decomp%xst(1)
          xe = decomp%xen(1)
          ys = decomp%xst(2) - level
          ye = decomp%xen(2) + level
          zs = decomp%xst(3) - level
          ze = decomp%xen(3) + level
       else
          xs = 1
          xe = s1
          ys = 1 - level
          ye = s2 + level
          zs = 1 - level
          ze = s3 + level
       end if
    else if (ipencil == 2) then  ! Y-pencil input
       if (global) then
          xs = decomp%yst(1) - level
          xe = decomp%yen(1) + level
          ys = decomp%yst(2)
          ye = decomp%yen(2)
          zs = decomp%yst(3) - level
          ze = decomp%yen(3) + level
       else
          xs = 1 - level
          xe = s1 + level
          ys = 1
          ye = s2
          zs = 1 - level
          ze = s3 + level
       end if
    else if (ipencil == 3) then  ! Z-pencil input
       if (global) then
          xs = decomp%zst(1) - level
          xe = decomp%zen(1) + level
          ys = decomp%zst(2) - level
          ye = decomp%zen(2) + level
          zs = decomp%zst(3)
          ze = decomp%zen(3)
       else
          xs = 1 - level
          xe = s1 + level
          ys = 1 - level
          ye = s2 + level
          zs = 1
          ze = s3
       end if
    else
       ! invalid input

       ! Set defaults to silence "uninitialised errors"
       xs = 1; xe = 1
       ys = 1; ye = 1
       zs = 1; ze = 1

       call decomp_2d_abort(__FILE__, __LINE__, 10, &
                            'Invalid data passed to update_halo')
    end if
    allocate (out(xs:xe, ys:ye, zs:ze))

    !    out = -1.0_mytype ! fill the halo for debugging

    !$acc enter data create(requests,neighbour)
    ! copy input data to output
    if (global) then
       ! using global coordinate
       ! note the input array passed in always has index starting from 1
       ! need to work out the corresponding global index
       if (ipencil == 1) then
          kst = decomp%xst(3); ken = decomp%xen(3)
          jst = decomp%xst(2); jen = decomp%xen(2)
          !$acc kernels default(present)
          do k = kst, ken
             do j = jst, jen
                do i = 1, s1  ! x all local
                   out(i, j, k) = in(i, j - decomp%xst(2) + 1, k - decomp%xst(3) + 1)
                end do
             end do
          end do
          !$acc end kernels
       else if (ipencil == 2) then
          kst = decomp%yst(3); ken = decomp%yen(3)
          ist = decomp%yst(1); ien = decomp%yen(1)
          !$acc kernels default(present)
          do k = kst, ken
             do j = 1, s2  ! y all local
                do i = ist, ien
                   out(i, j, k) = in(i - decomp%yst(1) + 1, j, k - decomp%yst(3) + 1)
                end do
             end do
          end do
          !$acc end kernels
       else if (ipencil == 3) then
          jst = decomp%zst(2); jen = decomp%zen(2)
          ist = decomp%xst(1); ien = decomp%xen(1)
          !$acc kernels default(present)
          do k = 1, s3  ! z all local
             do j = jst, jen
                do i = ist, ien
                   out(i, j, k) = in(i - decomp%zst(1) + 1, j - decomp%zst(2) + 1, k)
                end do
             end do
          end do
          !$acc end kernels
       end if
    else
       ! not using global coordinate
       !$acc kernels default(present)
       do k = 1, s3
          do j = 1, s2
             do i = 1, s1
                out(i, j, k) = in(i, j, k)
             end do
          end do
       end do
       !$acc end kernels
       !!! istat = cudaMemcpy(out,in,s1*s2*s3,cudaMemcpyDeviceToDevice)
    end if

    if (ipencil == 1) then
       call exchange_halo_x(out, opt_xlevel=(/0, level, level/))
    else if (ipencil == 2) then
       call exchange_halo_y(out, opt_ylevel=(/level, 0, level/))
    else if (ipencil == 3) then
       call exchange_halo_z(out, opt_zlevel=(/level, level, 0/))
    end if
