!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the divergence of a random field using
!   (1) global transposition
!   (2) halo-cell exchange
! The two method should give identical results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program halo_test

   use mpi

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_mpi
   use MPI
#if defined(_GPU)
   use cudafor
   use openacc
#endif

   implicit none

   integer, parameter :: nx_base = 65, ny_base = 48, nz_base = 33
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer :: nargin, arg, FNLength, status, DecInd
   character(len=80) :: InputFN

   real(mytype), allocatable, dimension(:, :, :) :: u1, v1, w1, div1
   real(mytype), allocatable, dimension(:, :, :) :: u2, v2, w2, div2
   real(mytype), allocatable, dimension(:, :, :) :: u3, v3, w3, div3
   real(mytype), allocatable, dimension(:, :, :) :: div, wk2, wk3

   integer :: i, j, k, ierror, n

   integer, allocatable, dimension(:) :: seed

   real(mytype) :: err
   integer :: xlast, ylast, zlast

   integer :: nx_expected, ny_expected, nz_expected

   logical :: passing, all_pass

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot/4) + 1
   nx = nx_base*resize_domain
   ny = ny_base*resize_domain
   nz = nz_base*resize_domain
   ! Now we can check if user put some inputs
   ! Handle input file like a boss -- GD
   nargin = command_argument_count()
   if ((nargin == 0) .or. (nargin == 2) .or. (nargin == 5)) then
      do arg = 1, nargin
         call get_command_argument(arg, InputFN, FNLength, status)
         read (InputFN, *, iostat=status) DecInd
         if (arg == 1) then
            p_row = DecInd
         elseif (arg == 2) then
            p_col = DecInd
         elseif (arg == 3) then
            nx = DecInd
         elseif (arg == 4) then
            ny = DecInd
         elseif (arg == 5) then
            nz = DecInd
         end if
      end do
   else
      ! nrank not yet computed we need to avoid write
      ! for every rank
      call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierror)
      if (nrank == 0) then
         print *, "This Test takes no inputs or 2 inputs as"
         print *, "  1) p_row (default=0)"
         print *, "  2) p_col (default=0)"
         print *, "or 5 inputs as"
         print *, "  1) p_row (default=0)"
         print *, "  2) p_col (default=0)"
         print *, "  3) nx "
         print *, "  4) ny "
         print *, "  5) nz "
         print *, "Number of inputs is not correct and the defult settings"
         print *, "will be used"
      end if
   end if
   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   xlast = xsize(1) - 1
   if (xend(2) == ny) then
      ylast = xsize(2) - 1
   else
      ylast = xsize(2)
   end if
   if (xend(3) == nz) then
      zlast = xsize(3) - 1
   else
      zlast = xsize(3)
   end if

   call initialise()
   !$acc data copyin(u1,v1,w1) create(u2,v2,w2,u3,v3,w3,div1,div2,div3,wk2,wk3) copy(div)
   call test_div_transpose()
   call test_div_haloX()
   call test_div_haloY()
   call test_div_haloZ()
   !$acc end data

   if (nrank == 0) then
      write (*, *) '-----------------------------------------------'
      write (*, *) "All pass: ", all_pass
      write (*, *) '==============================================='
   end if

   call finalize()

   call decomp_2d_finalize

   if (.not. all_pass) call decomp_2d_abort(1, "Error in halo_test")

   call MPI_FINALIZE(ierror)

contains

   !=====================================================================
   ! Initialize
   !=====================================================================
   subroutine initialise()

      use decomp_2d

      implicit none

#ifdef HALO_GLOBAL
      logical, parameter :: global = .true.
#else
      logical, parameter :: global = .false.
#endif

      ! initialise u,v,w with random numbers in X-pencil
      call alloc_x(u1, global)
      call alloc_x(v1, global)
      call alloc_x(w1, global)
      call alloc_x(div, global)
      call alloc_x(div1, global)
      call alloc_x(div2, global)
      call alloc_x(div3, global)

      call random_seed(size=n)
      allocate (seed(n))
      seed = nrank + 1
      call random_seed(put=seed)
      call random_number(u1)
      call random_number(v1)
      call random_number(w1)

      ! Working array used more than once
      call alloc_y(u2, global)
      call alloc_y(v2, global)
      call alloc_y(w2, global)
      call alloc_y(wk2, global)

      call alloc_z(u3, global)
      call alloc_z(v3, global)
      call alloc_z(w3, global)
      call alloc_z(wk3, global)

      all_pass = .true.

   end subroutine initialise

   !=====================================================================
   ! Finalize with deallocation of arrays
   !=====================================================================
   subroutine finalize()

      implicit none

      deallocate (u1, v1, w1)
      deallocate (u2, v2, w2)
      deallocate (u3, v3, w3)
      deallocate (wk2, wk3)

   end subroutine finalize
   !=====================================================================
   ! Calculate divergence using global transposition
   !=====================================================================
   subroutine test_div_transpose()

      implicit none

#ifdef HALO_GLOBAL
      logical, parameter :: global = .true.
#else
      logical, parameter :: global = .false.
#endif
      integer :: ifirst, ilast ! I loop start/end
      integer :: jfirst, jlast ! J loop start/end
      integer :: kfirst, klast ! K loop start/end

      ! du/dx calculated on X-pencil
#ifdef HALO_GLOBAL
      kfirst = xstart(3); klast = xend(3)
      jfirst = xstart(2); jlast = xend(2)
#else
      kfirst = 1; klast = xsize(3)
      jfirst = 1; jlast = xsize(2)
#endif
      ifirst = 2; ilast = xsize(1) - 1

      !$acc kernels default(present)
      div(:, :, :) = 0.0_mytype
      !$acc end kernels
      !$acc kernels default(present)
      do k = kfirst, klast
         do j = jfirst, jlast
            do i = ifirst, ilast
               div(i, j, k) = u1(i + 1, j, k) - u1(i - 1, j, k)
            end do
         end do
      end do
      !$acc end kernels

      ! dv/dy calculated on Y-pencil
#ifdef HALO_GLOBAL
      kfirst = ystart(3); klast = yend(3)
      ifirst = ystart(1); ilast = yend(1)
#else
      kfirst = 1; klast = ysize(3)
      ifirst = 1; ilast = ysize(1)
#endif
      jfirst = 2; jlast = ysize(2) - 1

      call transpose_x_to_y(v1, v2)
      call transpose_x_to_y(div, wk2)

      !$acc kernels default(present)
      do k = kfirst, klast
         do j = jfirst, jlast
            do i = ifirst, ilast
               wk2(i, j, k) = wk2(i, j, k) + v2(i, j + 1, k) - v2(i, j - 1, k)
            end do
         end do
      end do
      !$acc end kernels

      ! dw/dz calculated on Z-pencil
#ifdef HALO_GLOBAL
      jfirst = zstart(2); jlast = zend(2)
      ifirst = zstart(1); ilast = zend(1)
#else
      jfirst = 1; jlast = zsize(2)
      ifirst = 1; ilast = zsize(1)
#endif
      kfirst = 2; klast = zsize(3) - 1

      call transpose_x_to_y(w1, w2)
      call transpose_y_to_z(w2, w3)
      call transpose_y_to_z(wk2, wk3)

      !$acc kernels default(present)
      do k = kfirst, klast
         do j = jfirst, jlast
            do i = ifirst, ilast
               wk3(i, j, k) = wk3(i, j, k) + w3(i, j, k + 1) - w3(i, j, k - 1)
            end do
         end do
      end do
      !$acc end kernels

      ! result in X-pencil
      call transpose_z_to_y(wk3, wk2)
      call transpose_y_to_x(wk2, div)

#ifdef DEBUG
      if (nrank == 0) then
         write (*, *) 'Calculated via global transposition'
         !$acc update self(div)
         write (*, *) (div(i, i, i), i=2, 13)
      end if
#endif

   end subroutine test_div_transpose

   !=====================================================================
   ! Calculate divergence using halo-cell exchange (data in X-pencil)
   !=====================================================================
   subroutine test_div_haloX()

      implicit none

      real(mytype), allocatable, dimension(:, :, :) :: vh, wh
#if defined(_GPU)
      attributes(device) :: vh, wh
#endif

#ifdef HALO_GLOBAL
      logical, parameter :: global = .true.
#else
      logical, parameter :: global = .false.
#endif
      integer :: ifirst, ilast ! I loop start/end
      integer :: jfirst, jlast ! J loop start/end
      integer :: kfirst, klast ! K loop start/end

      ! Expected sizes
      nx_expected = nx
      ny_expected = xsize(2) + 2
      nz_expected = xsize(3) + 2

      ! Only global arrays defined in initialise needs to be ported
      ! Halo array are allocated in both host and device in update_halo
      ! Halo arrays are just removed before being deallocated
#ifdef HALO_GLOBAL
      call update_halo(v1, vh, 1, opt_global=.true., opt_pencil=1)
      call update_halo(w1, wh, 1, opt_global=.true., opt_pencil=1)

      kfirst = xstart(3); klast = xend(3)
      jfirst = xstart(2); jlast = xend(2)
#else
      call update_halo(v1, vh, 1, opt_pencil=1)
      call update_halo(w1, wh, 1, opt_pencil=1)

      kfirst = 1; klast = xsize(3)
      jfirst = 1; jlast = xsize(2)
#endif
      ifirst = 2; ilast = xsize(1) - 1

      call test_halo_size(vh, nx_expected, ny_expected, nz_expected, "X:v")
      call test_halo_size(wh, nx_expected, ny_expected, nz_expected, "X:w")

      !$acc kernels default(present)
      div1(:, :, :) = 0._mytype
      !$acc end kernels
      !$acc kernels default(present)
      do k = kfirst, klast
         do j = jfirst, jlast
            do i = ifirst, ilast
               div1(i, j, k) = (u1(i + 1, j, k) - u1(i - 1, j, k)) &
                               + (vh(i, j + 1, k) - vh(i, j - 1, k)) &
                               + (wh(i, j, k + 1) - wh(i, j, k - 1))
            end do
         end do
      end do
      !$acc end kernels

      ! Compute error
      call check_err(div1, "X")

      deallocate (vh, wh)

   end subroutine test_div_haloX

   !=====================================================================
   ! Calculate divergence using halo-cell exchange (data in Y-pencil)
   !=====================================================================
   subroutine test_div_haloY()

      implicit none
      real(mytype), allocatable, dimension(:, :, :) :: uh, wh
#if defined(_GPU)
      attributes(device) :: uh, wh
#endif

#ifdef HALO_GLOBAL
      logical, parameter :: global = .true.
#else
      logical, parameter :: global = .false.
#endif
      integer :: ifirst, ilast ! I loop start/end
      integer :: jfirst, jlast ! J loop start/end
      integer :: kfirst, klast ! K loop start/end

      ! Expected sizes
      nx_expected = ysize(1) + 2
      ny_expected = ny
      nz_expected = ysize(3) + 2

      call transpose_x_to_y(u1, u2)
      call transpose_x_to_y(v1, v2)
      call transpose_x_to_y(w1, w2)

      ! du/dx
#ifdef HALO_GLOBAL
      call update_halo(u2, uh, 1, opt_global=.true., opt_pencil=2)
      call update_halo(w2, wh, 1, opt_global=.true., opt_pencil=2)
      kfirst = ystart(3); klast = yend(3)
      ifirst = ystart(1); ilast = yend(1)
#else
      call update_halo(u2, uh, 1, opt_pencil=2)
      call update_halo(w2, wh, 1, opt_pencil=2)
      kfirst = 1; klast = ysize(3)
      ifirst = 1; ilast = ysize(1)
#endif
      jfirst = 2; jlast = ysize(2) - 1

      call test_halo_size(uh, nx_expected, ny_expected, nz_expected, "Y:u")
      call test_halo_size(wh, nx_expected, ny_expected, nz_expected, "Y:w")

      !$acc kernels default(present)
      do k = kfirst, klast
         do j = jfirst, jlast
            do i = ifirst, ilast
               wk2(i, j, k) = (uh(i + 1, j, k) - uh(i - 1, j, k)) &
                              + (v2(i, j + 1, k) - v2(i, j - 1, k)) &
                              + (wh(i, j, k + 1) - wh(i, j, k - 1))
            end do
         end do
      end do
      !$acc end kernels

      call transpose_y_to_x(wk2, div2)

      ! Compute error
      call check_err(div2, "Y")

      deallocate (uh, wh)

   end subroutine test_div_haloY

   !=====================================================================
   ! Calculate divergence using halo-cell exchange (data in Z-pencil)
   !=====================================================================
   subroutine test_div_haloZ()

      implicit none
      real(mytype), allocatable, dimension(:, :, :) :: uh, vh
#if defined(_GPU)
      attributes(device) :: vh, uh
#endif

#ifdef HALO_GLOBAL
      logical, parameter :: global = .true.
#else
      logical, parameter :: global = .false.
#endif
      integer :: ifirst, ilast ! I loop start/end
      integer :: jfirst, jlast ! J loop start/end
      integer :: kfirst, klast ! K loop start/end

      ! Expected sizes
      nx_expected = zsize(1) + 2
      ny_expected = zsize(2) + 2
      nz_expected = nz

      call transpose_y_to_z(u2, u3)
      call transpose_y_to_z(v2, v3)
      call transpose_y_to_z(w2, w3)

      ! du/dx
#ifdef HALO_GLOBAL
      call update_halo(u3, uh, 1, opt_global=.true., opt_pencil=3)
      call update_halo(v3, vh, 1, opt_global=.true., opt_pencil=3)
      ifirst = zstart(1); ilast = zend(1)
      jfirst = zstart(2); jlast = zend(2)
#else
      call update_halo(u3, uh, 1, opt_pencil=3)
      call update_halo(v3, vh, 1, opt_pencil=3)
      ifirst = 1; ilast = zsize(1)
      jfirst = 1; jlast = zsize(2)
#endif
      kfirst = 2; klast = zsize(3) - 1

      call test_halo_size(uh, nx_expected, ny_expected, nz_expected, "Z:u")
      call test_halo_size(vh, nx_expected, ny_expected, nz_expected, "Z:v")

      !$acc kernels default(present)
      do j = jfirst, jlast
         do i = ifirst, ilast
            do k = kfirst, klast
               wk3(i, j, k) = uh(i + 1, j, k) - uh(i - 1, j, k) &
                              + vh(i, j + 1, k) - vh(i, j - 1, k) &
                              + w3(i, j, k + 1) - w3(i, j, k - 1)
            end do
         end do
      end do
      !$acc end kernels

      call transpose_z_to_y(wk3, wk2)
      call transpose_y_to_x(wk2, div3)

      ! Compute error
      call check_err(div3, "Z")

      deallocate (uh, vh)
   end subroutine test_div_haloZ
   !=====================================================================
   ! Check the difference between halo and transpose divergence
   !=====================================================================
   subroutine check_err(divh, pencil)

      implicit none

      real(mytype), dimension(:, :, :), intent(in) :: divh
      character(len=*), intent(in) :: pencil
      real(mytype), dimension(:, :, :), allocatable :: tmp
      real(mytype) :: divmag, error
#if defined(_GPU)
      attributes(device) :: tmp
#endif
      ! XXX: The Intel compiler SEGFAULTs if the array difference is computed inplace
      !      i.e. mag(divh(2:xlast,2:ylast,2:zlast) - div1(2:xlast,2:ylast,2:zlast))
      !      causes a SEGFAULT. Explicitly computing the difference in a temporary
      !      array seems to be OK
      allocate (tmp(size(divh, 1), size(divh, 2), size(divh, 3)))

      !$acc kernels default(present)
      tmp(2:xlast, 2:ylast, 2:zlast) = divh(2:xlast, 2:ylast, 2:zlast) - div1(2:xlast, 2:ylast, 2:zlast)
      !$acc end kernels
      error = mag(tmp)
      !$acc kernels default(present)
      tmp(2:xlast, 2:ylast, 2:zlast) = div1(2:xlast, 2:ylast, 2:zlast)
      !$acc end kernels
      divmag = mag(tmp)

      if (err < epsilon(divmag)*divmag) then
         passing = .true.
      else
         passing = .false.
      end if
      all_pass = all_pass .and. passing

      if (nrank == 0) then
         write (*, *) '-----------------------------------------------'
         write (*, *) 'Calculated via halo exchange (data in '//pencil//'-pencil)'
#ifdef DEBUG
         write (*, *) (divh(i, i, i), i=2, 13)
#endif
         write (*, *) 'Error: ', err, '; Relative: ', err/divmag
         write (*, *) 'Pass: ', passing
      end if
      deallocate (tmp)

   end subroutine check_err
   !=====================================================================
   ! Compute the magnitude af the ayyays
   !=====================================================================
   real(mytype) function mag(a)

      implicit none

      real(mytype), dimension(:, :, :), intent(in) :: a
#if defined(_GPU)
      attributes(device) :: a
#endif

      real(mytype) :: lmag, gmag

      lmag = 0._mytype
      !$acc parallel loop default(present) collapse(3) reduction(+:lmag)
      do k = 2, zlast
         do j = 2, ylast
            do i = 2, xlast
               lmag = lmag + a(i, j, k)**2
            end do
         end do
      end do
      !$acc end parallel

      call MPI_Allreduce(lmag, gmag, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      if (ierror /= 0) then
         call decomp_2d_abort(__FILE__, __LINE__, ierror, &
                              "halo_test::mag::MPI_Allreduce")
      end if

      mag = sqrt(gmag/(nx - 2)/(ny - 2)/(nz - 2))

   end function mag
   !=====================================================================
   ! Check the dimensions of the halo arrays are the one expected
   !=====================================================================
   subroutine test_halo_size(arrh, nx_expected, ny_expected, nz_expected, tag)

      real(mytype), dimension(:, :, :), intent(in) :: arrh
#if defined(_GPU)
      attributes(device) :: arrh
#endif
      integer, intent(in) :: nx_expected, ny_expected, nz_expected
      character(len=*), intent(in) :: tag

      integer :: nx, ny, nz

      character(len=128) :: rank_lbl

      nx = size(arrh, 1)
      ny = size(arrh, 2)
      nz = size(arrh, 3)

      write (rank_lbl, "(A,I0,A)") "Rank", nrank, ":"

      if ((nx /= nx_expected) .or. &
          (ny /= ny_expected) .or. &
          (nz /= nz_expected)) then
         write (*, *) trim(rank_lbl), " ", tag, ":ERROR: halo size"
         write (*, *) trim(rank_lbl), " ", "+ Expected: ", nx_expected, " ", ny_expected, " ", nz_expected, " "
         write (*, *) trim(rank_lbl), " ", "+ Got:      ", nx, " ", ny, " ", nz, " "

         all_pass = .false.
      else
         write (*, *) trim(rank_lbl), " ", tag, ":PASS"
      end if

   end subroutine test_halo_size

end program halo_test
