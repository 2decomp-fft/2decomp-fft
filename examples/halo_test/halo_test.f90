!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the divergence of a random field using
!   (1) global transposition
!   (2) halo-cell exchange
! The two method should give identical results
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program halo_test

  use mpi

  use decomp_2d

  implicit none

  integer, parameter :: nx=171, ny=132, nz=113
  integer :: p_row=0, p_col=0

  real(mytype), allocatable, dimension(:,:,:) :: u1, u2, u3
  real(mytype), allocatable, dimension(:,:,:) :: v1, v2, v3
  real(mytype), allocatable, dimension(:,:,:) :: w1, w2, w3
  real(mytype), allocatable, dimension(:,:,:) :: wk2, wk3
  real(mytype), allocatable, dimension(:,:,:) :: uh, vh, wh
  real(mytype), allocatable, dimension(:,:,:) :: div1, div2, div3, div4

  integer :: i,j,k, ierror, n

  integer, allocatable, dimension(:) :: seed

  real(mytype) :: err, err_local
  integer :: xlast, ylast, zlast

  integer :: nx_expected, ny_expected, nz_expected
  
  logical :: passing, all_pass

  call MPI_INIT(ierror)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)

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
  call test_div_transpose()
  call test_div_haloX()
  call test_div_haloY()
  call test_div_haloZ()

  if (nrank == 0) then
     write(*,*) '-----------------------------------------------'
     write(*,*) "All pass: ", all_pass
     write(*,*) '==============================================='
  end if

  deallocate(u1,v1,w1,u2,v2,w2,u3,v3,w3)
  deallocate(div1,div2,div3,div4)

  call decomp_2d_finalize 

  if (.not. all_pass) call decomp_2d_abort(1, "Error in halo_test")

  call MPI_FINALIZE(ierror)

contains

  subroutine initialise()

    ! initialise u,v,w with random numbers in X-pencil
#ifdef HALO_GLOBAL
    allocate(u1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(v1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(w1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
#else
    allocate(u1(xsize(1), xsize(2), xsize(3)))
    allocate(v1(xsize(1), xsize(2), xsize(3)))
    allocate(w1(xsize(1), xsize(2), xsize(3)))
#endif

    call random_seed(size = n)
    allocate(seed(n))
    seed = nrank
    call random_seed(put=seed)
    call random_number(u1)
    call random_number(v1)
    call random_number(w1)

    all_pass = .true.

  end subroutine initialise

  !=====================================================================
  ! Calculate divergence using global transposition
  !=====================================================================
  subroutine test_div_transpose()

    integer :: i1, in ! I loop start/end
    integer :: j1, jn ! J loop start/end
    integer :: k1, kn ! K loop start/end

    ! du/dx calculated on X-pencil
#ifdef HALO_GLOBAL
    allocate(div1(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    k1 = xstart(3); kn = xend(3)
    j1 = xstart(2); jn = xend(2)
#else
    allocate(div1(xsize(1), xsize(2), xsize(3)))
    k1 = 1; kn = xsize(3)
    j1 = 1; jn = xsize(2)
#endif
    i1 = 2; in = xsize(1) - 1

    div1 = 0.0_mytype
    do k=k1,kn
       do j=j1,jn
          do i=i1,in
             div1(i,j,k) = u1(i+1,j,k)-u1(i-1,j,k)
          end do
       end do
    end do

    ! dv/dy calculated on Y-pencil
#ifdef HALO_GLOBAL
    allocate(v2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
    allocate(wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
    k1 = ystart(3); kn = yend(3)
    i1 = ystart(1); in = yend(1)
#else
    allocate(v2(ysize(1), ysize(2), ysize(3)))
    allocate(wk2(ysize(1), ysize(2), ysize(3)))
    k1 = 1; kn = ysize(3)
    i1 = 1; in = ysize(1)
#endif
    j1 = 2; jn = ysize(2) - 1

    call transpose_x_to_y(v1,v2)
    call transpose_x_to_y(div1,wk2)

    do k=k1,kn
       do j=j1,jn
          do i=i1,in
             wk2(i,j,k) = wk2(i,j,k) + v2(i,j+1,k)-v2(i,j-1,k)
          end do
       end do
    end do

    ! dw/dz calculated on Z-pencil
#ifdef HALO_GLOBAL
    allocate(w2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
    allocate(w3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
    allocate(wk3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
    j1 = zstart(2); jn = zend(2)
    i1 = zstart(1); in = zend(1)
#else
    allocate(w2(ysize(1), ysize(2), ysize(3)))
    allocate(w3(zsize(1), zsize(2), zsize(3)))
    allocate(wk3(zsize(1), zsize(2), zsize(3)))
    j1 = 1; jn = zsize(2)
    i1 = 1; in = zsize(1)
#endif
    k1 = 2; kn = zsize(3) - 1
    
    call transpose_x_to_y(w1,w2)
    call transpose_y_to_z(w2,w3)
    call transpose_y_to_z(wk2,wk3)

    do k=k1,kn
       do j=j1,jn
          do i=i1,in
             wk3(i,j,k) = wk3(i,j,k) + w3(i,j,k+1)-w3(i,j,k-1)
          end do
       end do
    end do

    ! result in X-pencil
    call transpose_z_to_y(wk3,wk2)
    call transpose_y_to_x(wk2,div1)

    if (nrank==0) then
       write(*,*) 'Calculated via global transposition'
#ifdef DEBUG
       write(*,*) (div1(i,i,i), i=2,13)
#endif
    end if

    deallocate(v2,w2,w3,wk2,wk3)

  end subroutine test_div_transpose

  !=====================================================================
  ! Calculate divergence using halo-cell exchange (data in X-pencil)
  !=====================================================================
  subroutine test_div_haloX()

    integer :: i1, in ! I loop start/end
    integer :: j1, jn ! J loop start/end
    integer :: k1, kn ! K loop start/end
    
    ! Expected sizes
    nx_expected = nx
    ny_expected = xsize(2) + 2
    nz_expected = xsize(3) + 2
    
#ifdef HALO_GLOBAL
    allocate(div2(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    call update_halo(v1,vh,1,opt_global=.true.,opt_pencil=1)
    call update_halo(w1,wh,1,opt_global=.true.,opt_pencil=1)
    
    k1 = xstart(3); kn = xend(3)
    j1 = xstart(2); jn = xend(2)
#else
    allocate(div2(xsize(1), xsize(2), xsize(3)))
    call update_halo(v1,vh,1,opt_pencil=1)
    call update_halo(w1,wh,1,opt_pencil=1)

    k1 = 1; kn = xsize(3)
    j1 = 1; jn = xsize(2)
#endif
    i1 = 2; in = xsize(1) - 1

    call test_halo_size(vh, nx_expected, ny_expected, nz_expected, "X:v")
    call test_halo_size(wh, nx_expected, ny_expected, nz_expected, "X:w")
    
    div2 = 0.0_mytype
    do k=k1,kn
       do j=j1,jn
          do i=i1,in
             div2(i,j,k) = (u1(i+1,j,k)-u1(i-1,j,k)) &
                  + (vh(i,j+1,k)-vh(i,j-1,k)) &
                  + (wh(i,j,k+1)-wh(i,j,k-1))
          end do
       end do
    end do

    ! Compute error
    call check_err(div2, "X")
           
    deallocate(vh,wh)

  end subroutine test_div_haloX

  !=====================================================================
  ! Calculate divergence using halo-cell exchange (data in Y-pencil)
  !=====================================================================
  subroutine test_div_haloY()

    integer :: i1, in ! I loop start/end
    integer :: j1, jn ! J loop start/end
    integer :: k1, kn ! K loop start/end

    ! Expected sizes
    nx_expected = ysize(1) + 2
    ny_expected = ny
    nz_expected = ysize(3) + 2
  
#ifdef HALO_GLOBAL
    allocate(u2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
    allocate(v2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
    allocate(w2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
    allocate(div3(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
    allocate(wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
#else
    allocate(u2(ysize(1), ysize(2), ysize(3)))
    allocate(v2(ysize(1), ysize(2), ysize(3)))
    allocate(w2(ysize(1), ysize(2), ysize(3)))
    allocate(div3(xsize(1), xsize(2), xsize(3)))
    allocate(wk2(ysize(1), ysize(2), ysize(3)))
#endif
    call transpose_x_to_y(u1,u2)
    call transpose_x_to_y(v1,v2)
    call transpose_x_to_y(w1,w2)

    ! du/dx
#ifdef HALO_GLOBAL
    call update_halo(u2,uh,1,opt_global=.true.,opt_pencil=2)
    call update_halo(w2,wh,1,opt_global=.true.,opt_pencil=2)
    k1 = ystart(3); kn = yend(3)
    i1 = ystart(1); in = yend(1)
#else
    call update_halo(u2,uh,1,opt_pencil=2)
    call update_halo(w2,wh,1,opt_pencil=2)
    k1 = 1; kn = ysize(3)
    i1 = 1; in = ysize(1)
#endif
    j1 = 2; jn = ysize(2) - 1

    call test_halo_size(uh, nx_expected, ny_expected, nz_expected, "Y:u")
    call test_halo_size(wh, nx_expected, ny_expected, nz_expected, "Y:w")

    do k=k1,kn
       do j=j1,jn
          do i=i1,in
             wk2(i,j,k) = (uh(i+1,j,k)-uh(i-1,j,k)) &
                  + (v2(i,j+1,k)-v2(i,j-1,k)) &
                  + (wh(i,j,k+1)-wh(i,j,k-1))
          end do
       end do
    end do

    call transpose_y_to_x(wk2,div3)

    ! Compute error
    call check_err(div3, "Y")

    deallocate(uh,wh,wk2)

  end subroutine test_div_haloY

  !=====================================================================
  ! Calculate divergence using halo-cell exchange (data in Z-pencil)
  !=====================================================================
  subroutine test_div_haloZ()

  ! Expected sizes
  nx_expected = zsize(1) + 2
  ny_expected = zsize(2) + 2
  nz_expected = nz

#ifdef HALO_GLOBAL
  allocate(u3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
  allocate(v3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
  allocate(w3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
#else
  allocate(u3(zsize(1), zsize(2), zsize(3)))
  allocate(v3(zsize(1), zsize(2), zsize(3)))
  allocate(w3(zsize(1), zsize(2), zsize(3)))
#endif
  call transpose_y_to_z(u2,u3)
  call transpose_y_to_z(v2,v3)
  call transpose_y_to_z(w2,w3)

#ifdef HALO_GLOBAL
  allocate(div4(xstart(1):xend(1), xstart(2):xend(2), xstart(3):xend(3)))
  allocate(wk2(ystart(1):yend(1), ystart(2):yend(2), ystart(3):yend(3)))
  allocate(wk3(zstart(1):zend(1), zstart(2):zend(2), zstart(3):zend(3)))
#else
  allocate(div4(xsize(1), xsize(2), xsize(3)))
  allocate(wk2(ysize(1), ysize(2), ysize(3)))
  allocate(wk3(zsize(1), zsize(2), zsize(3)))
#endif

  ! du/dx
#ifdef HALO_GLOBAL
  call update_halo(u3,uh,1,opt_global=.true.,opt_pencil=3)
#else
  call update_halo(u3,uh,1,opt_pencil=3)
#endif

  call test_halo_size(uh, nx_expected, ny_expected, nz_expected, "Z:u")
  
#ifdef HALO_GLOBAL
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
#else
  do j=1,zsize(2)
     do i=1,zsize(1)
#endif
        do k=2,zsize(3)-1
           wk3(i,j,k) = uh(i+1,j,k)-uh(i-1,j,k)
        end do
     end do
  end do

  ! dv/dy
#ifdef HALO_GLOBAL
  call update_halo(v3,vh,1,opt_global=.true.,opt_pencil=3)
#else
  call update_halo(v3,vh,1,opt_pencil=3)
#endif

  call test_halo_size(vh, nx_expected, ny_expected, nz_expected, "Z:v")
  
#ifdef HALO_GLOBAL
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
#else
  do j=1,zsize(2)
     do i=1,zsize(1)
#endif
        do k=2,zsize(3)-1
           wk3(i,j,k) = wk3(i,j,k) + vh(i,j+1,k)-vh(i,j-1,k)
        end do
     end do
  end do

  ! dw/dz
#ifdef HALO_GLOBAL
  do j=zstart(2),zend(2)
     do i=zstart(1),zend(1)
#else
  do j=1,zsize(2)
     do i=1,zsize(1)
#endif
        do k=2,zsize(3)-1
           wk3(i,j,k) = wk3(i,j,k) + w3(i,j,k+1)-w3(i,j,k-1)
        end do
     end do
  end do

  call transpose_z_to_y(wk3,wk2)
  call transpose_y_to_x(wk2,div4)

  ! Compute error
  call check_err(div4, "Z")

  deallocate(uh,vh,wk2,wk3)
  end subroutine test_div_haloZ

  subroutine check_err(divh, pencil)

    real(mytype), dimension(:,:,:), intent(in) :: divh
    character(len=*), intent(in) :: pencil

    real(mytype) :: divmag
    
    err = mag(divh(2:xlast,2:ylast,2:zlast) - div1(2:xlast,2:ylast,2:zlast))
    divmag = mag(div1(2:xlast,2:ylast,2:zlast))
    if (err < epsilon(divmag) * divmag) then
       passing = .true.
    else
       passing = .false.
    end if
    all_pass = all_pass .and. passing

    if (nrank==0) then
       write(*,*) '-----------------------------------------------'
       write(*,*) 'Calculated via halo exchange (data in '//pencil//'-pencil)'
#ifdef DEBUG
       write(*,*) (div2(i,i,i), i=2,13)
#endif
       write(*,*) 'Error: ', err, '; Relative: ', err / divmag
       write(*,*) 'Pass: ', passing
    end if
  end subroutine check_err

  real(mytype) function mag(a)

    real(mytype), dimension(:,:,:), intent(in) :: a
    
    real(mytype) :: lmag, gmag

    lmag = sum(a(:,:,:)**2)
    call MPI_Allreduce(lmag, gmag, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
    if (ierror /= 0) then
       call decomp_2d_abort(__FILE__, __LINE__, ierror, &
            "halo_test::mag::MPI_Allreduce")
    endif

    mag = sqrt(gmag / (nx - 2) / (ny - 2) / (nz - 2))

  end function mag

  subroutine test_halo_size(arrh, nx_expected, ny_expected, nz_expected, tag)

    real(mytype), dimension(:,:,:), intent(in) :: arrh
    integer, intent(in) :: nx_expected, ny_expected, nz_expected
    character(len=*), intent(in) :: tag

    integer :: nx, ny, nz

    character(len=128) :: rank_lbl

    nx = size(arrh, 1)
    ny = size(arrh, 2)
    nz = size(arrh, 3)

    write(rank_lbl, "(A,I0,A)") "Rank", nrank, ":"

    if ((nx /= nx_expected) .or. &
         (ny /= ny_expected) .or. &
         (nz /= nz_expected)) then
       write(*,*) trim(rank_lbl), " ", tag, ":ERROR: halo size"
       write(*,*) trim(rank_lbl), " ", "+ Expected: ", nx_expected, " ", ny_expected, " ", nz_expected, " "
       write(*,*) trim(rank_lbl), " ", "+ Got:      ", nx, " ", ny, " ", nz, " "

       all_pass = .false.
    else
       write(*,*) trim(rank_lbl), " ", tag, ":PASS"
    end if
    
  end subroutine test_halo_size
  
end program halo_test
