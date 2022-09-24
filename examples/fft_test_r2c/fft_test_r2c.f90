!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main test program for the FFT r2c/c2r interface
!  also demonstrate the use of the IO library 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program fft_test_r2c

  use decomp_2d
  use decomp_2d_fft
  use glassman
  ! use decomp_2d_io
  
  use MPI
  
  implicit none
  !include "fftw3.f"
  
  integer, parameter :: nx=4, ny=2, nz=3
  integer :: p_row=0, p_col=0
  
  real(mytype), allocatable, dimension(:,:,:) :: in, in2
  complex(mytype), allocatable, dimension(:,:,:) :: out
  
  integer, dimension(3) :: fft_start, fft_end, fft_size
  
  real(mytype), dimension(nx,ny,nz) :: in_global, in_g2
  complex(mytype), dimension(nx/2+1,ny,nz) :: out_global
  
  real(mytype) :: err
  !integer*8 :: plan
  integer :: ierror, i,j,k
  
  ! Init
  call MPI_INIT(ierror)
  if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_INIT")
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call decomp_2d_fft_init(PHYSICAL_IN_X)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Compute a small problem all on rank 0 as reference
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !call random_number(in_global)
  ! No Random input, but something global for everyone
  in_global(1,1,1) = (1.000)
  in_global(1,1,2) = (0.999)
  in_global(1,1,3) = (0.987)
  in_global(1,2,1) = (0.936)
  in_global(1,2,2) = (0.994)
  in_global(1,2,3) = (0.989)
  in_global(2,1,1) = (0.963)
  in_global(2,1,2) = (0.891)
  in_global(2,1,3) = (0.903)
  in_global(2,2,1) = (0.885)
  in_global(2,2,2) = (0.823)
  in_global(2,2,3) = (0.694)
  in_global(3,1,1) = (0.500)
  in_global(3,1,2) = (0.499)
  in_global(3,1,3) = (0.487)
  in_global(3,2,1) = (0.436)
  in_global(3,2,2) = (0.494)
  in_global(3,2,3) = (0.489)
  in_global(4,1,1) = (0.463)
  in_global(4,1,2) = (0.391)
  in_global(4,1,3) = (0.403)
  in_global(4,2,1) = (0.385)
  in_global(4,2,2) = (0.323)
  in_global(4,2,3) = (0.194)
  
  if (nrank==0) then
#ifdef DEBUG
     write(*,*) '*** Reference serial computation on rank 0 only'
     write(*,*) ' global real input'
     do i=1,nx
        write(*,20) ((in_global(i,j,k),j=1,ny),k=1,nz)
     end do
#endif
     
     ! Using a 3D FFT routine supplied by this library
     call glassman_3d_r2c(in_global,nx,ny,nz,out_global)
     
     ! If using FFTW library:
     !  - uncomment the FFTW include file & plan above
     !  - uncomment the follwing function calls
     !  - change names to dfftw... for double precision 
     !call sfftw_plan_dft_r2c_3d(plan,nx,ny,nz, &
     !     in_global,out_global,FFTW_ESTIMATE)
     !call sfftw_execute_dft_r2c(plan,in_global,out_global)
     
#ifdef DEBUG
     write(*,*) ' global complex output'
     do i=1,nx/2+1
        write(*,10) ((out_global(i,j,k),j=1,ny),k=1,nz)
     end do
#endif

  end if
10 format(1x,6(:,'(',F5.2,',',F5.2,')'))
20 format(1x,6F5.2)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c) 
  
  !  input is X-pencil real    data whose global size is nx*ny*nz
  ! output is Z-pencil complex data whose global size is (nx/2+1)*ny*nz
  allocate (in(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! each processor gets its local portion of global data
  do k=xstart(3),xend(3)
     do j=xstart(2),xend(2)
        do i=xstart(1),xend(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do
  
#ifdef DEBUG
  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (X-pencil input)'
     write(*,*) ' real input held by rank 0: ', nrank
     write(*,20) in
  end if
#endif
  
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)
  
#ifdef DEBUG
  if (nrank==0) then
     write(*,*) ' - after forward transform'
     write(*,*) ' complex output held by rank 0: ', nrank
     write(*,10) out
  end if
#endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the complex-to-real interface (c2r) 

  allocate (in2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
  
  ! compute c2r transform
  call decomp_2d_fft_3d(out,in2)
  
  ! normalisation
  in2 = in2 / real(nx, kind=mytype) / real(ny, kind=mytype) / real(nz, kind=mytype)

#ifdef DEBUG
  if (nrank==0) then
     write(*,*) ' - after backward transform and normalisation'
     write(*,*) ' real output held by rank 0: ', nrank
     write(*,20) in2
  end if
#endif

  ! Create the global array to verify the output
  call assemble_global(1,in2,in_g2,nx,ny,nz)

  deallocate(in,in2,out)
  call decomp_2d_fft_finalize

  ! check on rank 0 if input data is properly recovered
  if (nrank==0) then
     err = 0._mytype
     do k=1,nz
        do j=1,ny
           do i=1,nx
              err = err + (in_g2(i,j,k)-in_global(i,j,k))**2
           end do
        end do
     end do
     err = err / real(nx,mytype) / real(ny,mytype) / real(nz,mytype)
#ifdef DEBUG
     write(*,*)' Output on master node of the full result'
     write(*,*) ' error / mesh point: ', sqrt(err)
     write(*,*) '*** Output computation on rank 0 only'
     write(*,*) ' global real result'
     do i=1,nx
        write(*,20) ((in_g2(i,j,k),j=1,ny),k=1,nz)
     end do
#endif

     ! Abort in case of error
     if (maxval(abs(in_g2 - in_global)) > 10*epsilon(err)) &
        call decomp_2d_abort(1, "error in fft_r2c")

  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Repeat the above but using Z-pencil input
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  call decomp_2d_fft_init(PHYSICAL_IN_Z)
  ! reset the global output
  in_g2 = 0._mytype
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the real-to-complex interface (r2c) 
  
  allocate (in(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
  call decomp_2d_fft_get_size(fft_start,fft_end,fft_size)
  allocate (out(fft_start(1):fft_end(1), &
       fft_start(2):fft_end(2), &
       fft_start(3):fft_end(3)))
  
  ! each processor gets its local portion of global data
  do k=zstart(3),zend(3)
     do j=zstart(2),zend(2)
        do i=zstart(1),zend(1)
           in(i,j,k) = in_global(i,j,k)
        end do
     end do
  end do

#ifdef DEBUG
  if (nrank==0) then
     write(*,*) ' '
     write(*,*) '*** Distributed computation (Z-pencil input)'
     write(*,*) ' real input held by rank 0:'
     write(*,20) in
  end if
#endif
  
  ! compute r2c transform 
  call decomp_2d_fft_3d(in,out)

#ifdef DEBUG
  if (nrank==0) then
     write(*,*) ' - after forward transform'
     write(*,*) ' complex output held by rank 0:'
     write(*,10) out
  end if
#endif
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Test the complex-to-real interface (c2r) 
  
  allocate (in2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
  
  ! compute c2r transform
  call decomp_2d_fft_3d(out,in2)
  
  ! normalisation
  in2 = in2 / real(nx) / real(ny) / real(nz)

#ifdef DEBUG
  if (nrank==0) then
     write(*,*) ' - after backward transform and normalisation'
     write(*,*) ' real output held by rank 0:'
     write(*,20) in2
  end if
#endif

  ! Create the global array to verify the output
  call assemble_global(3,in2,in_g2,nx,ny,nz)
  ! check on rank 0 if input data is properly recovered
  ! this also tests the IO routines
  if (nrank==0) then
     err = 0._mytype
     do k=1,nz
        do j=1,ny
           do i=1,nx
              err = err + (in_g2(i,j,k)-in_global(i,j,k))**2
           end do
        end do
     end do
     err = err / real(nx,mytype) / real(ny,mytype) / real(nz,mytype)
#ifdef DEBUG
     write(*,*)' Output on master node of the full result'
     write(*,*) ' error / mesh point: ', sqrt(err)
     write(*,*) '*** Output computation on rank 0 only'
     write(*,*) ' global real result'
     do i=1,nx
        write(*,20) ((in_g2(i,j,k),j=1,ny),k=1,nz)
     end do
#endif

     ! Abort in case of error
     if (maxval(abs(in_g2 - in_global)) > 10*epsilon(err)) &
        call decomp_2d_abort(2, "error in fft_r2c")

  end if

  deallocate(in,in2,out)
  
  if (nrank == 0) then
     write(*,*) " "
     write(*,*) " fft_test_r2c completed"
     write(*,*) " "
  endif

  call decomp_2d_fft_finalize
  call decomp_2d_finalize
  call MPI_FINALIZE(ierror)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Collect data from each processor and assemble into a global array
! at the master rank
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine assemble_global(ndir,local,global,nx,ny,nz)
  
  use decomp_2d
  use MPI
  
  implicit none
  
  integer, intent(IN) :: ndir  ! 1 = X-pencil; 3 = Z-pencil
  integer, intent(IN) :: nx,ny,nz
  real(mytype), dimension(:,:,:), intent(IN) :: local
  real(mytype), dimension(nx,ny,nz), intent(OUT) :: global
  
  real(mytype), allocatable, dimension(:,:,:) :: rbuf
  integer, dimension(9) :: sbuf1, rbuf1
  
  integer :: ierror, i,j,k,m, i1,i2,j1,j2,k1,k2, count
  integer, dimension(MPI_STATUS_SIZE) :: status

  do k=1,nz
     do j=1,ny
        do i=1,nx
         global(i,j,k) = 0._mytype
        end do
     end do
  end do
  if (nrank==0) then
     ! master writes its own data to a global array
     if (ndir==3) then  ! Z-pencil
        i1 = zstart(1)
        i2 = zend(1)
        j1 = zstart(2)
        j2 = zend(2)
        k1 = zstart(3)
        k2 = zend(3)
     else if (ndir==1) then  ! X-pencil
        i1 = xstart(1)
        i2 = xend(1)
        j1 = xstart(2)
        j2 = xend(2)
        k1 = xstart(3)
        k2 = xend(3)
     end if
     do k=k1,k2
        do j=j1,j2
           do i=i1,i2
              ! 'local' is assumbed shape array
              ! but it is OK as starting index for rank 0 always 1
              global(i,j,k)=local(i,j,k)
           end do
        end do
     end do
     ! then loop through all other ranks to collect data
     do m=1,nproc-1
        CALL MPI_RECV(rbuf1,9,MPI_INTEGER,m,m,MPI_COMM_WORLD, &
             status,ierror)
        if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_RECV")
        allocate(rbuf(rbuf1(1):rbuf1(2),rbuf1(4):rbuf1(5), &
             rbuf1(7):rbuf1(8)))
        CALL MPI_RECV(rbuf,rbuf1(3)*rbuf1(6)*rbuf1(9),real_type,m, &
             m+nproc,MPI_COMM_WORLD,status,ierror)
        if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_RECV")
        do k=rbuf1(7),rbuf1(8)
           do j=rbuf1(4),rbuf1(5)
              do i=rbuf1(1),rbuf1(2)
                 global(i,j,k)=rbuf(i,j,k)
              end do
           end do
        end do
        deallocate(rbuf)
     end do
  else
     ! slaves send data to mater
     if (ndir==3) then  ! Z-pencil
        sbuf1(1) = zstart(1)
        sbuf1(2) = zend(1)
        sbuf1(3) = zsize(1)
        sbuf1(4) = zstart(2)
        sbuf1(5) = zend(2)
        sbuf1(6) = zsize(2)
        sbuf1(7) = zstart(3)
        sbuf1(8) = zend(3)
        sbuf1(9) = zsize(3)
        count = zsize(1)*zsize(2)*zsize(3)
     else if (ndir==1) then  ! X-pencil
        sbuf1(1) = xstart(1)
        sbuf1(2) = xend(1)
        sbuf1(3) = xsize(1)
        sbuf1(4) = xstart(2)
        sbuf1(5) = xend(2)
        sbuf1(6) = xsize(2)
        sbuf1(7) = xstart(3)
        sbuf1(8) = xend(3)
        sbuf1(9) = xsize(3)
        count = xsize(1)*xsize(2)*xsize(3)
     end if
     ! send partition information
     CALL MPI_SEND(sbuf1,9,MPI_INTEGER,0,nrank,MPI_COMM_WORLD,ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_SEND")
     ! send data array
     CALL MPI_SEND(local,count,real_type,0, &
          nrank+nproc,MPI_COMM_WORLD,ierror)
     if (ierror /= 0) call decomp_2d_abort(__FILE__, __LINE__, ierror, "MPI_SEND")
  end if
  
  return
end subroutine assemble_global

end program fft_test_r2c
