!!! SPDX-License-Identifier: BSD-3-Clause
!!!
!!! Example code to demonstrate visualisation of fields.

program visu

  use mpi

  use decomp_2d
  use decomp_2d_mpi
  use decomp_2d_constants
  use decomp_2d_io

  implicit none

  integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
  integer :: nx, ny, nz
  integer :: p_row = 0, p_col = 0

#ifdef COMPLEX_TEST
  complex(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
#else
  real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3
#endif

  integer, parameter :: output2D = 0 ! Which plane to write in 2D (0 for 3D)
  character(len=*), parameter :: io_name = "visu-io"
#ifndef ADIOS2
  logical :: dir_exists
#endif

  integer :: nargin, arg, FNLength, status, DecInd

  integer :: i, j, k

  integer :: ierr

  call init_example()
  call init_data()

  call write_data()
  call write_visu()
  
  call fin()
  
contains

  subroutine init_example()

    integer :: resize_domain
    integer :: nranks_tot
    character(len=80) :: InputFN

    call MPI_INIT(ierr)
    ! To resize the domain we need to know global number of ranks
    ! This operation is also done as part of decomp_2d_init
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierr)
    resize_domain = int(nranks_tot / 4) + 1
    nx = nx_base * resize_domain
    ny = ny_base * resize_domain
    nz = nz_base * resize_domain
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
       call MPI_COMM_RANK(MPI_COMM_WORLD, nrank, ierr)
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

    call decomp_2d_io_init()
    call decomp_2d_init_io(io_name)
#ifndef COMPLEX_TEST
    call decomp_2d_register_variable(io_name, "u1.dat", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "u2.dat", 2, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "u3.dat", 3, 0, output2D, mytype)
#else
    call decomp_2d_register_variable(io_name, "u1_re.dat", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "u2_re.dat", 2, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "u3_re.dat", 3, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "u1_im.dat", 1, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "u2_im.dat", 2, 0, output2D, mytype)
    call decomp_2d_register_variable(io_name, "u3_im.dat", 3, 0, output2D, mytype)
#endif

  end subroutine init_example

  subroutine init_data()

    call alloc_x(u1)
    call alloc_y(u2)
    call alloc_z(u3)
    
    do k = 1, xsize(3)
       do j = 1, xsize(2)
          do i = 1, xsize(1)
#ifndef COMPLEX_TEST
             u1(i, j, k) = real(nrank, mytype)
#else
             u1(i, j, k) = cmplx(real(nrank, mytype), -real(nrank, mytype))
#endif
          end do
       end do
    end do
    do k = 1, ysize(3)
       do j = 1, ysize(2)
          do i = 1, ysize(1)
#ifndef COMPLEX_TEST
             u2(i, j, k) = real(nrank, mytype)
#else
             u2(i, j, k) = cmplx(real(nrank, mytype), -real(nrank, mytype))
#endif
          end do
       end do
    end do
    do k = 1, zsize(3)
       do j = 1, zsize(2)
          do i = 1, zsize(1)
#ifndef COMPLEX_TEST
             u3(i, j, k) = real(nrank, mytype)
#else
             u3(i, j, k) = cmplx(real(nrank, mytype), -real(nrank, mytype))
#endif
          end do
       end do
    end do
    
  end subroutine init_data

  subroutine write_data()
    
    !! Write arrays + visu data from orientations 1, 2 and 3
#ifndef ADIOS2
    if (nrank == 0) then
       inquire (file="out", exist=dir_exists)
       if (.not. dir_exists) then
          call execute_command_line("mkdir out 2> /dev/null")
       end if
    end if
#endif
    
    ! Standard I/O pattern - file per field
#ifdef ADIOS2
    call decomp_2d_open_io(io_name, "out", decomp_2d_write_mode)
    call decomp_2d_start_io(io_name, "out")
#endif

#ifndef COMPLEX_TEST
    call decomp_2d_write_one(1, u1, 'out', 'u1.dat', 0, io_name)
    call decomp_2d_write_one(2, u2, 'out', 'u2.dat', 0, io_name)
    call decomp_2d_write_one(3, u3, 'out', 'u3.dat', 0, io_name)
#else
    call decomp_2d_write_one(1, real(u1, mytype), 'out', 'u1_re.dat', 0, io_name)
    call decomp_2d_write_one(2, real(u1, mytype), 'out', 'u2_re.dat', 0, io_name)
    call decomp_2d_write_one(3, real(u1, mytype), 'out', 'u3_re.dat', 0, io_name)
    call decomp_2d_write_one(1, aimag(u1, mytype), 'out', 'u1_im.dat', 0, io_name)
    call decomp_2d_write_one(2, aimag(u2, mytype), 'out', 'u2_im.dat', 0, io_name)
    call decomp_2d_write_one(3, aimag(u3, mytype), 'out', 'u3_im.dat', 0, io_name)
#endif

#ifdef ADIOS2
    call decomp_2d_end_io(io_name, "out")
    call decomp_2d_close_io(io_name, "out")
#endif

  end subroutine write_data

  subroutine write_visu()
    ! This subroutine is based on the xdmf writers in Xcompact3d.
    ! Copyright (c) 2012-2022, Xcompact3d
    ! SPDX-License-Identifier: BSD 3-Clause

    integer :: ioxdmf

    character(len=:), allocatable :: fmt

    integer :: precision

    integer :: varctr
    character(len=16) :: filename
    character(len=2) :: varname
    
    if (nrank == 0) then
       OPEN(newunit=ioxdmf, file="./out.xdmf")

       write(ioxdmf,'(A22)')'<?xml version="1.0" ?>'
       write(ioxdmf,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
       write(ioxdmf,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
       write(ioxdmf,*)'<Domain>'

       write(ioxdmf,'(A)')'    <Topology name="topo" TopologyType="3DCoRectMesh"'
       fmt = "(A, I0, A, I0, A, I0, A)"
       write(ioxdmf,fmt)'        Dimensions="', nz, " ", ny, " ", nx, '">'
       write(ioxdmf,'(A)')'    </Topology>'

       write(ioxdmf,*)'    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
       write(ioxdmf,*)'        <!-- Origin -->'
       write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
       write(ioxdmf,*)'          0.0 0.0 0.0'
       write(ioxdmf,*)'        </DataItem>'
       write(ioxdmf,*)'        <!-- DxDyDz -->'
       write(ioxdmf,*)'        <DataItem Format="XML" Dimensions="3">'
       if (mytype == kind(0.0d0)) then
          fmt = "(A, E24.17, A, E24.17, A, E24.17)"
       else
          fmt = "(A, E16.9, A, E16.9, A, E16.9)"
       end if
       write(ioxdmf,fmt)'        ', 1.0_mytype, " ", 1.0_mytype, " ", 1.0_mytype
       write(ioxdmf,*)'        </DataItem>'
       write(ioxdmf,*)'    </Geometry>'

       write(ioxdmf, *)'   <Grid Name="1" GridType="Uniform">'
       write(ioxdmf, *)'       <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
       write(ioxdmf, *)'       <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

       do varctr = 1, 3
#ifndef COMPLEX_TEST
          write(varname, "(A, I0)") "u", varctr
          write(filename, '(A, I0, A)') "./out/u", varctr, ".dat"
#else
          write(varname, "(A, I0, A)") "u", varctr, "_re"
          write(filename, '(A, I0, A)') "./out/u", varctr, "_re.dat"
          write(varname, "(A, I0, A)") "u", varctr, "_im"
          write(filename, '(A, I0, A)') "./out/u", varctr, "_im.dat"
#endif
          write(ioxdmf,*)'       <Attribute Name="'//trim(varname)//'" Center="Node">'
#ifndef ADIOS2
          write(ioxdmf,*)'          <DataItem Format="Binary"'
#else
          write(ioxdmf,*)'          <DataItem Format="HDF"'
#endif

#ifdef DOUBLE_PREC
          print *, "Double precision build"
#ifdef SAVE_SINGLE
          if (output2D.eq.0) then
             precision = 4
          else
             precision = 8
          endif
#else
          precision = 8
#endif
#else
          precision = 4
#endif
          write(ioxdmf,"(A,I0,A)")'            DataType="Float" Precision="', precision, '" Endian="little" Seek="0"'

          fmt = "(A, I0, A, I0, A, I0, A)"
          write(ioxdmf,fmt)'            Dimensions="', nz, " ", ny, " ", nx, '">'

          write(ioxdmf,*)'              '//trim(filename)

          write(ioxdmf,*)'           </DataItem>'
          write(ioxdmf,*)'        </Attribute>'
       end do
       
       write(ioxdmf,'(/)')
       write(ioxdmf,*)'    </Grid>'
       write(ioxdmf,*)'</Domain>'
       write(ioxdmf,'(A7)')'</Xdmf>'
       close(ioxdmf)
    endif
    
  end subroutine write_visu
  
  subroutine fin()

    integer :: ierr
    
    call decomp_2d_finalize
    call MPI_FINALIZE(ierr)

 end subroutine fin
  
end program visu

