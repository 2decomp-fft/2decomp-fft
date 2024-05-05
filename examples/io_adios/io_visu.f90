!!! SPDX-License-Identifier: BSD-3-Clause
!!!
!!! Example code to demonstrate visualisation of fields.

program visu

   use mpi

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io
   use decomp_2d_io_adios
   use decomp_2d_io_object_adios
   use decomp_2d_mpi
   use decomp_2d_testing

   implicit none

   integer, parameter :: nx_base = 17, ny_base = 13, nz_base = 11
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0

   real(mytype), allocatable, dimension(:, :, :) :: u1, u2, u3

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

      call MPI_INIT(ierr)
      ! To resize the domain we need to know global number of ranks
      ! This operation is also done as part of decomp_2d_init
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierr)
      resize_domain = int(nranks_tot / 4) + 1
      nx = nx_base * resize_domain
      ny = ny_base * resize_domain
      nz = nz_base * resize_domain
      ! Now we can check if user put some inputs
      call decomp_2d_testing_init(p_row, p_col, nx, ny, nz)

      call decomp_2d_init(nx, ny, nz, p_row, p_col)

      call decomp_2d_testing_log()

      call decomp_2d_io_init()
      call decomp_2d_register_var("u1.dat", 1, real_type)
      call decomp_2d_register_var("u2.dat", 2, real_type)
      call decomp_2d_register_var("u3.dat", 3, real_type)

   end subroutine init_example

   subroutine init_data()

      integer :: xen1, xen2, xen3
      integer :: yen1, yen2, yen3
      integer :: zen1, zen2, zen3

      call alloc_x(u1)
      call alloc_y(u2)
      call alloc_z(u3)
      xen1 = xsize(1)
      xen2 = xsize(2)
      xen3 = xsize(3)
      yen1 = ysize(1)
      yen2 = ysize(2)
      yen3 = ysize(3)
      zen1 = zsize(1)
      zen2 = zsize(2)
      zen3 = zsize(3)

      ! distribute the data
      !$acc data copy(u1,u2,u3)

      !$acc parallel loop default(present)
      do k = 1, xen3
         do j = 1, xen2
            do i = 1, xen1
               u1(i, j, k) = real(100 + nrank, mytype)
            end do
         end do
      end do
      !$acc end loop
      !$acc parallel loop default(present)
      do k = 1, yen3
         do j = 1, yen2
            do i = 1, yen1
               u2(i, j, k) = real(100 + nrank, mytype)
            end do
         end do
      end do
      !$acc end loop
      !$acc parallel loop default(present)
      do k = 1, zen3
         do j = 1, zen2
            do i = 1, zen1
               u3(i, j, k) = real(100 + nrank, mytype)
            end do
         end do
      end do
      !$acc end loop
      !$acc update self (u1)
      !$acc update self (u2)
      !$acc update self (u3)
      !$acc end data

   end subroutine init_data

   subroutine write_data()

      type(d2d_io_adios) :: io

      ! Standard I/O pattern - file per field
      call io%open_start(decomp_2d_write_mode)
      call decomp_2d_adios_write_var(io, u1, 'u1.dat')
      call decomp_2d_adios_write_var(io, u2, 'u2.dat')
      call decomp_2d_adios_write_var(io, u3, 'u3.dat')
      call io%end_close

   end subroutine write_data

   subroutine write_visu()
      ! This subroutine is based on the xdmf writers in Xcompact3d.
      ! Copyright (c) 2012-2022, Xcompact3d
      ! SPDX-License-Identifier: BSD 3-Clause
      use, intrinsic :: iso_fortran_env, only: real64

      integer :: ioxdmf

      character(len=:), allocatable :: fmt

      integer :: precision

      integer :: varctr
      character(len=16) :: filename
      character(len=2) :: varname

      if (nrank == 0) then
         OPEN (newunit=ioxdmf, file="./out.xdmf")

         write (ioxdmf, '(A22)') '<?xml version="1.0" ?>'
         write (ioxdmf, *) '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
         write (ioxdmf, *) '<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
         write (ioxdmf, *) '<Domain>'

         write (ioxdmf, '(A)') '    <Topology name="topo" TopologyType="3DCoRectMesh"'
         fmt = "(A, I0, A, I0, A, I0, A)"
         write (ioxdmf, fmt) '        Dimensions="', nz, " ", ny, " ", nx, '">'
         write (ioxdmf, '(A)') '    </Topology>'

         write (ioxdmf, *) '    <Geometry name="geo" Type="ORIGIN_DXDYDZ">'
         write (ioxdmf, *) '        <!-- Origin -->'
         write (ioxdmf, *) '        <DataItem Format="XML" Dimensions="3">'
         write (ioxdmf, *) '          0.0 0.0 0.0'
         write (ioxdmf, *) '        </DataItem>'
         write (ioxdmf, *) '        <!-- DxDyDz -->'
         write (ioxdmf, *) '        <DataItem Format="XML" Dimensions="3">'
         if (mytype == kind(0._real64)) then
            fmt = "(A, E24.17, A, E24.17, A, E24.17)"
         else
            fmt = "(A, E16.9, A, E16.9, A, E16.9)"
         end if
         write (ioxdmf, fmt) '        ', 1.0_mytype, " ", 1.0_mytype, " ", 1.0_mytype
         write (ioxdmf, *) '        </DataItem>'
         write (ioxdmf, *) '    </Geometry>'

         write (ioxdmf, *) '   <Grid Name="1" GridType="Uniform">'
         write (ioxdmf, *) '       <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
         write (ioxdmf, *) '       <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'

         do varctr = 1, 3
            write (varname, "(A, I0)") "u", varctr
            write (filename, '(A, I0, A)') "./out/u", varctr, ".dat"
            write (ioxdmf, *) '       <Attribute Name="'//trim(varname)//'" Center="Node">'
            write (ioxdmf, *) '          <DataItem Format="HDF"'

#ifdef DOUBLE_PREC
            print *, "Double precision build"
#ifdef SAVE_SINGLE
            precision = 4
#else
            precision = 8
#endif
#else
            precision = 4
#endif
            write (ioxdmf, "(A,I0,A)") '            DataType="Float" Precision="', precision, '" Endian="little" Seek="0"'

            fmt = "(A, I0, A, I0, A, I0, A)"
            write (ioxdmf, fmt) '            Dimensions="', nz, " ", ny, " ", nx, '">'

            write (ioxdmf, *) '              '//trim(filename)

            write (ioxdmf, *) '           </DataItem>'
            write (ioxdmf, *) '        </Attribute>'
         end do

         write (ioxdmf, '(/)')
         write (ioxdmf, *) '    </Grid>'
         write (ioxdmf, *) '</Domain>'
         write (ioxdmf, '(A7)') '</Xdmf>'
         close (ioxdmf)
      end if

   end subroutine write_visu

   subroutine fin()

      integer :: ierr

      call decomp_2d_io_fin
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)

   end subroutine fin

end program visu

