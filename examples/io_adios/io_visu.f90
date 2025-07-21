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
      call decomp_2d_register_var("u1", 1, real_type)
      call decomp_2d_register_var("u2", 2, real_type)
      call decomp_2d_register_var("u3", 3, real_type)

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
      call decomp_2d_adios_write_var(io, u1, 'u1')
      call decomp_2d_adios_write_var(io, u2, 'u2')
      call decomp_2d_adios_write_var(io, u3, 'u3')
      call io%end_close

   end subroutine write_data

   subroutine write_visu()
      ! This subroutine is based on the xdmf writers in Xcompact3d.
      ! Copyright (c) 2012-2022, Xcompact3d
      ! SPDX-License-Identifier: BSD 3-Clause
      use, intrinsic :: iso_fortran_env, only: real64

      type(d2d_io_family), pointer :: family
      integer :: ioxml
      real(real64) :: lx, ly, lz, dx, dy, dz

      ! Size of the domain
      lx = 3.
      ly = 2.
      lz = 1.

      ! Uniform grid
      dx = lx / (nx - 1)
      dy = ly / (ny - 1)
      dz = lz / (nz - 1)

      if (nrank == 0) then

         ! Open the file
         family => decomp_2d_adios_get_default_family()
         open (newunit=ioxml, file=family%get_folder()//"/vtk.xml")
         nullify (family)

         ! Header for a uniform grid
         write (ioxml, *) '<?xml version="1.0"?>'
         write (ioxml, *) '<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">'
         ! Size of the domain : [3, 2, 1]
         ! Extent should be in reversed order
         write (ioxml, *) '  <ImageData WholeExtent="1 ', nz, ' 1 ', ny, ' 1 ', nx, &
            '" Origin="0 0 0" Spacing="', dx, ' ', dy, ' ', dz, '">'
         write (ioxml, *) '    <Piece Extent="1 ', nz, ' 1 ', ny, ' 1 ', nx, '">'

         ! Data
         write (ioxml, *) '      <PointData>'
         write (ioxml, *) '        <DataArray Name="u1" />'
         write (ioxml, *) '        <DataArray Name="u2" />'
         write (ioxml, *) '        <DataArray Name="u3" />'
         write (ioxml, *) '      </PointData>'

         ! Footer
         write (ioxml, *) '    </Piece>'
         write (ioxml, *) '  </ImageData>'
         write (ioxml, *) '</VTKFile>'

      end if

   end subroutine write_visu

   subroutine fin()

      integer :: ierr

      call decomp_2d_io_fin
      call decomp_2d_finalize
      call MPI_FINALIZE(ierr)

   end subroutine fin

end program visu

