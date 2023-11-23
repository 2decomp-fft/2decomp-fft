!! SPDX-License-Identifier: BSD-3-Clause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This example calculates the gradient of a periodic field using global
! transposition
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program grad3d

   use mpi

   use decomp_2d
   use decomp_2d_constants
   use decomp_2d_io
   use decomp_2d_mpi
   use decomp_2d_testing
   use MPI
#if defined(_GPU)
   use cudafor
   use openacc
#endif

   implicit none

   integer, parameter :: nx_base = 64, ny_base = 64, nz_base = 64
   integer :: nx, ny, nz
   integer :: p_row = 0, p_col = 0
   integer :: resize_domain
   integer :: nranks_tot
   integer :: ierror
   logical :: all_pass

   real(mytype), parameter :: lx = 1.0_mytype
   real(mytype), parameter :: ly = 1.0_mytype
   real(mytype), parameter :: lz = 1.0_mytype

   real(mytype) :: dx, dy, dz
   real(mytype) :: error_ref

   real(mytype), allocatable, dimension(:, :, :) :: phi1, phi2, phi3
   real(mytype), allocatable, dimension(:, :, :) :: dphiX, dphiY, dphiz
   real(mytype), allocatable, dimension(:, :, :) :: wk2, wk3

   integer, parameter :: data_type = d2d_real

   call MPI_INIT(ierror)
   ! To resize the domain we need to know global number of ranks
   ! This operation is also done as part of decomp_2d_init
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nranks_tot, ierror)
   resize_domain = int(nranks_tot / 4) + 1
   nx = nx_base * resize_domain
   ny = ny_base * resize_domain
   nz = nz_base * resize_domain
   ! Now we can check if user put some inputs
   call decomp_2d_testing_init(p_row, p_col, nx, ny, nz)

   call decomp_2d_init(nx, ny, nz, p_row, p_col)

   call decomp_2d_testing_log()

   dx = lx / real(nx, mytype)
   dy = ly / real(ny, mytype)
   dz = lz / real(nz, mytype)

   if (nrank == 0) then
      write (*, *) '-----------------------------------------------'
      write (*, *) "Mesh Resolution ", nx, ny, nz
      write (*, *) '-----------------------------------------------'
   end if

   call allocate_var()
   !$acc data create(phi2,phi3,wk2,wk3) copy(phi1,dphiX, dphiY, dphiZ)
   call init_phi()
   call compute_grad()
   call test_derX(dphiX)
   call test_derY(dphiY)
   call test_derZ(dphiZ)
   call write_data()
   !$acc end data

   if (nrank == 0) then
      write (*, *) '-----------------------------------------------'
      write (*, *) "End GRAD calculation check all pass: ", all_pass
      write (*, *) '==============================================='
   end if

   call finalize()

   call decomp_2d_finalize

   if (.not. all_pass) call decomp_2d_abort(1, "Error in grad3d")

   call MPI_FINALIZE(ierror)

contains

   !=====================================================================
   ! Initialize
   !=====================================================================
   subroutine allocate_var()

      use decomp_2d

      implicit none

      logical, parameter :: glob = .false.

      ! Allocate main variables in X-pencil
      call alloc_x(phi1, glob)
      call alloc_x(dphiX, glob)
      call alloc_x(dphiY, glob)
      call alloc_x(dphiZ, glob)

      ! Working array used more than once
      call alloc_y(phi2, glob)
      call alloc_y(wk2, glob)

      call alloc_z(phi3, glob)
      call alloc_z(wk3, glob)

      all_pass = .true.

   end subroutine allocate_var

   !=====================================================================
   ! Initialize the scalar field
   !=====================================================================
   subroutine init_phi()

      implicit none

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Scalar field
      !$acc kernels default(present)
      do concurrent(k=1:xe3, j=1:xe2, i=1:xe1)
         z = (k + xs3 - 2) * dz
         y = (j + xs2 - 2) * dy
         x = (i + xs1 - 2) * dx
         phi1(i, j, k) = -2._mytype * cos(twopi * (x / lx)) * cos(twopi * (y / ly)) * sin(twopi * (z / lz))
      end do
      !$acc end kernels

   end subroutine init_phi
   !=====================================================================
   ! Finalize with deallocation of arrays
   !=====================================================================
   subroutine finalize()

      implicit none

      deallocate (phi1, phi2, phi3)
      deallocate (dphiX, dphiY, dphiZ)
      deallocate (wk2, wk3)

   end subroutine finalize
   !=====================================================================
   ! Calculate gradient using global transposition
   !=====================================================================
   subroutine compute_grad()

      implicit none

      ! Compute X derivative
      call derx(dphiX, phi1, dx, xsize(1), xsize(2), xsize(3))

      ! Compute Y derivative
      call transpose_x_to_y(phi1, phi2)
      call dery(wk2, phi2, dy, ysize(1), ysize(2), ysize(3))
      call transpose_y_to_x(wk2, dphiY)

      ! Compute Z derivative
      call transpose_y_to_z(phi2, phi3)
      call derz(wk3, phi3, dz, zsize(1), zsize(2), zsize(3))
      call transpose_z_to_y(wk3, wk2)
      call transpose_y_to_x(wk2, dphiZ)

   end subroutine compute_grad
   !=====================================================================
   ! Calculate gradient in X (data in X-pencil)
   !=====================================================================
   subroutine derx(df, ff, delta, nx, ny, nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx, ny, nz) :: df
      real(mytype), intent(in), dimension(nx, ny, nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff / delta

      !$acc kernels default(present)
      do concurrent(k=1:nz, j=1:ny)
         df(1, j, k) = coeff * (ff(2, j, k) - ff(nx, j, k))
         do concurrent(i=2:nx - 1)
            df(i, j, k) = coeff * (ff(i + 1, j, k) - ff(i - 1, j, k))
         end do
         df(nx, j, k) = coeff * (ff(1, j, k) - ff(nx - 1, j, k))
      end do
      !$acc end kernels

   end subroutine derx
   !=====================================================================
   ! Calculate gradient in Y (data in Y-pencil)
   !=====================================================================
   subroutine dery(df, ff, delta, nx, ny, nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx, ny, nz) :: df
      real(mytype), intent(in), dimension(nx, ny, nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff / delta

      !$acc kernels default(present)
      do concurrent(k=1:nz)
         do concurrent(i=1:nx)
            df(i, 1, k) = coeff * (ff(i, 2, k) - ff(i, ny, k))
         end do
         do concurrent(j=2:ny - 1, i=1:nx)
            df(i, j, k) = coeff * (ff(i, j + 1, k) - ff(i, j - 1, k))
         end do
         do concurrent(i=1:nx)
            df(i, ny, k) = coeff * (ff(i, 1, k) - ff(i, ny - 1, k))
         end do
      end do
      !$acc end kernels

   end subroutine dery
   !=====================================================================
   ! Calculate gradient in Z (data in Z-pencil)
   !=====================================================================
   subroutine derz(df, ff, delta, nx, ny, nz)

      implicit none

      ! Arguments
      integer, intent(in) :: nx, ny, nz
      real(mytype), intent(in) :: delta
      real(mytype), intent(out), dimension(nx, ny, nz) :: df
      real(mytype), intent(in), dimension(nx, ny, nz) :: ff

      ! Local variables
      integer :: i, j, k
      real(mytype) :: coeff = 0.5_mytype

      coeff = coeff / delta

      !$acc kernels default(present)
      do concurrent(j=1:ny, i=1:nx)
         df(i, j, 1) = coeff * (ff(i, j, 2) - ff(i, j, nz))
      end do
      do concurrent(k=2:nz - 1, j=1:ny, i=1:nx)
         df(i, j, k) = coeff * (ff(i, j, k + 1) - ff(i, j, k - 1))
      end do
      do concurrent(j=1:ny, i=1:nx)
         df(i, j, nz) = coeff * (ff(i, j, 1) - ff(i, j, nz - 1))
      end do
      !$acc end kernels

   end subroutine derz
   !=====================================================================
   ! Test derivatives against analytical solution (data in X-pencil)
   !=====================================================================
   subroutine test_derX(df)

      implicit none
      ! Arguments
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: df

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      real(mytype) :: dphi, dphi_num
      real(mytype) :: error = 0._mytype
      real(mytype) :: err_all = 0._mytype
      real(mytype) :: dphi2 = 0._mytype
      real(mytype) :: sum_dphi2 = 0._mytype
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Compute the error against analytical solution
      !$acc parallel loop default(present) reduction(+:error)
      do k = 1, xe3
         do j = 1, xe2
            do i = 1, xe1
               z = (k + xs3 - 2) * dz
               y = (j + xs2 - 2) * dy
               x = (i + xs1 - 2) * dx
               dphi = -2._mytype * (twopi / lx) &
                      * sin(twopi * (x / lx)) * cos(twopi * (y / ly)) * sin(twopi * (z / lz))
               dphi2 = dphi2 + dphi * dphi
               dphi_num = df(i, j, k)
               error = error + (dphi - dphi_num) * (dphi - dphi_num)
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(dphi2, sum_dphi2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      err_all = sqrt(err_all / sum_dphi2) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))
      error_ref = err_all

      if (nrank == 0) then
         write (*, *) 'DX error / mesh point: ', err_all
      end if

   end subroutine test_derX

   !=====================================================================
   ! Test derivatives against analytical solution (data in Y-pencil)
   !=====================================================================
   subroutine test_derY(df)

      implicit none
      ! Arguments
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: df

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      real(mytype) :: dphi, dphi_num
      real(mytype) :: error = 0._mytype
      real(mytype) :: err_all = 0._mytype
      real(mytype) :: dphi2 = 0._mytype
      real(mytype) :: sum_dphi2 = 0._mytype
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Compute the error against analytical solution
      !$acc parallel loop default(present) reduction(+:error)
      do k = 1, xe3
         do j = 1, xe2
            do i = 1, xe1
               z = (k + xs3 - 2) * dz
               y = (j + xs2 - 2) * dy
               x = (i + xs1 - 2) * dx
               dphi = -2._mytype * (twopi / ly) &
                      * cos(twopi * (x / lx)) * sin(twopi * (y / ly)) * sin(twopi * (z / lz))
               dphi2 = dphi2 + dphi * dphi
               dphi_num = df(i, j, k)
               error = error + (dphi - dphi_num) * (dphi - dphi_num)
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(dphi2, sum_dphi2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      err_all = sqrt(err_all / sum_dphi2) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

      if (nrank == 0) then
         write (*, *) 'DY error / mesh point: ', err_all
      end if

      if (abs(err_all - error_ref) > 1.0e-5_mytype) all_pass = .false.

   end subroutine test_derY

   !=====================================================================
   ! Test derivatives against analytical solution (data in Z-pencil)
   !=====================================================================
   subroutine test_derZ(df)

      implicit none
      ! Arguments
      real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: df

      integer :: i, j, k
      real(mytype), parameter :: twopi = 2._mytype * acos(-1._mytype)
      real(mytype) :: x, y, z
      real(mytype) :: dphi, dphi_num
      real(mytype) :: error = 0._mytype
      real(mytype) :: err_all = 0._mytype
      real(mytype) :: dphi2 = 0._mytype
      real(mytype) :: sum_dphi2 = 0._mytype
      integer :: xe1, xe2, xe3
      integer :: xs1, xs2, xs3

      xe1 = xsize(1)
      xe2 = xsize(2)
      xe3 = xsize(3)
      xs1 = xstart(1)
      xs2 = xstart(2)
      xs3 = xstart(3)

      ! Compute the error against analytical solution
      !$acc parallel loop default(present) reduction(+:error)
      do k = 1, xe3
         do j = 1, xe2
            do i = 1, xe1
               z = (k + xs3 - 2) * dz
               y = (j + xs2 - 2) * dy
               x = (i + xs1 - 2) * dx
               dphi = 2._mytype * (twopi / lz) &
                      * cos(twopi * (x / lx)) * cos(twopi * (y / ly)) * cos(twopi * (z / lz))
               dphi2 = dphi2 + dphi * dphi
               dphi_num = df(i, j, k)
               error = error + (dphi - dphi_num) * (dphi - dphi_num)
            end do
         end do
      end do
      !$acc end loop
      call MPI_ALLREDUCE(error, err_all, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      call MPI_ALLREDUCE(dphi2, sum_dphi2, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      err_all = sqrt(err_all / sum_dphi2) / (real(nx, mytype) * real(ny, mytype) * real(nz, mytype))

      if (nrank == 0) then
         write (*, *) 'DZ error / mesh point: ', err_all
      end if

      if (abs(err_all - error_ref) > 1.0e-5_mytype) all_pass = .false.

   end subroutine test_derZ
   !=====================================================================
   ! Write of the results (all data are in X-pencil)
   !=====================================================================
   subroutine write_data()

      use decomp_2d_io

      implicit none

      character(len=*), parameter :: io_name = "grad-io"
#ifndef ADIOS2
      logical :: dir_exists
#endif

      !$acc update self (phi1)
      !$acc update self (dphiX)
      !$acc update self (dphiY)
      !$acc update self (dphiZ)
#ifndef ADIOS2
      if (nrank == 0) then
         inquire (file="out", exist=dir_exists)
         if (.not. dir_exists) then
            call execute_command_line("mkdir out 2> /dev/null")
         end if
      end if
#endif

      call decomp_2d_io_init()
      call decomp_2d_init_io(io_name)

      call decomp_2d_register_variable(io_name, "phi1.dat", 1, 0, 0, data_type, mytype)
      call decomp_2d_register_variable(io_name, "dphiX.dat", 1, 0, 0, data_type, mytype)
      call decomp_2d_register_variable(io_name, "dphiY.dat", 1, 0, 0, data_type, mytype)
      call decomp_2d_register_variable(io_name, "dphiZ.dat", 1, 0, 0, data_type, mytype)

      ! Standard I/O pattern - file per field
#ifdef ADIOS2
      call decomp_2d_open_io(io_name, "out", decomp_2d_write_mode)
      call decomp_2d_start_io(io_name, "out")
#endif
      call decomp_2d_write_one(1, phi1, 'out', 'phi1.dat', 0, io_name)
      call decomp_2d_write_one(1, dphiX, 'out', 'dphiX.dat', 0, io_name)
      call decomp_2d_write_one(1, dphiY, 'out', 'dphiY.dat', 0, io_name)
      call decomp_2d_write_one(1, dphiZ, 'out', 'dphiZ.dat', 0, io_name)
#ifdef ADIOS2
      call decomp_2d_end_io(io_name, "out")
      call decomp_2d_close_io(io_name, "out")
#else
      call write_xdmf()
#endif

   end subroutine write_data
   !=====================================================================
   ! Write of the xdmf file to visualise in paraview
   !=====================================================================
   subroutine write_xdmf()
      ! This subroutine is based on the xdmf writers in Xcompact3d.
      ! Copyright (c) 2012-2022, Xcompact3d
      ! SPDX-License-Identifier: BSD 3-Clause

      integer :: ioxdmf

      character(len=:), allocatable :: fmt

      integer :: precision
      integer, parameter :: output2D = 0 ! Which plane to write in 2D (0 for 3D)

      integer :: varctr
      character(len=16) :: filename
      character(len=5) :: varname
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
         if (mytype == kind(0.0d0)) then
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
         do varctr = 1, 4
            select case (varctr)
            case (1)
               write (varname, "(A)") "phi1"
               write (filename, '(A)') "./out/phi1.dat"
            case (2)
               write (varname, "(A)") "dphiX"
               write (filename, '(A)') "./out/dphiX.dat"
            case (3)
               write (varname, "(A)") "dphiY"
               write (filename, '(A)') "./out/dphiY.dat"
            case (4)
               write (varname, "(A)") "dphiZ"
               write (filename, '(A)') "./out/dphiZ.dat"
            end select
            write (ioxdmf, *) '       <Attribute Name="'//trim(varname)//'" Center="Node">'
#ifndef ADIOS2
            write (ioxdmf, *) '          <DataItem Format="Binary"'
#else
            write (ioxdmf, *) '          <DataItem Format="HDF"'
#endif

#ifdef DOUBLE_PREC
            print *, "Double precision build"
#ifdef SAVE_SINGLE
            if (output2D == 0) then
               precision = 4
            else
               precision = 8
            end if
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

   end subroutine write_xdmf
end program grad3d
