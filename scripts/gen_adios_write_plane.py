#!/usr/bin/env python3

# Name of the output file
outfile = "adios_write_plane.f90"

# Number of variables
nformat = 4
# Suffix for the subroutine name
ext = ["freal", "fcplx", "dreal", "dcplx"]


#
# Open the file
#
f = open(outfile, "w")

#
# Write the interface
#
f.write("\n")
f.write("   interface decomp_2d_adios_write_plane\n")
for i in range(nformat):
    f.write("      module procedure write_plane_"+ext[i]+"\n")

f.write("   end interface decomp_2d_adios_write_plane\n")
f.write("\n")

for i in range(nformat):
    #
    # Header
    #
    f.write("   subroutine write_plane_"+ext[i]+"(io, var, varname, &\n")
    f.write("                                opt_mode, &\n")
    f.write("                                opt_family, &\n")
    f.write("                                opt_nplanes, &\n")
    f.write("                                opt_iplane, &\n")
    f.write("                                opt_reduce_prec, &\n")
    f.write("                                opt_ipencil)\n")
    f.write("\n")
    f.write("      implicit none\n")
    f.write("\n")
    #
    # Arguments
    #
    f.write("      ! Arguments\n")
    f.write("      type(d2d_io_adios), intent(inout) :: io\n")
    if (i==0):
        f.write("      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    elif (i==2):
        f.write("      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    elif (i==1):
        f.write("      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    elif (i==3):
        f.write("      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    f.write("      character(len=*), intent(in) :: varname\n")
    f.write("      integer, intent(in), optional :: opt_mode\n")
    f.write("      type(d2d_io_family), intent(inout), optional :: opt_family\n")
    f.write("      integer, intent(in), optional :: opt_nplanes\n")
    f.write("      integer, intent(in), optional :: opt_iplane\n")
    f.write("      logical, intent(in), optional :: opt_reduce_prec\n")
    f.write("      integer, intent(in), optional :: opt_ipencil\n")
    f.write("\n")
    #
    # Local variables
    #
    f.write("      ! Local variables\n")
    if (i==2 or i==3):
        f.write("      logical :: reduce\n")
    if (i==0):
        f.write("      real(real32), allocatable, dimension(:, :, :) :: var2d\n")
    elif (i==1):
        f.write("      complex(real32), allocatable, dimension(:, :, :) :: var2d\n")
    elif (i==2):
        f.write("      real(real64), allocatable, dimension(:, :, :) :: var2d\n")
        f.write("      real(real32), allocatable, dimension(:, :, :) :: var2dbis\n")
    elif (i==3):
        f.write("      complex(real64), allocatable, dimension(:, :, :) :: var2d\n")
        f.write("      complex(real32), allocatable, dimension(:, :, :) :: var2dbis\n")
    f.write("      integer :: iplane, ipencil\n")
    f.write("\n")
    #
    # Start profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_start(\"adios_write_plane\")\n")
    f.write("\n")
    #
    # Deal with optional arguments
    #
    f.write("      if (present(opt_iplane)) then\n")
    f.write("         iplane = opt_iplane\n")
    f.write("      else\n")
    f.write("         iplane = 1\n")
    f.write("      end if\n")
    f.write("      ipencil = 0\n")
    f.write("      if (present(opt_ipencil)) ipencil = opt_ipencil\n")
    f.write("\n")
    f.write("      ! Safety check\n")
    f.write("      if (.not. (present(opt_nplanes) .or. present(opt_ipencil))) then\n")
    f.write("         call decomp_2d_abort(__FILE__, __LINE__, -1, \"Invalid arguments\")\n")
    f.write("      end if\n")
    f.write("\n")
    if (i==2 or i==3):
        f.write("      ! One can write to single precision using opt_reduce_prec\n")
        f.write("      if (present(opt_reduce_prec)) then\n")
        f.write("         reduce = opt_reduce_prec\n")
        f.write("      else\n")
        f.write("         reduce = default_opt_reduce_prec\n")
        f.write("      end if\n")
        f.write("\n")
    #
    # Call the lower level IO subroutine
    #
    if (i==0 or i==1):
        f.write("      if (present(opt_nplanes)) then\n")
        f.write("         call adios_write(io, varname, &\n")
        f.write("                          opt_mode=opt_mode, &\n")
        f.write("                          opt_family=opt_family, &\n")
        if (i==0):
            f.write("                          freal=var)\n")
        elif (i==1):
            f.write("                          fcplx=var)\n")
        f.write("      else\n")
        f.write("         if (opt_ipencil == 1) then\n")
        f.write("            allocate (var2d(1, size(var, 2), size(var, 3)))\n")
        f.write("            if (iplane < 1) then\n")
        f.write("               var2d(1, :, :) = sum(var(:, :, :), dim=ipencil)\n")
        f.write("            else\n")
        f.write("               var2d(1, :, :) = var(iplane, :, :)\n")
        f.write("            end if\n")
        f.write("         else if (opt_ipencil == 2) then\n")
        f.write("            allocate (var2d(size(var, 1), 1, size(var, 3)))\n")
        f.write("            if (iplane < 1) then\n")
        f.write("               var2d(:, 1, :) = sum(var(:, :, :), dim=ipencil)\n")
        f.write("            else\n")
        f.write("               var2d(:, 1, :) = var(:, iplane, :)\n")
        f.write("            end if\n")
        f.write("         else if (opt_ipencil == 3) then\n")
        f.write("            allocate (var2d(size(var, 1), size(var, 2), 1))\n")
        f.write("            if (iplane < 1) then\n")
        f.write("               var2d(:, :, 1) = sum(var(:, :, :), dim=ipencil)\n")
        f.write("            else\n")
        f.write("               var2d(:, :, 1) = var(:, :, iplane)\n")
        f.write("            end if\n")
        f.write("         end if\n")
        f.write("         call adios_write(io, varname, &\n")
        f.write("                          opt_mode=decomp_2d_io_sync, &\n")
        f.write("                          opt_family=opt_family, &\n")
        if (i==0):
            f.write("                          freal=var2d)\n")
        elif (i==1):
            f.write("                          fcplx=var2d)\n")
        f.write("         deallocate (var2d)\n")
        f.write("      end if\n")
    #
    if (i==2 or i==3):
        f.write("      if (present(opt_nplanes)) then\n")
        f.write("         if (reduce) then\n")
        f.write("            allocate (var2dbis(size(var, 1), &\n")
        f.write("                               size(var, 2), &\n")
        f.write("                               size(var, 3)))\n")
        if (i==2):
            f.write("            var2dbis = real(var, kind=real32)\n")
        if (i==3):
            f.write("            var2dbis = cmplx(var, kind=real32)\n")
        f.write("            call adios_write(io, varname, &\n")
        f.write("                             opt_mode=decomp_2d_io_sync, &\n")
        f.write("                             opt_family=opt_family, &\n")
        if (i==2):
            f.write("                             freal=var2dbis)\n")
        elif (i==3):
            f.write("                             fcplx=var2dbis)\n")
        f.write("         else\n")
        f.write("            call adios_write(io, varname, &\n")
        f.write("                             opt_mode=opt_mode, &\n")
        f.write("                             opt_family=opt_family, &\n")
        if (i==2):
            f.write("                             dreal=var)\n")
        elif (i==3):
            f.write("                             dcplx=var)\n")
        f.write("         end if\n")
        f.write("      else\n")
        f.write("         if (reduce .and. opt_ipencil == 1) then\n")
        f.write("            allocate (var2dbis(1, size(var, 2), size(var, 3)))\n")
        f.write("            if (iplane < 1) then\n")
        if (i==2):
            f.write("               var2dbis(1, :, :) = real(sum(var(:, :, :), dim=ipencil), kind=real32)\n")
        else:
            f.write("               var2dbis(1, :, :) = cmplx(sum(var(:, :, :), dim=ipencil), kind=real32)\n")
        f.write("            else\n")
        if (i==2):
            f.write("               var2dbis(1, :, :) = real(var(iplane, :, :), kind=real32)\n")
        else:
            f.write("               var2dbis(1, :, :) = cmplx(var(iplane, :, :), kind=real32)\n")
        f.write("            end if\n")
        f.write("         else if (reduce .and. opt_ipencil == 2) then\n")
        f.write("            allocate (var2dbis(size(var, 1), 1, size(var, 3)))\n")
        f.write("            if (iplane < 1) then\n")
        if (i==2):
            f.write("               var2dbis(:, 1, :) = real(sum(var(:, :, :), dim=ipencil), kind=real32)\n")
        else:
            f.write("               var2dbis(:, 1, :) = cmplx(sum(var(:, :, :), dim=ipencil), kind=real32)\n")
        f.write("            else\n")
        if (i==2):
            f.write("               var2dbis(:, 1, :) = real(var(:, iplane, :), kind=real32)\n")
        else:
            f.write("               var2dbis(:, 1, :) = cmplx(var(:, iplane, :), kind=real32)\n")
        f.write("            end if\n")
        f.write("         else if (reduce .and. opt_ipencil == 3) then\n")
        f.write("            allocate (var2dbis(size(var, 1), size(var, 2), 1))\n")
        f.write("            if (iplane < 1) then\n")
        if (i==2):
            f.write("               var2dbis(:, :, 1) = real(sum(var(:, :, :), dim=ipencil), kind=real32)\n")
        else:
            f.write("               var2dbis(:, :, 1) = cmplx(sum(var(:, :, :), dim=ipencil), kind=real32)\n")
        f.write("            else\n")
        if (i==2):
            f.write("               var2dbis(:, :, 1) = real(var(:, :, iplane), kind=real32)\n")
        else:
            f.write("               var2dbis(:, :, 1) = cmplx(var(:, :, iplane), kind=real32)\n")
        f.write("            end if\n")
        f.write("         else if (opt_ipencil == 1) then\n")
        f.write("            allocate (var2d(1, size(var, 2), size(var, 3)))\n")
        f.write("            if (iplane < 1) then\n")
        f.write("               var2d(1, :, :) = sum(var(:, :, :), dim=ipencil)\n")
        f.write("            else\n")
        f.write("               var2d(1, :, :) = var(iplane, :, :)\n")
        f.write("            end if\n")
        f.write("         else if (opt_ipencil == 2) then\n")
        f.write("            allocate (var2d(size(var, 1), 1, size(var, 3)))\n")
        f.write("            if (iplane < 1) then\n")
        f.write("               var2d(:, 1, :) = sum(var(:, :, :), dim=ipencil)\n")
        f.write("            else\n")
        f.write("               var2d(:, 1, :) = var(:, iplane, :)\n")
        f.write("            end if\n")
        f.write("         else if (opt_ipencil == 3) then\n")
        f.write("            allocate (var2d(size(var, 1), size(var, 2), 1))\n")
        f.write("            if (iplane < 1) then\n")
        f.write("               var2d(:, :, 1) = sum(var(:, :, :), dim=ipencil)\n")
        f.write("            else\n")
        f.write("               var2d(:, :, 1) = var(:, :, iplane)\n")
        f.write("            end if\n")
        f.write("         end if\n")
        f.write("         if (reduce) then\n")
        f.write("            call adios_write(io, varname, &\n")
        f.write("                             opt_mode=decomp_2d_io_sync, &\n")
        f.write("                             opt_family=opt_family, &\n")
        if (i==2):
            f.write("                             freal=var2dbis)\n")
        elif (i==3):
            f.write("                             fcplx=var2dbis)\n")
        f.write("         else\n")
        f.write("            call adios_write(io, varname, &\n")
        f.write("                             opt_mode=decomp_2d_io_sync, &\n")
        f.write("                             opt_family=opt_family, &\n")
        if (i==2):
            f.write("                             dreal=var2d)\n")
        elif (i==3):
            f.write("                             dcplx=var2d)\n")
        f.write("         end if\n")
        f.write("      end if\n")
        f.write("      if (allocated(var2d)) deallocate (var2d)\n")
        f.write("      if (allocated(var2dbis)) deallocate (var2dbis)\n")
    #
    # Cleaning + unused variables
    #
    f.write("\n")
    if (i==0 or i==1):
        f.write("      associate (p => opt_reduce_prec); end associate\n")
        f.write("\n")
    #
    # End profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_end(\"adios_write_plane\")\n")
    f.write("\n")
    #
    # Footer
    #
    f.write("   end subroutine write_plane_"+ext[i]+"\n")
    f.write("   !\n")

#
# Closing the file
#
f.close()
