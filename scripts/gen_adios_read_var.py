#!/usr/bin/env python3

# Name of the output file
outfile = "adios_read_var.f90"

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
f.write("   interface decomp_2d_adios_read_var\n")
for i in range(nformat):
    f.write("      module procedure read_var_"+ext[i]+"\n")

f.write("   end interface decomp_2d_adios_read_var\n")
f.write("\n")

for i in range(nformat):
    #
    # Header
    #
    f.write("   subroutine read_var_"+ext[i]+"(var, varname, &\n")
    f.write("                             opt_reduce_prec, &\n")
    f.write("                             opt_ipencil, &\n")
    f.write("                             opt_decomp, &\n")
    f.write("                             opt_family, &\n")
    f.write("                             opt_io)\n")
    f.write("\n")
    f.write("      implicit none\n")
    f.write("\n")
    #
    # Arguments
    #
    f.write("      ! Arguments\n")
    if (i==0):
        f.write("      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==2):
        f.write("      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==1):
        f.write("      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==3):
        f.write("      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    #
    f.write("      character(len=*), intent(in) :: varname\n")
    f.write("      logical, intent(in), optional :: opt_reduce_prec\n")
    f.write("      integer, intent(in), optional :: opt_ipencil\n")
    f.write("      type(decomp_info), target, intent(in), optional :: opt_decomp\n")
    f.write("      type(d2d_io_family), target, intent(in), optional :: opt_family\n")
    f.write("      type(d2d_io_adios), intent(inout), optional :: opt_io")
    f.write("\n")
    #
    # Local variables
    #
    if (i==2 or i==3):
        f.write("      ! Local variable(s)\n")
        f.write("      logical :: reduce\n")
        f.write("      type(decomp_info), pointer :: decomp\n")
        if (i==2):
            f.write("      real(real32), dimension(:, :, :), allocatable :: tmp\n")
        elif (i==3):
            f.write("      complex(real32), dimension(:, :, :), allocatable :: tmp\n")
    f.write("\n")
    #
    # Start profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_start(\"adios_read_var\")\n")
    f.write("\n")
    #
    # Deal with optional arguments
    #
    if (i==2 or i==3):
        f.write("      ! One can write to single precision using opt_reduce_prec\n")
        f.write("      if (present(opt_reduce_prec)) then\n")
        f.write("         reduce = opt_reduce_prec\n")
        f.write("      else\n")
        f.write("         reduce = default_opt_reduce_prec\n")
        f.write("      end if\n")
        f.write("\n")
        f.write("      if (reduce) then\n")
        f.write("         if (present(opt_decomp)) then\n")
        f.write("            decomp => opt_decomp\n")
        f.write("         else\n")
        f.write("            decomp => decomp_main\n")
        f.write("         end if\n")
        f.write("         ! Safety check\n")
        f.write("         if (.not.present(opt_ipencil)) then\n")
        f.write("            call decomp_2d_abort(__FILE__, __LINE__, 1, \"Invalid arguments\")\n")
        f.write("         else\n")
        f.write("            if (opt_ipencil==1) then\n")
        f.write("               call alloc_x(tmp, decomp)\n")
        f.write("            elseif (opt_ipencil==2) then\n")
        f.write("               call alloc_y(tmp, decomp)\n")
        f.write("            elseif (opt_ipencil==3) then\n")
        f.write("               call alloc_z(tmp, decomp)\n")
        f.write("            else\n")
        f.write("               call decomp_2d_abort(__FILE__, __LINE__, opt_ipencil, \"Invalid arguments\")\n")
        f.write("            end if\n")
        f.write("            nullify(decomp)\n")
        f.write("         end if\n")
        f.write("      end if\n")
    #
    # Call the lower level IO subroutine
    #
    if (i==0 or i==1):
        f.write("      call adios_read(varname, &\n")
        f.write("                      opt_family=opt_family, &\n")
        f.write("                      opt_io=opt_io, &\n")
        if (i==0):
            f.write("                      freal=var)\n")
        elif (i==1):
            f.write("                      fcplx=var)\n")
    #
    if (i==2 or i==3):
        f.write("      if (reduce) then\n")
        f.write("         call adios_read(varname, &\n")
        f.write("                         opt_family=opt_family, &\n")
        f.write("                         opt_io=opt_io, &\n")
        if (i==2):
            f.write("                         freal=tmp)\n")
        elif (i==3):
            f.write("                         fcplx=tmp)\n")
        f.write("         var = tmp\n")
        f.write("         deallocate(tmp)\n")
        f.write("      else\n")
        f.write("         call adios_read(varname, &\n")
        f.write("                         opt_family=opt_family, &\n")
        f.write("                         opt_io=opt_io, &\n")
        if (i==2):
            f.write("                         dreal=var)\n")
        elif (i==3):
            f.write("                         dcplx=var)\n")
        f.write("      end if\n")
    #
    # Cleaning + unused variables
    #
    f.write("\n")
    if (i==0 or i==1):
        f.write("      associate (p => opt_reduce_prec); end associate\n")
        f.write("      associate (p => opt_ipencil); end associate\n")
        f.write("      associate (p => opt_decomp); end associate\n")
        f.write("\n")
    #
    # End profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_end(\"adios_read_var\")\n")
    f.write("\n")
    #
    # Footer
    #
    f.write("   end subroutine read_var_"+ext[i]+"\n")
    f.write("   !\n")

#
# Closing the file
#
f.close()
