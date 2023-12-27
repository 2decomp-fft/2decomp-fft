#!/usr/bin/env python3

# Name of the output file
outfile = "write_one.f90"

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
f.write("   interface decomp_2d_write_one\n")
for i in range(nformat):
    f.write("      module procedure write_one_"+ext[i]+"\n")

f.write("   end interface decomp_2d_write_one\n")
f.write("\n")

for i in range(nformat):
    #
    # Header
    #
    f.write("   subroutine write_one_"+ext[i]+"(ipencil, var, varname, opt_mode, opt_family, &\n")
    f.write("                              opt_io, opt_dirname, opt_reduce_prec, opt_decomp)\n")
    f.write("\n")
    f.write("      implicit none\n")
    f.write("\n")
    #
    # Arguments
    #
    f.write("      ! Arguments\n")
    f.write("      integer, intent(IN) :: ipencil\n")
    if (i==0):
        f.write("      real(real32), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    elif (i==2):
        f.write("      real(real64), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    elif (i==1):
        f.write("      complex(real32), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    elif (i==3):
        f.write("      complex(real64), contiguous, dimension(:, :, :), intent(IN) :: var\n")
    #
    f.write("      character(len=*), intent(in) :: varname\n")
    f.write("      integer, intent(in), optional :: opt_mode\n")
    f.write("      type(d2d_io_family), target, intent(in), optional :: opt_family\n")
    f.write("      type(d2d_io), intent(inout), optional :: opt_io\n")
    f.write("      character(len=*), intent(in), optional :: opt_dirname\n")
    f.write("      logical, intent(in), optional :: opt_reduce_prec\n")
    f.write("      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp\n")
    f.write("\n")
    #
    # Local variables
    #
    f.write("      ! Local variable(s)\n")
    if (i==2 or i==3):
        f.write("      logical :: reduce\n")
    #
    f.write("      TYPE(DECOMP_INFO), pointer :: decomp\n")
    f.write("\n")
    #
    # Start profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_start(\"io_write_one\")\n")
    f.write("\n")
    #
    # Deal with optional arguments
    #
    f.write("      ! Use the provided decomp_info or the default one\n")
    f.write("      if (present(opt_decomp)) then\n")
    f.write("         decomp => opt_decomp\n")
    f.write("      else\n")
    f.write("         decomp => decomp_main\n")
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
        f.write("      call write_one(ipencil, varname, decomp, &\n")
        f.write("                     opt_mode=opt_mode, &\n")
        f.write("                     opt_family=opt_family, &\n")
        f.write("                     opt_io=opt_io, &\n")
        f.write("                     opt_dirname=opt_dirname, &\n")
        if (i==0):
            f.write("                     freal=var)\n")
        elif (i==1):
            f.write("                     fcplx=var)\n")
    #
    if (i==2 or i==3):
        f.write("      if (reduce) then\n")
        f.write("         call write_one(ipencil, varname, decomp, &\n")
        f.write("                        opt_mode=decomp_2d_io_sync, &\n")
        f.write("                        opt_family=opt_family, &\n")
        f.write("                        opt_io=opt_io, &\n")
        f.write("                        opt_dirname=opt_dirname, &\n")
        if (i==2):
            f.write("                        freal=real(var, kind=real32)) ! Warning, implicit memory allocation\n")
        elif (i==3):
            f.write("                        fcplx=cmplx(var, kind=real32)) ! Warning, implicit memory allocation\n")
        f.write("      else\n")
        f.write("         call write_one(ipencil, varname, decomp, &\n")
        f.write("                        opt_mode=opt_mode, &\n")
        f.write("                        opt_family=opt_family, &\n")
        f.write("                        opt_io=opt_io, &\n")
        f.write("                        opt_dirname=opt_dirname, &\n")
        if (i==2):
            f.write("                        dreal=var)\n")
        elif (i==3):
            f.write("                        dcplx=var)\n")
        f.write("      end if\n")
    #
    # Cleaning + unused variables
    #
    f.write("\n")
    f.write("      nullify (decomp)\n")
    f.write("\n")
    if (i==0 or i==1):
        f.write("      associate (p => opt_reduce_prec)\n")
        f.write("      end associate\n")
        f.write("\n")
    #
    # End profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_end(\"io_write_one\")\n")
    f.write("\n")
    #
    # Footer
    #
    f.write("   end subroutine write_one_"+ext[i]+"\n")
    f.write("   !\n")

#
# Closing the file
#
f.close()
