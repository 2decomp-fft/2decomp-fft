#!/usr/bin/env python3

# Name of the output file
outfile = "read_var.f90"

# Number of variables
nformat = 6
# Suffix for the subroutine name
ext = ["freal", "fcplx", "dreal", "dcplx", "ints", "logs"]


#
# Open the file
#
f = open(outfile, "w")

#
# Write the interface
#
f.write("\n")
f.write("   interface decomp_2d_read_var\n")
for i in range(nformat):
    f.write("      module procedure read_var_"+ext[i]+"\n")

f.write("   end interface decomp_2d_read_var\n")
f.write("\n")

for i in range(nformat):
    #
    # Header
    #
    f.write("   subroutine read_var_"+ext[i]+"(io, ipencil, var, &\n")
    f.write("                             opt_reduce_prec, &\n")
    f.write("                             opt_decomp, &\n")
    f.write("                             opt_nb_req)\n")
    f.write("\n")
    f.write("      implicit none\n")
    f.write("\n")
    #
    # Arguments
    #
    f.write("      ! Arguments\n")
    f.write("      type(d2d_io_mpi), intent(INOUT) :: io\n")
    f.write("      integer, intent(IN) :: ipencil\n")
    if (i==0):
        f.write("      real(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==2):
        f.write("      real(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==1):
        f.write("      complex(real32), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==3):
        f.write("      complex(real64), contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==4):
        f.write("      integer, contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    elif (i==5):
        f.write("      logical, contiguous, dimension(:, :, :), intent(OUT) :: var\n")
    #
    f.write("      logical, intent(in), optional :: opt_reduce_prec\n")
    f.write("      TYPE(DECOMP_INFO), target, intent(IN), optional :: opt_decomp\n")
    f.write("      integer, intent(inout), optional :: opt_nb_req\n")
    f.write("\n")
    #
    # Local variables
    #
    f.write("      ! Local variable(s)\n")
    if (i==2 or i==3):
        f.write("      logical :: reduce\n")
    #
    f.write("      TYPE(DECOMP_INFO), pointer :: decomp\n")
    if (i==2):
        f.write("      real(real32), allocatable, dimension(:, :, :) :: tmp\n")
    elif (i==3):
        f.write("      complex(real32), allocatable, dimension(:, :, :) :: tmp\n")
    f.write("\n")
    #
    # Start profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_start(\"io_read_var\")\n")
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
    if (i==0 or i==1 or i==4 or i==5):
        f.write("      call read_var(io, ipencil, decomp, &\n")
        f.write("                    opt_nb_req=opt_nb_req, &\n")
        if (i==0):
            f.write("                    freal=var)\n")
        elif (i==1):
            f.write("                    fcplx=var)\n")
        elif (i==4):
            f.write("                    ints=var)\n")
        elif (i==5):
            f.write("                    logs=var)\n")
    #
    if (i==2 or i==3):
        f.write("      if (reduce) then\n")
        f.write("         if (ipencil == 1) then\n")
        f.write("            call alloc_x(tmp, decomp)\n")
        f.write("         else if (ipencil == 2) then\n")
        f.write("            call alloc_y(tmp, decomp)\n")
        f.write("         else\n")
        f.write("            call alloc_z(tmp, decomp)\n")
        f.write("         end if\n")
        f.write("         call read_var(io, ipencil, decomp, &\n")
        if (i==2):
            f.write("                      freal=tmp)\n")
            f.write("         var = real(tmp, kind=real32)\n")
        elif (i==3):
            f.write("                      fcplx=tmp)\n")
            f.write("         var = cmplx(tmp, kind=real64)\n")
        f.write("         deallocate(tmp)\n")
        f.write("      else\n")
        f.write("         call read_var(io, ipencil, decomp, &\n")
        f.write("                       opt_nb_req=opt_nb_req, &\n")
        if (i==2):
            f.write("                           dreal=var)\n")
        elif (i==3):
            f.write("                           dcplx=var)\n")
        f.write("      end if\n")
    #
    # Cleaning + unused variables
    #
    f.write("\n")
    f.write("      nullify (decomp)\n")
    f.write("\n")
    if (i==0 or i==1 or i==4 or i==5):
        f.write("      associate (p => opt_reduce_prec); end associate\n")
        f.write("\n")
    #
    # End profiling
    #
    f.write("      if (decomp_profiler_io) call decomp_profiler_end(\"io_read_var\")\n")
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
