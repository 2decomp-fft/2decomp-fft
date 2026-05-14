# Compilers Flags for AOCC
set(D2D_FFLAGS "-cpp -g")
#set(D2D_FFLAGS_RELEASE "-O3")
set(D2D_FFLAGS_RELEASE "-O2 -fopenmp -cpp -I/shared/apps/ubuntu/opt/rocm-7.2.0/include/hipfort/amdgcn --offload-arch=gfx942")
set(D2D_FFLAGS_DEBUG   "-g -O0 -DDEBUG")
set(D2D_FFLAGS_DEV     "${D2D_FFLAGS_DEBUG}")
