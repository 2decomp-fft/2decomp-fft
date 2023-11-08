---
title: 'The 2DECOMP&FFT library: an update with new CPU/GPU capabilities'
tags:
  - Fortran
  - FFT
  - Finite Difference Discretisation
  - Pencil Decomposition
authors:
  - name: Stefano Rolfo
    orcid: 0000-0001-6325-7629
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Cédric Flageul
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Paul Bartholomew
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 3
  - name: Filippo Spiga
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 4
  - name: Sylvain Laizet
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 5
affiliations:
 - name: STFC Daresbury Laboratory, Scientific Computing Department, UKRI, UK 
   index: 1
 - name: PPRIME institute, Curiosity Group, Université de Poitiers, CNRS, ISAE-ENSMA, Poitiers, France
   index: 2
 - name: The University of Edinburgh, EPCC, Edinburgh, UK
   index: 3
 - name: NVIDIA Corporation, Cambridge, UK
   index: 4
 - name: Department of Aeronautics, Imperial College London, London, UK
   index: 5
date: 31 July 2023
bibliography: paper.bib

---

# Summary

The 2DECOMP&FFT library is a software framework written in modern Fortran to build large-scale parallel applications. 
It is designed for applications using three-dimensional structured meshes with a particular focus on 
spatially implicit numerical algorithms. However, the library can be easily used also with other discretisation schemes 
based on a structured layout and where pencil decomposition can apply. 
It is based on a general-purpose 2D pencil decomposition for data distribution and data Input Output (I/O). 
A 1D slab decomposition is also available as a special case of the 2D pencil decomposition.
The library includes a highly scalable and efficient interface to perform three-dimensional 
Fast Fourier Transforms (FFTs). 
The library has been designed to be user-friendly, with a clean application programming interface 
hiding most communication details from application developers, and portable with support for modern CPUs 
and NVIDIA GPUs (support for AMD and Intel GPUs to follow).


# Statement of need

The 2DECOMP&FFT library [@li_2010_2decomp] was originally designed for CPU hardware 
and is now used by many research groups worldwide. 
The library is based on a 2D-pencil decomposition for data distribution on distributed memory systems 
and is used as the core of many CFD solvers such as Xcompact3d [@BARTHOLOMEW2020100550] or CaNS [@COSTA20181853], 
with excellent strong scaling performance up to hundreds of thousands of CPU cores. 
2DECOMP&FFT mainly relies on MPI, and it offers a user-friendly interface that hides the complexity 
of the communication. 
The version 2.0.0 of the library offers also a 1D slab decomposition, which is implemented as a special case of the 
2D decomposition. Two alternatives are possible:

- Initial slabs orientation in the `XY` plane; 
- Initial slabs orientation in the `XZ` plane.

In many configurations the slabs decomposition gives some gain in performances with respect to the 
2D-pencil decomposition. This is a consequence of having data already in memory when transposing between 
the two dimensions of the slab. Therefore, it is possible to perform a simple memory copy between 
input and output arrays instead of the full MPI communication.

The library also offers also a very efficient and flexible 3D tool to perform 
Fast Fourier Transform (FFT) for distributed memory systems. However, 2DECOMP&FFT is mainly designed to perform 
data management, communication and the actual computation of the 1D FFT is delegated to 3rd-parly libraries. 
The supported FFT backends are: FFTW [@FFTW05], the Intel Math Kernel Library (MKL) and the CUDA FFT (cuFFT) 
which is used for FFT on NVIDIA GPUs. A Generic FFT backend, based on the 
Glassman's general N Fast Fourier Transform [@FERGUSON1982401], 
is also available to make the library more portable.   

While 2DECOMP&FFT library has been designed with high order compact schemes in mind, it is possible 
that some derivatives can be evaluated using explicit formulation using local stencils. For this reason a 
halo support API is also provided to support explicit message passing between neighbouring pencils. 

Finally, the library provides an infrastructure to perform data I/O using MPI I/O or ADIOS2 [@godoy2020adios]. 
The API provide several features such as: write of a single or multiple 3D arrays into a file, write of 2D slices 
from a 2D array, write of the results using a lower precision (i.e. double precision into single) 
and with lower resolution. 

The first version of the library was released in 2010 as a tar.gz package, with a Makefile approach, 
and could only make use of CPUs. It has not been modified since its release. 
The new version of the library can now leverage NVIDIA GPUs, modern CPUs and various compilers 
(GNU, Intel, NVHPC, CRAY). 
It has CMAKE capabilities as well as a proper continuous integration framework with automated tests. 
The new library was designed to be more appealing to the scientific community,
and can now be easily implemented as an independent library for use by other software.

# GPU porting

An initial port of 2DECOMP&FFT to GPU has been performed within the solver AFiD-GPU [@ZHU2018199],
which was mainly based on CUDA-Fortran for some kernels and CUDA-aware-MPI for communications.
A second library, named cuDECOMP, which is directly inspired by 2DECOMP&FFT, 
takes full advantages of CUDA and uses NVIDIA most recent libraries for communications 
such as NVIDIA Collective Communication Library (NCCL), is presented in [@Romero_2022_cuDecomp].
Indeed, cuDECOMP only targets NVIDIA GPUs.
The updated 2DECOMP&FFT mainly uses a mix of CUDA-Fortran and openACC for the GPU porting 
together with CUDA-aware-MPI and NCCL for the communications. In addition to previous work,
the FFT module is ported to GPUs using cuFFT. 
The next step is also to implement OpenMP for GPU porting to support both AMD and Intel GPU hardware.

# How to use 2DECOMP&FFT
The 2D Pencil Decomposition API is defined with three Fortran module which should be used by applications as:
```
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
```
where ``use decomp_2d_constants`` defines all the parameters, ``use decomp_2d_mpi`` introduces all the MPI 
related interfaces and ``use decomp_2d`` contains the main decomposition and transposition APIs.
The library is initialized using:
```
  call decomp_2d_init(nx, ny, nz, p_row, p_col)
```
where ``nx``, ``ny`` and ``nz`` are the spatial dimensions of the problem, to be distributed over
a 2D processor grid $p_{row} \times p_{col}$.
Note that none of the dimensions need to be divisible by ``p_row`` or ``p_col``.
In case of ``p_row=p_col=0`` an automatic decomposition is selected among all possible combination available.
A key element of this library is a set of communication routines that perform the data transpositions.
As mentioned, one needs to perform 4 global transpositions to go through all 3 pencil orientations
(i.e. one has to go from x-pencils to y-pencils to z-pencils to y-pencils to x-pencils).
Correspondingly, the library provides 4 communication subroutines:
```
  call transpose_x_to_y(var_in,var_out)
  call transpose_y_to_z(var_in,var_out)
  call transpose_z_to_y(var_in,var_out)
  call transpose_y_to_x(var_in,var_out)
```
The input array ``var_in`` and output array ``var_out`` are defined by the code using the library
and contain distributed data for the correct pencil orientations.

Note that the library is written using Fortran's generic interface so different data types are supported
without user input. That means in and out above can be either real or complex arrays,
the latter being useful for applications involving 3D Fast Fourier Transforms.
Finally, before exit, applications should clean up the memory by:
```
  call decomp_2d_finalize
```
Detailed information about the decomposition API are available [here](https://2decomp-fft.github.io/pages/api_domain.html) 

# Acknowledgements

The first version of the library was initially designed thanks to several projects funded under the 
HECToR Distributed Computational Science and Engineering (CSE) Service operated by NAG Ltd. 
The new library has been designed thanks to the support of EPSRC via the CCP Turbulence (EP/T026170/1) and work funded under
the embedded CSE programme of the ARCHER2 UK National Supercomputing Service (http://www.archer2.ac.uk) (ARCHER2 eCSE03-2).

# References

