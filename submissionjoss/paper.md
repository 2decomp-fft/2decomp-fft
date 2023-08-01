---
title: 'The new version of 2DECOMP&FFT'
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
  - name: C\'edric Flageul
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Paul Bartholomew
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
  - name: Filippo Spiga
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 4
  - name: Sylvain Laizet
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 5
affiliations:
 - name: STFC Daresbury Laboratory, Scientific Computing Department, UKRI, UK 
   index: 1
 - name: Universit\'e de Poitiers, CNRS,ISAE-ENSMA, CNRS, PPRIME, Poitiers, France 
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

The 2DECOMP&FFT library is a software framework written in modern Fortran to build large-scale parallel applications. It is designed for applications using three-dimensional structured mesh and spatially implicit numerical algorithms. It is based on a general-purpose 2D pencil decomposition for data distribution. It includes a highly scalable and efficient interface to perform three-dimensional Fast Fourier Transforms (FFTs). It relies on MPI but provides a user-friendly programming interface that hides communication details from application developers. This library has been designed to be user friendly, with a clean application programming interface hiding most communication details from applications, and portable (modern CPUs and NVIDIA GPUs; support for AMD and Intel GPUs to follow).


# Statement of need

The 2DECOMP&FFT library [@li_2010_2decomp] is a software framework written in Fortran targeting large-scale parallel applications on structured meshes. It was originally designed for CPU hardware and is now used by many research groups worldwide. The library is based on a 2D-pencil decomposition for data distribution on distributed memory systems and is used as the core of many CFD solvers such as Xcompact3d [@BARTHOLOMEW2020100550], with excellent strong scaling performance up to hundreds of thousands of CPU cores. 
2DECOMP&FFT mainly relies on MPI, and it offers a user-friendly interface that hides the complexity of the communication. The library also offers also a very efficient and flexible 3D tool to perform Fast Fourier Transform (FFT) using popular FFT libraries like FFTW. 

The first version of the library was released in 2010 as a tar.gz package, with a Makefile approach, and could only make use of CPUs. It has not been modified since its release. The new version of the library can now leverage NVIDIA GPUs, modern CPUs and various compilers. It has CMAKE capabilities as well as a proper continuous integration framework with automated tests. The new library was designed to be more appealing to the scientific community,  and can now be easily implemented as an independant library for use by other software. 

# Acknowledgements

The first version of the library was initially designed thanks to several projects funded under the HECToR Distributed Computational Science and Engineering (CSE) Service operated by NAG Ltd. The new library has been designed thanks to the support of EPSRC via the CCP Turbulence (EP/T026170/1). 

# References

