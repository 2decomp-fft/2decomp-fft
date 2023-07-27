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
 - name: STFC Daresbury Laboratory,Scientific Computing Department, UKRI, UK 
   index: 1
 - name: Universit\'e de Poitiers, CNRS,ISAE-ENSMA, CNRS, PPRIME, Poitiers, France 
   index: 2
 - name: The University of Edinburgh, EPCC, Edinburgh, UK
   index: 3
 - name: NVIDIA Corporation, Cambridge, UK
   index: 4
 - name: Imperial College London, Department of Aeronautics, London, UK
   index: 5
date: 31 July 2023
bibliography: paper.bib

---

# Summary

The 2DECOMP&FFT library is a software framework written in modern Fortran 
to build large-scale parallel applications. 
It is designed for applications using three-dimensional structured mesh 
and spatially implicit numerical algorithms. At the foundation it implements 
a general-purpose 2D pencil decomposition for data distribution. 
On top it provides a highly scalable and efficient interface to 
perform three-dimensional distributed Fast Fourier Transforms. 


# Statement of need

The 2DECOMP&FFT library [@li_2010_2decomp] is a software framework written 
in Fortran targeting large-scale 
parallel applications. 
The library is based on a 2D-pencil decomposition for data distribution 
on distributed memory systems. The library is at the core of many CFD solver such as 
Xcompact3d [@BARTHOLOMEW2020100550], 
where it has been shown to scale up to hundreds of thousands of cores. 
2DECOMP&FFT mainly relies on MPI, but it offers a user-friendly 
interface that hides the complexity of the communication. 
Indeed, the library offers also a very efficient and flexible 3D tool 
to perform Fast Fourier Transform (FFT) that also 
exposes other popular FFT libraries like FFTW. 

# Acknowledgements

The library was initially designed thanks to several projects funded under the 
HECToR Distributed Computational Science and Engineering (CSE) Service operated by NAG Ltd. 
The new library has been designed thanks to the support of EPSRC via the CCP Turbulence (EP/T026170/1).

# References

