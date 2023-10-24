# How To use 2DECOMP&FFT
Full informations on how to use the 2DECOMP&FFT library together with some related theory are available [here](https://2decomp-fft.github.io/).

The 2D Pencil Decomposition API is defined in three Fortran module which should be used by applications as:
```
  use decomp_2d_constants
  use decomp_2d_mpi
  use decomp_2d
```
where ``use decomp_2d_constants`` defines all the parameters, ``use decomp_2d_mpi`` introduces all the MPI 
related interfaces and ``use decomp_2d`` cointains the main decomposition and transposition APIs. 
```
  call decomp_2d_init(nx, ny, nz, p_row, p_col)
```
where ``nx``, ``ny`` and ``nz`` are the size of 3D global data to be distributed over
a 2D processor grid :math:`p_row \times p_col`.
Note that none of the dimensions need to be divisible by ``p_row`` or ``p_col``, i.e. the library can handle non-evenly distributed data.
In case of ``p_row=p_col=0`` an automatic decomposition is selected among all possible combination available.
A key element of this library is a set of communication routines that actually perform the data transpositions.
As mentioned, one needs to perform 4 global transpositions to go through all 3 pencil orientations.
Correspondingly, the library provides 4 communication subroutines:
```
  call transpose_x_to_y(var_in,var_out)
  call transpose_y_to_z(var_in,var_out)
  call transpose_z_to_y(var_in,var_out)
  call transpose_y_to_x(var_in,var_out)
```
The input array ``var_in`` and ``var_output`` array out should have been defined
and contain distributed data for the correct pencil orientations.

Note that the library is written using Fortran's generic interface so different data types are supported
without user intervention. That means in and out above can be either real arrays or complex arrays,
the latter being useful for FFT-type of applications.
Finally, before exit, applications should clean up the memory by:
```
  call decomp_2d_finalize
```
Detailed information about the decomposition API are available [here](https://2decomp-fft.github.io/pages/api_domain.html) 
### Use of the FFT module
To use the FFT programming interface, first of all, one additional Fortran module has to be used:
```
  use decomp_2d_fft
```
Then one needs to initialise the FFT interface by:
```
  call decomp_2d_fft_init(pencil, n1, n2, n3)
```
where ``pencil=PHYSICAL_IN_X`` or ``PHYSICAL_IN_Z`` and ``n1, n2, n3`` is an arbitrary problem size
that can be different from :math:`nx\times ny\times nz`.
For complex-to-complex (c2c) FFTs, the user interface is:
```
  call decomp_2d_fft_3d(input, output, direction)
```
where ``direction`` can be either ``DECOMP_2D_FFT_FORWARD == -1`` for forward transforms, 
or ``DECOMP_2D_FFT_BACKWARD == 1`` for backward transforms.
The input array ``input`` and ``output`` array out are both complex
and have to be either a X-pencil/Z-pencil combination or vice-versa.
The interface for the the real-to-complex and complex-to-real transform is
```
  call decomp_2d_fft_3d(input, output)
```
If the ``input`` data are real type a forward transform is assumed obtaining a complex ``output``.
Similarly a backward FFT is computed if ``input`` is a complex array and ``output`` a real array.
Finally, to release the memory used by the FFT interface:
```
  call decomp_2d_fft_finalize
```
Detailed information about the FFT API are available [here](https://2decomp-fft.github.io/pages/api_fft.html) 
### Use of the IO module
All the I/O functions have been packed in a Fortran module:
```
  use decomp_2d_io
```
To write a single three-dimensional array to a file
```
  call decomp_2d_write_one(ipencil,var,directory,filename,icoarse,io_name)
```
where ``ipencil`` describes how the data is distributed (valid values are: 1 for X-pencil; 2 for
Y-pencil and 3 for Z-pencil); ``var`` is the reference to the data array, which can be either real or
complex; ``directory`` is the path to where I/O should be written; ``filename`` is the name of the
file to be written; ``icoarse`` indicates whether the I/O should be coarsend (valid values are: 0
for no; 1 for the ``nstat`` and 2 for the ``nvisu`` coarsenings); ``io_name`` is the name of the I/O
group to be used. 

There are two ways of writing multiple variables to a single file which may
be used for check-pointing purposes, for example. The newer interface is described first and allows
codes to use the ADIOS2 and MPI-IO backends, the older interface is supported for backwards
compatibility.

When ``decomp_2d_write_one`` is called, the ``directory`` and ``io_name`` are combined to check
whether a particular output location is already opened, if not then a new file will be opened and
written to - this is the "standard" use.  If, however, a file is opened first then the call to
``decomp_2d_write_one`` will append to the current file, resulting in a single file with multiple
fields.  Once the check-pointing is complete the file can then be closed.

The original interface for writing multiple variables to a file and is only
supported by the MPI-IO backend and takes the following form:
```
   call decomp_2d_write_var(fh,disp,ipencil,var)
```
where ``fh`` is a MPI-IO file handle provided by the application (file opened using ``MPI_FILE_OPEN``);
``ipencil`` describes the distribution of the input data (valid values are: 1 for X-pencil; 2 for
Y-pencil and 3 for Z-pencil); ``disp`` (meaning displacement) is a variable of ``kind MPI_OFFSET_KIND``
and of ``intent INOUT``. 
Detailed information about the IO API are available [here](https://2decomp-fft.github.io/pages/api_io.html) 
### Use of the halo exchange
The halo-cell support API provides data structures and nearest-neighbour communication routines 
that support explicit message passing between neighbouring pencils. 
```
  call update_halo(var, var_halo, level)
```
Here the first parameter ``var``, a 3D input array, contains the normal pencil-distributed data as defined by the decomposition library. 
After the subroutine call, the second parameter ``var_halo`` returns all original data plus halo data from the neighbouring processes.
The third parameter level defines how many layers of overlapping are required. 
Detailed information about the halo module are available [here](https://2decomp-fft.github.io/pages/api_halo.html) 
