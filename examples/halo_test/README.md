# Halo test

List of the tests:
- [HALO test](halo_test.f90): Test for the halo exchange capability of the library. 

This example demonstrates the use of the halo-cell support API. It calculates
the divergence of an arbitrary field, which contains evaluation of spatial
derivatives in all three dimensions. The calculation was first implemented via
the global transposition routines, then via halo-cell exchanges. Identical
results are to be expected regardless of the communication algorithm. The 
computation is based on an explicit finite difference method so clearly using 
the halo-cell support API is more efficient.

What to input: The program takes max 5 inputs as: 

1. p_row [optional]
1. p_col [optional] 
1. nx    [optional]
1. ny    [optional]
1. nz    [optional]

If the decomposition is imposed both (1) and (2) are necessary. 
If the resolution is imposed (1-5) are necessary.

What to expect: the output using different communication algorithms should be the same up to machine precision. 
