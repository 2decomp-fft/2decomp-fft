test2d, timing2d_real, timing2d_complex 
------

The program test2d is to validate the 2D pencil decomposition library. It transposes
a set of data into different storage formats so that all communication routines
are tested.
The programs timing2d_real and timing2d_complex perform the same tests several times 
and are used to report timing and performances for real and complex transforms. 

What to input: The program takes max 6 inputs as : 
               1) p_row [optional]
	       2) p_col [optional] 
	       3) nx    [optional]
	       4) ny    [optional]
	       5) nz    [optional]
	       6) nt    [optional]
	       In case the decomposition is imposed both (1) and (2) are necessary. 
	       If also the resolution is imposed (1-5) are necessary

What to expect: 
- For test2d the output is e success message or an error message with the direction of error in swapping. 
- For the timing, beside the error/success message timing for all transpostions
  (X->Y, Y->Z, Z->Y, Y->X) together with the sum are reported. 
