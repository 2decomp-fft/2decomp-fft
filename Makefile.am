if HAVE_FFTW3
AM_FCFLAGS=@FFTW3FCFLAGS@
endif
if HAVE_MKL
AM_FCFLAGS=@MKLFCFLAGS@
endif

# Decomp2d library
lib_LTLIBRARIES = libdecomp2d.la
libdecomp2d_la_SOURCES = \
decomp_2d.f90 \
glassman.f90 \
fft_@FFTENGINE@.f90
libdecomp2d_la_LDFLAGS = -no-undefined
if HAVE_MKL
libdecomp2d_la_SOURCES += mkl_dfti.f90
endif

# Allow parallel make
glassman.lo: decomp_2d.lo
if HAVE_MKL
fft_@FFTENGINE@.lo: glassman.lo mkl_dfti.lo
else
fft_@FFTENGINE@.lo: glassman.lo
endif

clean-local:
	-rm -f *.mod
