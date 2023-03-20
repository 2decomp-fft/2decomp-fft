#!/usr/bin/env sh
#
# format.sh
#
# Runs the prettifier over 2decomp sources.

# Find 2decomp root
SCRDIR=$( dirname -- "$0"; )
D2DDIR=${SCRDIR}/../

# Run fprettify with config file scripts/.fprettifyrc
fprettify -c ${SCRDIR}/.fprettifyrc ${D2DDIR}/src/*f90 ${D2DDIR}/examples/*/*f90
