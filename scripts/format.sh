#!/usr/bin/env sh
#
# format.sh
#
# Runs the prettifier over 2decomp sources.

# Find 2decomp root
SCRDIR=$( dirname -- "$0"; )
D2DDIR=${SCRDIR}/../

# fprettify options
# --enable-replacements : Replaces relational operators, e.g. .lt. -> <
# --c-relations : Use C-style relational operators, e.g. <=
# -w 3 : Whitespace preset - operators, print/read, *, /, +, -
fprettify --enable-replacements --c-relations -w 3 ${D2DDIR}/src/*f90 ${D2DDIR}/examples/*/*f90
