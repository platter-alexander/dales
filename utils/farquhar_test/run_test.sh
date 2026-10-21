#!/bin/bash
# Build and run the leaf-level Farquhar test against the DALES debug build.
# Usage: ./run_test.sh   (after building DALES in build-debug)
set -e
HERE="$(cd "$(dirname "$0")" && pwd)"
BUILD="${BUILD:-$HERE/../../build-debug}"
OBJ="$BUILD/src/CMakeFiles/dales4.dir"
OUT="${OUT:-$BUILD/test_farquhar_leaf}"   # executable goes into the (git-ignored) build directory

# all DALES objects except the main program
objs=$(ls "$OBJ"/*.o | grep -v '/program.f90.o$')

mpif90 -fdefault-real-8 -ffree-line-length-none -fbounds-check -fbacktrace -O0 -g \
       -ffpe-trap=invalid,zero,overflow -I"$BUILD/src" \
       "$HERE/test_farquhar_leaf.f90" $objs -lnetcdff -lnetcdf -o "$OUT"
"$OUT"
