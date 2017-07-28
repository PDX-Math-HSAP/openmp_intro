#! /bin/bash

gcc impls.c -o impls_release -O3 -lm -fopenmp
gcc impls.c -o impls_debug   -O0 -lm -fopenmp

rm *.res
./impls_release > release_results.res
./impls_debug > debug_results.res
