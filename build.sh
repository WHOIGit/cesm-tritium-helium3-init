#!/bin/bash
gfortran -I/usr/include -L/usr/lib -fconvert=big-endian -o tritium_dep_test tritium_dep_test.f90 -lnetcdff -lnetcdf
