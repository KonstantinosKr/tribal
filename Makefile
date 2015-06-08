# C files
C_SRC=loba.c tribal.c

# ISPC files
ISPC_SRC=tribal.ispc

# ISPC targets
ISPC_TARGETS=avx2

# MPI
MPICC=mpicc
MPICXX=mpic++

# Zoltan paths
ZOLTANINC = -I/usr/local/include
ZOLTANLIB = -L/usr/local/lib -lzoltan

# Program name
EXE=tribal

# Floating point type
REAL=double

# Debug version
DEBUG=yes

# Do the rest
include common.mk
