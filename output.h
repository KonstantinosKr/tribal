#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "error.h"

void write_pointsVTK(int myrank, unsigned int nt, iREAL *t[3][3], iREAL *v[3], iREAL lo[3], iREAL hi[3], unsigned int timesteps);

void postProcessing(int nranks, unsigned int long long size, unsigned int timesteps);
