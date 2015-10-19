#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "error.h"
#include "loba.h"
#include <zoltan.h>

void output_state(struct loba *lb, int myrank, unsigned long long int nt, iREAL *t[3][3], iREAL *v[3], unsigned long long int timesteps);

void postProcessing(int nranks, unsigned long long int size, unsigned long long int timesteps);
