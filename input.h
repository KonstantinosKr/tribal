#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "error.h"

unsigned int load_pointsVTK(iREAL *t[3][3], unsigned int tid[], iREAL *mint, iREAL *maxt);

void normalize(unsigned int nt, iREAL *t[3][3], iREAL mint, iREAL maxt); 
