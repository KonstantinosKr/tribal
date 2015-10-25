#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "error.h"

void load_points(int ptype, unsigned int *nt, unsigned int bodyID, unsigned int startIDX, iREAL *t[3][3], unsigned int tid[], unsigned int pid[], iREAL *mint, iREAL *maxt);

void load_enviroment(int ptype[], unsigned int *nt, unsigned int nParticles, iREAL *t[3][3], unsigned int tid[], unsigned int pid[], iREAL *mint, iREAL *maxt);
void normalize(unsigned int nt, iREAL *t[3][3], iREAL mint, iREAL maxt); 
