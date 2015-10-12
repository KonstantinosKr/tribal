#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "error.h"

unsigned int load_points(unsigned int ptype, unsigned int bodyID, unsigned int startIDX, iREAL *t[3][3], unsigned int tid[], unsigned int pid[], iREAL *mint, iREAL *maxt);

unsigned int load_enviroment(unsigned int ptype[], unsigned int nParticles, iREAL *t[3][3], unsigned int tid[], unsigned int pid[], iREAL *mint, iREAL *maxt);
void normalize(unsigned int nt, iREAL *t[3][3], iREAL mint, iREAL maxt); 
