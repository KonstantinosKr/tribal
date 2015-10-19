#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "error.h"

unsigned long long int load_points(unsigned int ptype, unsigned int bodyID, unsigned int startIDX, iREAL *t[3][3], unsigned long long int tid[], unsigned long long int pid[], iREAL *mint, iREAL *maxt);

unsigned long long int load_enviroment(unsigned int ptype[], unsigned int nParticles, iREAL *t[3][3], unsigned long long int tid[], unsigned long long int pid[], iREAL *mint, iREAL *maxt);
void normalize(unsigned long long int nt, iREAL *t[3][3], iREAL mint, iREAL maxt); 
