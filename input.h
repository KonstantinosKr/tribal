#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include "error.h"

void load_points(int ptype, unsigned long long int *nt, unsigned long long int bodyID, unsigned long long int startIDX, iREAL *t[3][3], unsigned long long int tid[], unsigned long long int pid[], iREAL *mint, iREAL *maxt);

void load_enviroment(int ptype[], unsigned long long int *nt, unsigned long long int nParticles, iREAL *t[3][3], unsigned long long int tid[], unsigned long long int pid[], iREAL *mint, iREAL *maxt);
void normalize(unsigned long long int nt, iREAL *t[3][3], iREAL mint, iREAL maxt); 
