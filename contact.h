#ifndef HEADER_FILE
#define HEADER_FILE

#include "bf_ispc.h"

/* calculate distances */
void contact_detection (unsigned long long int s1, unsigned long long int e1, unsigned long long int s2, unsigned long long int e2,  unsigned long long int size, iREAL *t[3][3], iREAL *p[3], iREAL *q[3], iREAL *distance);

#endif
