#ifndef HEADER_FILE
#define HEADER_FILE

#include "bf_ispc.h"

/* calculate distances */
unsigned int contact_distance (unsigned int start, unsigned int end, unsigned long long int size, iREAL *t[3][3], iREAL *p[3], iREAL *q[3], iREAL *distance);

#endif
