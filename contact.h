#ifndef HEADER_FILE
#define HEADER_FILE

#include "bf_ispc.h"

/* calculate distances */
void contact_distance (unsigned int nt, REAL *t[3][3], REAL *p[3], REAL *q[3], REAL **distance);

#endif
