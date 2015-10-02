#ifndef HEADER_FILE
#define HEADER_FILE

#include "bf_ispc.h"

/* calculate distances */
void contact_distance (unsigned int nt, iREAL *t[3][3], iREAL *p[3], iREAL *q[3], iREAL *distance);

#endif
