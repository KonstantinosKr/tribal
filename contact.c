#include "contact.h"
#include "stdio.h"
#include "stdlib.h"
#include "algo.h"
#include "math.h"


void contact_distance (unsigned int nt, REAL *t[3][3], REAL *p[3], REAL *q[3], REAL **distance)
{
    REAL **a;
    a = t[0];
    
    REAL **b;
    b = t[1];
    
    REAL **c;
    c = t[2];
    
    
    REAL **d;
    d = t[0];
    
    REAL **e;
    e = t[1];
    
    REAL **f;
    f = t[2];
    
    for(unsigned int i =0; i<nt;i++)
    {
        printf("A0:%f, A1:%f, A2:%f\n", t[0][0][i], t[0][1][i], t[0][2][i]);
        printf("alpha0:%f, alpha1:%f, alpha2:%f\n", a[0][i], a[1][i], a[2][i]);
    }
    
    REAL dist[1000];
    ispc_bf (nt, a, b, c, d, e, f, p, q, distance[0]);
    for(unsigned int i=0; i<nt;i++)
    {
        printf("distance: %f, p: %f\n",distance[0][i], p[i]);
    }
}

