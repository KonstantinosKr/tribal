#include "contact.h"
#include "stdio.h"
#include "stdlib.h"
#include "algo.h"
#include "math.h"


void contact_distance (unsigned int nt, REAL *t[3][3], REAL *p[3], REAL *q[3], REAL *distance)
{

    REAL *d[3], *e[3], *f[3];

    for(int i=0;i<3;i++)
    {
      d[i] = (REAL *) malloc(nt * sizeof(REAL));
      e[i] = (REAL *) malloc(nt * sizeof(REAL));
      f[i] = (REAL *) malloc(nt * sizeof(REAL));
    }
    
    for(unsigned int i=0;i<nt;i++)
    {
      for(int z=0;z<3;z++)
      {
        d[z][i] = t[0][z][i];
        e[z][i] = t[1][z][i];
        f[z][i] = t[2][z][i];
      }
    }
    
    REAL a[3], b[3], c[3];

    for(unsigned int i=0;i<nt;i++)
    { 
      for(int z=0;z<3;z++)
      {
        a[z] = t[0][z][i];
        b[z] = t[1][z][i];
        c[z] = t[2][z][i];

        ispc_bf (nt, a, b, c, d, e, f, p, q); 
      }
    }
/*
    unsigned int counter=0;
    for(unsigned int i=0;i<nt;i++)
    {
      for(unsigned int j=0;j<nt;j++)
      {
        printf("counter:%u, distance: %f, %f\n",counter, distance[counter], p[0][counter]);
        counter++;
      }
    }*/
}

