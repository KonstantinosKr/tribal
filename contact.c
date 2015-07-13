#include "contact.h"
#include "stdio.h"
#include "stdlib.h"
#include "algo.h"
#include "math.h"


void contact_distance (unsigned int nt, REAL *t[3][3], REAL *p[3], REAL *q[3], REAL *distance)
{
    /*REAL **a, **b, **c;
    a = t[0];
    b = t[1];
    c = t[2];
    */

    REAL *a[3], *b[3], *c[3];
    for(int i=0;i<3;i++)
    {
      a[i] = (REAL *) malloc(nt * nt * nt * sizeof(REAL));
      b[i] = (REAL *) malloc(nt * nt * nt * sizeof(REAL));
      c[i] = (REAL *) malloc(nt * nt * nt * sizeof(REAL));
    }

    for(unsigned int i=0;i<nt;i++)
    { 
      for(unsigned int x=0;x<nt;x++)
      {
        for(int z=0;z<3;z++)
        {
          a[z][i*nt+x] = t[0][z][i];
          b[z][i*nt+x] = t[1][z][i];
          c[z][i*nt+x] = t[2][z][i];
        }
      }
      printf("alpha0:%f, alpha1:%f, alpha2:%f\n", a[0][i*nt], a[1][i*nt], a[2][i*nt]); 
      printf("beta0:%f, beta1:%f, beta2:%f\n", b[0][i*nt], b[1][i*nt], b[2][i*nt]);
      printf("ceta0:%f, ceta1:%f, ceta2:%f\n", c[0][i*nt], c[1][i*nt], c[2][i*nt]); 
    }
    printf("-----------\n\n");

    REAL *d[3], *e[3], *f[3];

    for(int i=0;i<3;i++)
    {
      d[i] = (REAL *) malloc(nt * nt * sizeof(REAL));
      e[i] = (REAL *) malloc(nt * nt * sizeof(REAL));
      f[i] = (REAL *) malloc(nt * nt * sizeof(REAL));
    }
    
    for(unsigned int i=0;i<nt;i++)
    {
      for(unsigned int x=0;x<nt;x++)
      { 
        for(int z=0;z<3;z++)
        {
          d[z][i*nt+x] = t[0][z][x];
          e[z][i*nt+x] = t[1][z][x];
          f[z][i*nt+x] = t[2][z][x];
          //printf("fvalue:%f\n", f[z][x]);
        }
      }
      
      printf("delta0:%f, delta1:%f, delta2:%f\n", d[0][i*nt], d[1][i*nt], d[2][i*nt]); 
      printf("epsilon0:%f, epsilon1:%f, epsilon2:%f\n", e[0][i*nt], e[1][i*nt], e[2][i*nt]);
      printf("feta0:%f, feta1:%f, feta2:%f\n", f[0][i*nt], f[1][i*nt], f[2][i*nt]); 
    }
    
    ispc_bf (nt*nt, a, b, c, d, e, f, p, q, distance); 

    unsigned int counter=0;
    for(unsigned int i=0;i<nt;i++)
    {
      for(unsigned int j=0;j<nt;j++)
      {
        printf("counter:%u, distance: %f, %f\n",counter, distance[counter], p[0][counter]);
        counter++;
      }
    }
}

