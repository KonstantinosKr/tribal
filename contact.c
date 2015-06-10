#include "contact.h"
#include "stdio.h"
#include "stdlib.h"
#include "algo.h"
#include "math.h"


void contact_distance (unsigned int nt, REAL *t[3][3], REAL *p[3], REAL *q[3], REAL **distance)
{
    REAL **a, **b, **c;
    a = t[0];
    b = t[1];
    c = t[2];

    for(unsigned int i=0;i<nt;i++)
    {
      //printf("A0:%f, A1:%f, A2:%f\n", t[0][0][i], t[0][1][i], t[0][2][i]);
      printf("alpha0:%f, alpha1:%f, alpha2:%f\n", a[0][i], a[1][i], a[2][i]); 
      printf("beta0:%f, beta1:%f, beta2:%f\n", b[0][i], b[1][i], b[2][i]);
      printf("ceta0:%f, ceta1:%f, ceta2:%f\n", c[0][i], c[1][i], c[2][i]); 
    }
    printf("-----------\n\n");

    REAL **d, **e, **f;

    d = (REAL **) malloc(nt * sizeof(REAL *));
    e = (REAL **) malloc(nt * sizeof(REAL *));
    f = (REAL **) malloc(nt * sizeof(REAL *));

    for(int i=0;i<3;i++)
    {
      d[i] = (REAL *) malloc(nt * sizeof(REAL));
      e[i] = (REAL *) malloc(nt * sizeof(REAL));
      f[i] = (REAL *) malloc(nt * sizeof(REAL));
    }
    
    for(unsigned int i=0;i<nt;i++)
    {
      
      for(unsigned int x=0;x<nt;x++)
      {
        for(int z=0;z<3;z++)
        {
          d[z][x] = t[0][z][i];
          e[z][x] = t[1][z][i];
          f[z][x] = t[2][z][i];
        }
      }
       
      printf("delta0:%f, delta1:%f, delta2:%f\n", d[0][i], d[1][i], d[2][i]); 
      printf("epsilon0:%f, epsilon1:%f, epsilon2:%f\n", e[0][i], e[1][i], e[2][i]);
      printf("feta0:%f, feta1:%f, feta2:%f\n", f[0][i], f[1][i], f[2][i]); 
      
      ispc_bf (nt, a, b, c, d, e, f, p, q, distance[i]);
      
      for(unsigned int j=0;j<nt;j++)
      {
        printf("distance: %f, p: %f\n",distance[i][j], p[0][j]);
      }
    }
}

