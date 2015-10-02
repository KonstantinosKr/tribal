#include "contact.h"
#include "stdio.h"
#include "stdlib.h"
#include "algo.h"
#include "math.h"


void contact_distance (unsigned int nt, iREAL *t[3][3], iREAL *p[3], iREAL *q[3], iREAL *distance)
{

    iREAL *d[3], *e[3], *f[3];
    iREAL *pp[3], *qq[3];

    for(int i=0;i<3;i++)
    {
      d[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      e[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      f[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      pp[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      qq[i] = (iREAL *) malloc(nt * sizeof(iREAL));
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
    
    iREAL a[3], b[3], c[3];

    for(unsigned int i=0;i<nt;i++)
    { 
      for(int z=0;z<3;z++)
      {
        a[z] = t[0][z][i];
        b[z] = t[1][z][i];
        c[z] = t[2][z][i];
      }

//      ispc_bf (nt, a, b, c, d, e, f, pp, qq); 
        
      for(unsigned int x=0;x<nt;x++)
      {
        for(int z=0;z<3;z++)
        {
          p[z][(i*nt)+x]= pp[z][x]; 
          q[z][(i*nt)+x] = qq[z][x];
        }
        distance[(i*nt)+x] = sqrt(pow((qq[0]-pp[0]),2)+pow((qq[1]-pp[1]),2)+pow((qq[2]-pp[1]),2));
      }
    }
}

