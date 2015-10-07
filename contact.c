#include "contact.h"
#include "stdio.h"
#include "stdlib.h"
#include "algo.h"
#include "math.h"

void contact_distance (unsigned int nt, iREAL *t[3][3], iREAL *p[3], iREAL *q[3], iREAL *distance)
{
    iREAL *d[3], *e[3], *f[3];
    iREAL *pp[3], *qq[3];
    
    //allocate memory
    for(int i=0;i<3;i++)
    {
      d[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      e[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      f[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      pp[i] = (iREAL *) malloc(nt * sizeof(iREAL));
      qq[i] = (iREAL *) malloc(nt * sizeof(iREAL));
    }
   
    //Set triangle 2 points D,E,F
    for(unsigned int i=0;i<nt;i++)
    {
      d[0][i] = t[0][0][i];
      d[1][i] = t[0][1][i];
      d[2][i] = t[0][2][i];
      
      e[0][i] = t[1][0][i];
      e[1][i] = t[1][1][i];
      e[2][i] = t[1][2][i];
      
      f[0][i] = t[2][0][i];
      f[1][i] = t[2][1][i];
      f[2][i] = t[2][2][i];
    }
    
    iREAL a[3], b[3], c[3];
    //Set triangle 1 points A,B,C
    for(unsigned int i=0;i<nt;i++)
    { 
      a[0] = t[0][0][i];
      a[1] = t[0][1][i];
      a[2] = t[0][2][i];
      
      b[0] = t[1][0][i];
      b[1] = t[1][1][i];
      b[2] = t[1][2][i];
      
      c[0] = t[2][0][i];
      c[1] = t[2][1][i];
      c[2] = t[2][2][i];

      ispc_bf (nt, a, b, c, d, e, f, pp, qq); 
        
      /*
      for(unsigned int j=0;j<nt;j++)
      {
        p[0][(i*nt)+j]= pp[0][j]; 
        p[1][(i*nt)+j]= pp[1][j]; 
        p[2][(i*nt)+j]= pp[2][j]; 
          
        q[0][(i*nt)+j] = qq[0][j];
        q[1][(i*nt)+j] = qq[1][j];
        q[2][(i*nt)+j] = qq[2][j];

        distance[(i*nt)+j] = sqrt(pow((qq[0][j]-pp[0][j]),2)+pow((qq[1][j]-pp[1][j]),2)+pow((qq[2][j]-pp[1][j]),2));
      }*/
    }
}

