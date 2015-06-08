#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include "loba.h"
#include "error.h"
#include "tribal_ispc.h"

/* migrate triangles "in-place" to new ranks */
static void migrate_triangles (int size, int nt, REAL *t[3][3], REAL *v[3], int *rank)
{
  /* TODO */
}

unsigned int load_pointsVTK(double *t[3][3], unsigned int tid[])
{
    FILE *fp1 = fopen("input1.vtk", "r");
    
    unsigned int nt=0;
    
    if( fp1 == NULL )
    {
        perror("Error while opening the file.\n");
        exit(EXIT_FAILURE);
    }
    
    char ch, word[100];
    double *point[3];
    
    do {
        ch = fscanf(fp1,"%s",word);
        if(strcmp(word, "POINTS")==0)
        {
            printf("found!");
            ch = fscanf(fp1,"%s",word);
            unsigned int n = atoi(word);
            //get points
            ch = fscanf(fp1,"%s",word);
            
            ERRMEM (point[0] = (double *)malloc (sizeof(double[n])));
            ERRMEM (point[1] = (double *)malloc (sizeof(double[n])));
            ERRMEM (point[2] = (double *)malloc (sizeof(double[n])));
            
            for(unsigned int i=0;i<n;i++)
            {
                ch = fscanf(fp1, "%s", word);
                point[0][i] = atof(word);
                ch = fscanf(fp1, "%s", word);
                point[1][i] = atof(word);
                ch = fscanf(fp1, "%s", word);
                point[2][i] = atof(word);
            }
        }
        if(strcmp(word, "CELLS")==0)
        { 
            ch = fscanf(fp1,"%s",word);
            unsigned int n = atoi(word);
            nt = n;
            ch = fscanf(fp1,"%s",word);
            printf(":::%u::\n",n);
            for(unsigned int i=0;i<n;i++)
            {
                ch = fscanf(fp1,"%s",word);
                ch = fscanf(fp1,"%s",word);
                
                unsigned int index = atoi(word);
                t[0][0][i] = point[0][index];
                t[0][1][i] = point[1][index];
                t[0][2][i] = point[2][index];
                
                ch = fscanf(fp1,"%s",word);
                index = atoi(word);
                t[1][0][i] = point[0][index];
                t[1][1][i] = point[1][index];
                t[1][2][i] = point[2][index];
                
                index = atoi(word);
                t[2][0][i] = point[0][index];
                t[2][1][i] = point[1][index];
                t[2][2][i] = point[2][index];
            }
        }
    } while (ch != EOF);
    return nt;
}


void write_pointsVTK(unsigned int nt, REAL *t[3][3], REAL *v[3])
{
    FILE *fp = fopen("output.vtk", "w+");
    
    fprintf(fp,"# vtk DataFile Version 2.0\nOutput vtk file\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS %i float\n", nt*3);
    
    for(unsigned i = 0; i < nt; i++)
    {
        fprintf(fp,"%.5f %.5f %.5f\n%.5f %.5f %.5f\n%.5f %.5f %.5f\n", t[0][0][i], t[0][1][i], t[0][2][i], t[1][0][i], t[1][1][i], t[1][2][i], t[2][0][i], t[2][1][i], t[2][2][i]);
    }
    fprintf(fp,"\nCELLS %i %i\n", nt, nt+nt*3);
    for(unsigned i = 0; i < nt*3; i=i+3)
    {
        fprintf(fp,"3 %i %i %i\n", i, i+1, i+2);
    }
    
    fprintf(fp,"\nCELL_TYPES %i\n", nt);
    for(unsigned i = 0; i<nt; i++)
    {
        fprintf(fp,"5\n");
    }
    fclose(fp);
}

int main (int argc, char **argv)
{
  REAL *t[3][3]; /* triangles */
  REAL *v[3]; /* velocities */
  unsigned int *tid; /* triangle identifiers */
  REAL lo[3] = {0, 0, 0}; /* lower corner */
  REAL hi[3] = {1, 1, 1}; /* upper corner */
  unsigned int nt; /* number of triangles */
  int *rank; /* migration ranks */
  unsigned int size; /* buffers size */
  int i;
  int myrank;

  /* init */ 
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
  {
    /* set nt */
    if (argc > 1) nt = atoi (argv[1]);
    else nt = 10000;

    /* buffers */
    size = 4*nt;

#if 1
    for (i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = (REAL *) malloc (sizeof(REAL[size])));
    }
    ERRMEM (rank = (int *) malloc (sizeof(int[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));

//    nt = load_pointsVTK(t, tid); 
//    write_pointsVTK(nt, t, v);
    
    /* generate triangles and velocities */
    generate_triangles_and_velocities (lo, hi, nt, t, v, tid); 
#endif
  }
  else
  {
    /* set nt */
    nt = 0;

    /* buffers */
    if (argc > 1) size = atoi (argv[1])*4;
    else size = 10000*4;

#if 1
    for (i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = (REAL *) malloc (sizeof(REAL[size])));
    }
    ERRMEM (rank = (int *) malloc (sizeof(int[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    
    /* generate triangles and velocities */  
    generate_triangles_and_velocities (lo, hi, nt, t, v, tid);
#endif
  }
  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  REAL step = 1E-3, time;
  
  //for (time = 0.0; time < 1.0; time += step)
  {
    loba_balance (lb, nt, t[0], tid, 1.1, rank);

    //migrate_triangles (size, nt, t, v, rank);

    //integrate_triangles (step, lo, hi, nt, t, v);
  }

  //write_pointsVTK(nt, t, v);
  /* finalise */
  loba_destroy (lb);

  for (i = 0; i < 3; i ++)
  {
    free (t[0][i]);
    free (t[1][i]);
    free (t[2][i]);
    free (v[i]);
  }
  free (rank);

  MPI_Finalize ();

  return 0;
}
