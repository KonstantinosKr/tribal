#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include "loba.h"
#include "error.h"
#include "tribal_ispc.h"
#include "contact.h"

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

void migrate_status(int myrank, unsigned int nt, REAL *t[3][3], unsigned int timesteps, unsigned int *importn, unsigned int *exportn)
{   
  //sort 

  if(myrank != 0)
  {
    //MPI_send();
    //MPI_send();
    //MPI_send();
  
  
  
  } else if(myrank == 0)
  {
    //MPI_receive();
    //MPI_receive();

    unsigned int **nparticles, **import, **export;
    ERRMEM (nparticle = (unsigned int **) malloc(timesteps * sizeof(unsigned int *)));
    ERRMEM (import = (unsigned int **) malloc(timesteps * sizeof(unsigned int *)));
    ERRMEM (export = (unsigned int **) malloc(timesteps * sizeof(unsigned int *)));
    
    for(unsigned int i=0;i<nt;i++)
    {
      ERRMEM (nparticles[i] = (unsigned int *) malloc(nt * sizeof(REAL))); 
      ERRMEM (import[i] = (unsigned int *) malloc(nt * sizeof(REAL)));
      ERRMEM (export[i] = (unsigned int *) malloc(nt * sizeof(REAL)));
    }

    int nnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
  
    FILE *fp = fopen("mpi.csv", "w+");
    
    for(int i=0;i<nnodes;i++)
    {
      fprintf(fp,"npar_node[%u], imp_dep[%u], exp_dep[%u]", i,i,i);
      if(i!=nnodes-1)
      {
        fprintf(fp,", ");
      }
    }
    fprintf(fp,"\n");
  
    //replace with time/step for all timesteps
    for(unsigned int i=0;i<timesteps;i++)
    { 
      for(int j=0;j<nnodes;i++)
      {
        fprintf(fp,"%u, %u, %u", npar[j][i], import_dependencies[j][i], export_dependencies[j][i]);
      }
      if(i!=timesteps-1)
      {
        fprintf(fp,", ");
      }else
      {
        fprintf(fp,"\n");
      }

    } 
    fclose(fp);
  }
}

/* migrate triangles "in-place" to new ranks */
static void migrate_triangles (int myrank, unsigned int size, unsigned int nt, REAL *t[3][3], REAL *v[3], int *rank, unsigned int *importn, unsigned int *exportn)
{
  *importn=0;
  *exportn=0;
  
  int nnodes;
  MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
  
  unsigned int **send_idx, *pivot; REAL **tbuffer[3][3], **vbuffer[3];

  ERRMEM (tbuffer = (REAL **) malloc(nnodes*sizeof(REAL*)));
  ERRMEM (vbuffer = (REAL **) malloc(nnodes*sizeof(REAL*)));
  ERRMEM (send_idx = (unsigned int **) malloc(nnodes*sizeof(unsigned int*)));
  ERRMEM (pivot = (unsigned int *) malloc(nnodes*sizeof(unsigned int)));
  for(int i=0;i<nnodes;i++)
  {
    ERRMEM (tbuffer[i] = (REAL *) malloc(nt*sizeof(REAL*)));
    ERRMEM (vbuffer[i] = (REAL *) malloc(nt*sizeof(REAL*)));
    pivot[i] = 0;
    ERRMEM (send_idx[i] = (unsigned int *) malloc(nt*sizeof(unsigned int)));
  }
  
  for (unsigned int i = 0; i < nt; i++) 
  {
    send_idx[pivot[rank[i]]++][rank[i]] = i; 
  }

  for(int i=0;i<nnodes;i++)
  {
    for(unsigned int j=0;j<pivot[i];j++)
    {
      for(int k=0;k<3;k++)
      {
        tbuffer[0][k][send_idx[j][i]][i] = t[0][k][send_idx[j][i]]; 
        tbuffer[1][k][send_idx[j][i]][i] = t[1][k][send_idx[j][i]];
        tbuffer[2][k][send_idx[j][i]][i] = t[2][k][send_idx[j][i]];

        vbuffer[k][send_idx[j][i]][i] = v[k][send_idx[j][i]];
      }
    }
  
    //store number of exports
    if(i != myrank)
    {
      *exportn += pivot[i];
    }
  }//can merge to next loop if correct since <nnodes limit;

  nt = nt - *exportn;

  unsigned int *pivot_buffer;
  ERRMEM (pivot_buffer = (unsigned int *) malloc(nnodes*sizeof(unsigned int)));
  pivot_buffer = pivot;

  MPI_Status status;

  //send buffers
  for(int i=0;i<nnodes;i++)
  {
    //MPI_send(&pivot[i], pivot[i], MPI_INT, i, myrank, MPI_COMM_WORLD);
    //MPI_send(&tbuffer[i], pivot[i], MPI_DOUBLE, i, myrank, MPI_COMM_WORLD);
    //MPI_send(&vbuffer[i], pivot[i], MPI_DOUBLE, i, myrank, MPI_COMM_WORLD);

    MPI_Sendrecv(&pivot_buffer[i], pivot[i], MPI_INT, i, myrank, MPI_COMM_WORLD, &pivot_buffer[i], pivot[i], MPI_INT, i, myrank, MPI_COMM_WORLD, &status);
    
    *importn += pivot[myrank];
    
    MPI_Sendrecv(&tbuffer[i], pivot[i], MPI_DOUBLE, i, myrank, MPI_COMM_WORLD, &tbuffer[i], pivot[i], MPI_DOUBLE, i, myrank, MPI_COMM_WORLD, &status);
    MPI_Sendrecv(&vbuffer[i], pivot[i], MPI_DOUBLE, i, myrank, MPI_COMM_WORLD, &vbuffer[i], pivot[i], MPI_DOUBLE, i, myrank, MPI_COMM_WORLD, &status);  
  }
  nt = nt + *importn;  
}

int main (int argc, char **argv)
{
  REAL *t[3][3]; /* triangles */
  REAL *v[3]; /* velocities */
  REAL **d; /*distance */
  REAL *p[3],*q[3];
  unsigned int *tid; /* triangle identifiers */
  REAL lo[3] = {0, 0, 0}; /* lower corner */
  REAL hi[3] = {1, 1, 1}; /* upper corner */
  unsigned int nt; /* number of triangles */
  int *rank; /* migration ranks */
  unsigned int *pid; /*particle identifier */
  unsigned int size; /* buffers size */
  int myrank;
  
  /*load balancing n dependencies */
  unsigned int *importn, *exportn;
  
  /* init */ 
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);

  if (myrank == 0)
  {
    /* set nt */
    if (argc > 1) nt = atoi (argv[1]);
    else nt = 10;

    /* buffers */
    size = 4*nt;

    for (int i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = (REAL *) malloc (sizeof(REAL[size])));
    
      ERRMEM (p[i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (q[i] = (REAL *) malloc (sizeof(REAL[size])));
    }
    ERRMEM (d = (REAL **) malloc (nt * sizeof(REAL *)));
    for(unsigned int j=0;j<nt;j++)
    {
      ERRMEM (d[j] = (REAL *) malloc (nt * sizeof(REAL)));
    }
    ERRMEM (rank = (int *) malloc (sizeof(int[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (importn = (unsigned int *) malloc(sizeof(unsigned int[size])));
    ERRMEM (exportn = (unsigned int *) malloc(sizeof(unsigned int[size])));
    //nt = load_pointsVTK(t, tid); 
    /* generate triangles and velocities */
    generate_triangles_and_velocities (lo, hi, nt, t, v, tid); 
  }
  else
  {
    /* set nt */
    nt = 0;

    /* buffers */
    if (argc > 1) size = atoi (argv[1])*4;
    else size = 10*4;

    for (int i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = (REAL *) malloc (sizeof(REAL[size])));
      
      ERRMEM (p[i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (q[i] = (REAL *) malloc (sizeof(REAL[size])));
    }
    ERRMEM (d = (REAL **) malloc (nt * sizeof(REAL *)));
    for(unsigned int j=0;j<nt;j++)
    {
      ERRMEM (d[j] = (REAL *) malloc (nt * sizeof(REAL)));
    }
    ERRMEM (importn = (unsigned int *) malloc(sizeof(unsigned int[size])));
    ERRMEM (exportn = (unsigned int *) malloc(sizeof(unsigned int[size])));
    ERRMEM (rank = (int *) malloc (sizeof(int[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
  }
  
  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  REAL step = 1E-3, time; unsigned int timesteps=0;
  
  //for (time = 0.0; time < 1.0; time += step)
  {
    loba_balance (lb, nt, t[0], tid, 1.1, rank);
    
    for (unsigned int j = 0; j < nt; j++) 
    {
      printf ("rank[%u] = %d\n", j, rank[j]);
    }
    
    migrate_triangles (myrank, size, nt, t, v, rank, &importn[timesteps], &exporn[timesteps]);

    contact_distance(nt, t, p, q, d); 

    //integrate_triangles (step, lo, hi, nt, t, v);
    timesteps++;
  }
  
  if(myrank == 0)
  {
    migrate_status(myrank, nt, t, timesteps, import, export);
    //write_pointsVTK(nt, t, v);
  }

  /* finalise */
  loba_destroy (lb);

  for (int i = 0; i < 3; i ++)
  {
    free (t[0][i]);
    free (t[1][i]);
    free (t[2][i]);
    free (v[i]);
    free (p[i]);
    free (q[i]);
  }
  free (rank);

  MPI_Finalize ();

  return 0;
}
