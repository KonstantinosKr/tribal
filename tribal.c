#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "input.h"
#include "output.h"
#include "loba.h"
#include "error.h"
#include "tribal_ispc.h"
#include "contact.h"
#include "motion.h"
#include "migration.h"

int main (int argc, char **argv)
{
  iREAL *t[3][3]; /* triangles */
  iREAL *v[3]; /* velocities */
  iREAL *distance; /*distance */
  iREAL *p[3],*q[3];//p and q points
  unsigned long long int nt; /* number of triangles */
  unsigned long long int *pid; /*particle identifier */
  unsigned long long int *tid; /* triangle identifiers */
  iREAL lo[3] = {-255, -255, -255}; /* lower corner */
  iREAL hi[3] = {255, 255, 255}; /* upper corner */
  
  unsigned long long int nParticles = 5;
  unsigned long long int size = 4000000; /* memory buffer size */
  int nprocs;
  int myrank;

  /* init */ 
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   
  if (myrank == 0)
  {
    /* set nt */
    if (argc > 1) nt = atol (argv[1]);
    else nt = 0;

    for (int i = 0; i < 3; i ++)
    { 
      (t[0][i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (t[1][i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (t[2][i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (v[i] = (iREAL *) malloc (size*sizeof(iREAL)));
    
      (p[i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (q[i] = (iREAL *) malloc (size*sizeof(iREAL)));
    }
    (distance = (iREAL *) malloc (size*sizeof(iREAL)));
    (tid = (unsigned long long int *) malloc (size*sizeof(unsigned long long int)));
    (pid = (unsigned long long int *) malloc (size*sizeof(unsigned long long int)));
    
    for(unsigned long long int i=0;i<size;i++) tid[i] = ULLONG_MAX; 
    
    iREAL mint, maxt;
    
    //Input Type
    //0: Triangulated Mesh
    //1: Triangle
    //2: Sphere
    //3: Square
    //4: Hexahedron
  
    int ptype[nParticles];
    ptype[0] = 0;
    ptype[1] = 0;
    ptype[2] = 0;
    ptype[3] = 0;
    ptype[4] = 0;
    ptype[5] = 0;
    ptype[6] = 0;
    ptype[7] = 0;
    ptype[8] = 0;
    //ptype[9] = 0;
    //ptype[10] = 0;

    load_enviroment(ptype, &nt, nParticles, t, tid, pid, &mint, &maxt);
    
    gen_velocities(lo, hi, nt, v);
  }
  else
  {
    /* set nt */
    nt = 0;

    /* buffers */
    if (argc > 1) size = atol (argv[1])*4;

    for (int i = 0; i < 3; i ++)
    {
      (t[0][i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (t[1][i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (t[2][i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (v[i] = (iREAL *) malloc (size*sizeof(iREAL)));
      
      (p[i] = (iREAL *) malloc (size*sizeof(iREAL)));
      (q[i] = (iREAL *) malloc (size*sizeof(iREAL)));
    }
    (distance = (iREAL*) malloc (size*sizeof(iREAL)));
    (tid = (unsigned long long int *) malloc (sizeof(unsigned long long int[size])));
    (pid = (unsigned long long int *) malloc (sizeof(unsigned long long int[size])));
      
    //for(unsigned long long int i=0;i<size;i++) tid[i] = ULLONG_MAX;
  }
  
  int num_import, num_export;
  int *import_procs, *export_procs;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids, export_local_ids;
  

  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  iREAL step = 1E-3, time; unsigned long long int timesteps=0;
  
  //for (time = 0.0; time < 1.0; time += step)
  for(time = 0; time < 100; time++)
  {
    if(myrank == 0){printf("TIMESTEP: %i\n", timesteps);}
    
    loba_balance (lb, nt, t[0], tid, 1.1,
                  &num_import, &import_procs, 
                  &num_export, &export_procs, 
                  &import_global_ids, &import_local_ids, 
                  &export_global_ids, &export_local_ids);
  if(myrank == 0) 
  printf("passed load balance\n"); 
    migrate_triangles (size, &nt, t, v, tid, pid, 
                        num_import, import_procs, 
                        num_export, export_procs, 
                        import_global_ids, import_local_ids, 
                        export_global_ids, export_local_ids);
  if(myrank == 0) 
   printf("passed migration\n"); 
   // loba_migrateGhosts(lb, myrank, size, &nt, t, v, p, q, distance, tid, pid);
    
  if(myrank == 0) 
   printf("passed data exchange\n"); 
    integrate (step, lo, hi, nt, t, v);

  if(myrank == 0) 
   printf("passed data exchange\n"); 
    output_state(lb, myrank, nt, t, v, timesteps);
    timesteps++;
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0)//have to make sure all ranks finished
  {
    printf("Computation Finished.\n");
    //postProcessing(nprocs, size, timesteps);
    //printf("Post-Processing Finished.\n");
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

  free (distance);
  
  Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs, &export_global_ids, &export_local_ids, &export_procs);
  
  MPI_Finalize ();

  return 0;
}
