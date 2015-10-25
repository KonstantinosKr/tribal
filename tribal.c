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
#include "tmr.h"

int main (int argc, char **argv)
{
  iREAL *t[3][3]; /* triangles */
  iREAL *v[3]; /* velocities */
  iREAL *distance; /*distance */
  iREAL *p[3],*q[3];//p and q points
  unsigned int nt = 0; /* number of triangles */
  unsigned int *pid; /*particle identifier */
  unsigned int *tid; /* triangle identifiers */
  iREAL lo[3] = {-255, -255, -255}; /* lower corner */
  iREAL hi[3] = {255, 255, 255}; /* upper corner */
  
  unsigned int nParticles = 12;
  unsigned int size = 20000000; /* memory buffer size */
  int nprocs, myrank;

  /* init */ 
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   
  if (myrank == 0)
  {
    size = 20000000;
    for (int i = 0; i < 3; i ++)
    { 
      t[0][i] = (iREAL *) malloc (size*sizeof(iREAL));
      t[1][i] = (iREAL *) malloc (size*sizeof(iREAL));
      t[2][i] = (iREAL *) malloc (size*sizeof(iREAL));
      v[i] = (iREAL *) malloc (size*sizeof(iREAL));
    
      p[i] = (iREAL *) malloc (size*sizeof(iREAL));
      q[i] = (iREAL *) malloc (size*sizeof(iREAL));
    }
    distance = (iREAL *) malloc (size*sizeof(iREAL));
    tid = (unsigned int *) malloc (size*sizeof(unsigned int));
    pid = (unsigned int *) malloc (size*sizeof(unsigned int));
    
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX; 
    
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
    ptype[9] = 0;
    ptype[10] = 0;
    ptype[11] = 0;
    ptype[12] = 0;
    ptype[13] = 0;
    ptype[14] = 0;
    ptype[15] = 0;
    
    load_enviroment(ptype, &nt, nParticles, t, tid, pid, &mint, &maxt);
    
    iREAL velo[3] = {50, 50, 50};
    gen_velocities(lo, velo, nt, v);
  }
  else
  {
    size = 20000000;
    for (int i = 0; i < 3; i ++)
    {
      t[0][i] = (iREAL *) malloc (size*sizeof(iREAL));
      t[1][i] = (iREAL *) malloc (size*sizeof(iREAL));
      t[2][i] = (iREAL *) malloc (size*sizeof(iREAL));
      v[i] = (iREAL *) malloc (size*sizeof(iREAL));
      
      p[i] = (iREAL *) malloc (size*sizeof(iREAL));
      q[i] = (iREAL *) malloc (size*sizeof(iREAL));
    }
    distance = (iREAL*) malloc (size*sizeof(iREAL));
    tid = (unsigned int *) malloc (size*sizeof(unsigned int));
    pid = (unsigned int *) malloc (size*sizeof(unsigned int));
      
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX;
  }
  
  int num_import, num_export, *import_procs, *export_procs;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids, export_local_ids;
  
  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  iREAL step = 1E-3, time; unsigned int timesteps=0;
  
  TIMING tbalance[100];
  TIMING tmigration[100];
  TIMING tdataExchange[100];
  TIMING tintegration[100];
  iREAL tTimer1[100];
  iREAL tTimer2[100];
  iREAL tTimer3[100];
  iREAL tTimer4[100];

  iREAL timer1, timer2, timer3;
  timer1 = 0.0;
  timer2 = 0.0;
  timer3 = 0.0;
  
  //for (time = 0.0; time < 1.0; time += step)
  for(time = 0; time < 0.1; time+=step)
  {
    if(myrank == 0){printf("TIMESTEP: %i\n", timesteps);}
    
    timerstart(&tbalance[timesteps]);
    loba_balance (lb, nt, t[0], tid, 1.1,
                  &num_import, &import_procs, 
                  &num_export, &export_procs, 
                  &import_global_ids, &import_local_ids, 
                  &export_global_ids, &export_local_ids);
    timerend (&tbalance[timesteps]);
   
    printf("RANK[%i]: load balance:%f\n", myrank, tbalance[timesteps].total);
    
    timerstart(&tmigration[timesteps]);
    migrate_triangles (size, &nt, t, v, tid, pid, 
                        num_import, import_procs, 
                        num_export, export_procs, 
                        import_global_ids, import_local_ids, 
                        export_global_ids, export_local_ids);
    timerend (&tmigration[timesteps]);
    
    printf("RANK[%i]: migration:%f\n", myrank, tmigration[timesteps].total);
    
    timer1 = 0.0;
    timer2 = 0.0;
    timer3 = 0.0;
    
    timerstart (&tdataExchange[timesteps]);
    loba_migrateGhosts(lb, myrank, size, &nt, t, v, p, q, distance, tid, pid, &timer1, &timer2, &timer3);
    timerend (&tdataExchange[timesteps]);
   
    tTimer1[timesteps] = timer1;
    tTimer2[timesteps] = timer2;
    tTimer3[timesteps] = timer3;
 
    printf("RANK[%i]: data exchange:%f\n", myrank, tdataExchange[timesteps].total);
    
    timerstart (&tintegration[timesteps]);
    integrate (step, lo, hi, nt, t, v);
    timerend (&tintegration[timesteps]);
    
    printf("RANK[%i]: integration:%f\n", myrank, tintegration[timesteps].total);

    //output_state(lb, myrank, nt, t, v, timesteps);
    
    timesteps++;
  }

  iREAL subtotal = 0;
  iREAL bal = 0;
  iREAL mig = 0;
  iREAL de = 0;
  iREAL in = 0;
  iREAL dt1 = 0;
  iREAL dt2 = 0;
  iREAL dt3 = 0;
  for(int i = 0; i<timesteps;i++)
  {
    subtotal = subtotal + tbalance[i].total + tmigration[i].total + tdataExchange[i].total + tintegration[i].total;
    bal = bal + tbalance[i].total;
    mig = mig + tmigration[i].total;
    de = de + tdataExchange[i].total;
    in = in + tintegration[i].total;
    dt1 = dt1 + tTimer1[i];
    dt2 = dt2 + tTimer2[i];
    dt3 = dt3 + tTimer3[i];
  }  

  printf("RANK[%i]: TOTAL:%f Z-BALANCE:%f, MIGRATION:%f, DATAXCHANGE:%f, DT1:%f, DT2:%f, DT3:%f, INTEGRATION:%f\n", myrank, subtotal, bal, mig, de, dt1, dt2, dt3, in);
  printf("RANK[%i]: FIRST MIGRATION:%f\n", myrank, tmigration[0].total);
  MPI_Barrier(MPI_COMM_WORLD);
  if(myrank == 0)//have to make sure all ranks finished
  {
    printf("\nComputation Finished.\n");
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
