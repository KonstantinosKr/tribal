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

#define LARGENUM 100000000000 //100 billion

//using namespace ispc;

/* migrate triangles "in-place" to new ranks */
static void migrate_triangles (unsigned long long int size, unsigned int *nt, iREAL *t[3][3], iREAL *v[3],  
                              iREAL *p[3], iREAL *q[3], unsigned int *tid, unsigned int *pid,  
                              int num_import, int *import_procs, int num_export, int *export_procs, 
                              ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
                              ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids)
{
  int nproc, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  
  //allocate memory for tmp buffers
  int **send_idx, *pivot, *tid_buffer; 
  iREAL *tbuffer[3], *vbuffer, *pbuffer, *qbuffer;
  tbuffer[0] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  tbuffer[1] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  tbuffer[2] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL)); 
  vbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  //initially there is nothing in p,q buffers
  pbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  qbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));

  send_idx = (int **) malloc(nproc*sizeof(int*));
  pivot = (int *) malloc(nproc*sizeof(int));

  for(int i=0;i<nproc;i++)
  {
    pivot[i] = 0;
    send_idx[i] = (int *) malloc(size*sizeof(int));
  }

  //prepare export buffers
  int *export_unique_procs = (int*) malloc(nproc*sizeof(int));
  for(int i=0;i<nproc;i++)
  {
    export_unique_procs[i] = -1;
  }
  int num_export_unique=0;//number of unique ids to export
  int idx=0;
  for (unsigned int i = 0; i < num_export; i++)//loop through export data/ids 
  {
    int proc = export_procs[i]; //proc is the export process for data id[i]
    int exists = 0; //set to 0 to mean doesn't exist
    for(unsigned int j = 0; j < nproc;j++)
    {
      if(proc == export_unique_procs[j])//search list of unique export for duplicates 
      {
        exists = 1; //if exist skip following code to not create a duplicate
        break;
      }
    }

    //if proc is not in export_unique_procs then list it
    if(exists == 0)
    {
      export_unique_procs[num_export_unique] = proc; ///save proc id to array of unique export processes 
      num_export_unique++;//increase number of unique export processes
      idx++;//increase index number of processed unique processes
      
      //sort unique procs increamentally
      if(num_export_unique > 1)
      { //sort binary-style
        if(export_unique_procs[idx-2] > export_unique_procs[idx-1])
        {//swap
          int tmp = export_unique_procs[idx-2];
          export_unique_procs[idx-2] = export_unique_procs[idx-1];
          export_unique_procs[idx-1] = tmp;
        }
      }
    }

    //set send indices and pivot for buffers for each export process
    send_idx[proc][pivot[proc]] = export_local_ids[i];
    pivot[proc]++;
    
    //mark tid that will be exported thus deleted
    tid[export_local_ids[i]] = UINT_MAX; 
  }
  
  //assign values to tmp export buffers
  for(int i=0;i<nproc;i++)//n processes to prepare buffers for
  {
    for(unsigned int j=0;j<pivot[i];j++)//pivot gives n number of ids to loop through
    {
      for(int k=0;k<3;k++)//loop through the xyz axis
      {
        tbuffer[0][(i*size*3)+(j*3)+k] = t[0][k][send_idx[i][j]]; //point 0/A        
        tbuffer[1][(i*size*3)+(j*3)+k] = t[1][k][send_idx[i][j]]; //point 1/B
        tbuffer[2][(i*size*3)+(j*3)+k] = t[2][k][send_idx[i][j]]; //point 2/C

        ///printf("POSITION:%i\n\n", (j*3)+k);
        vbuffer[(i*size*3)+(j*3)+(k)] = v[k][send_idx[i][j]];
        pbuffer[(i*size*3)+(j*3)+(k)] = p[k][send_idx[i][j]];
        qbuffer[(i*size*3)+(j*3)+(k)] = q[k][send_idx[i][j]];
      }
    }
  }

  ///////////////////////////////////////////////
  //refine local arrays and ids (memory gaps)
  unsigned int pv=*nt;
  for(unsigned int i=0;i<num_export;i++)
  {//be cautious bug may be hidden here;
  
    for(unsigned int j=pv;j>export_local_ids[i];j--)//from last towards first but only until gap of exported
    {
      if(tid[j] != UINT_MAX)//if not marked as to be exported switch fill gaps
      {
        tid[export_local_ids[i]] = tid[j]; //send from 'last to first' the tids to 'first to last' in tid array
        tid[j] = UINT_MAX; //mark moved tid
        for(int k=0;k<3;k++)
        {
          t[0][k][export_local_ids[i]] = t[0][k][j];
          t[1][k][export_local_ids[i]] = t[1][k][j];
          t[2][k][export_local_ids[i]] = t[2][k][j];
//          printf("EXECUTED:t[2][i][local_id]:%f\n", t[2][k][export_local_ids[i]]);
          v[k][export_local_ids[i]] = v[k][j]; 
          p[k][export_local_ids[i]] = p[k][j]; 
          q[k][export_local_ids[i]] = q[k][j];
        }
        break;//break loop/search from 'last to first' and go to next id that was exported
      }
    }
  }

  //////////////////////////////////////////////////////////////////////////
  //prepare import buffers
  int *import_unique_procs = (int*) malloc(nproc*sizeof(int));
  for(int i=0;i<nproc;i++)
  {
    import_unique_procs[i] = -1;
  }
  int num_import_unique=0;
  idx=0;
  unsigned int receive_idx=0;
  if(*nt > 0 && num_export > 0) 
  {
    receive_idx = *nt - num_export; //set to last id
  } else if(*nt >= 0 && num_export <= 0){
    receive_idx = *nt;
  }

  for(unsigned int i=0; i < num_import; i++) //loop throught imports
  {
    int proc = import_procs[i]; //get process of import id i
    int exists = 0; //set exists to unknown
    for(unsigned int j = 0; j < nproc; j++) //loop through until idx reached
    {
      if(proc == import_unique_procs[j])
      {
        exists = 1;
        break;
      }
    }

    if(exists == 0) //if doesn't exist or first id i
    {
      import_unique_procs[num_import_unique] = proc; //assign proc to import unique process
      num_import_unique++; //increase unique import process
      idx++; //increase number of processes
      if(num_import_unique > 1) //if number of import processes are bigger than 2
      { //sort
        if(import_unique_procs[idx-2] > import_unique_procs[idx-1])
        {//swap
          int tmp = import_unique_procs[idx-2];
          import_unique_procs[idx-2] = import_unique_procs[idx-1];
          import_unique_procs[idx-1] = tmp;
        }
      }
    }
    tid[receive_idx] = import_global_ids[i];//set import id to tid
    receive_idx++;//increase receive index
  }

  printf("RANK[%i]: NUM_IMPORT:%i\n", myrank, num_import);
  printf("RANK[%i]: NUM_EXPORT:%i\n", myrank, num_export);
  for(int i=0;i < num_import;i++)
  {
    printf("%d:RANK[%d] - import_ID:%d, local_id:%d\n", i, myrank, import_procs[i], import_global_ids[i]);
  }
  for(int i=0;i < num_export;i++)
  {
    printf("%d:RANK[%d] - export_ID:%d, local_id:%d\n", i, myrank, export_procs[i], export_global_ids[i]);
  }
  
  if(*nt > 0 && num_export > 0) 
  {
    receive_idx = *nt - num_export; //set to last id
  } else if(*nt >= 0 && num_export <= 0){
    receive_idx = *nt;
  }
  MPI_Request *myRequest = (MPI_Request*) malloc(num_export_unique*7*sizeof(MPI_Request));//7 sends
  MPI_Request *myrvRequest = (MPI_Request*) malloc(num_import_unique*7*sizeof(MPI_Request));//7 sends 
  
  for(int x=0;x<num_export_unique;x++)
  {
    int i = export_unique_procs[x];
    MPI_Send(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
  }
  
  for(int x=0;x<num_import_unique;x++)
  {
    int i = import_unique_procs[x];  
    MPI_Recv(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  for(int x=0;x<num_import_unique;x++)
  {
    int i = import_unique_procs[x];  
    printf("RANK[%d]: receive from rank %d\n", myrank, i);
   
    //Asychronous Communication
    MPI_Irecv(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myrvRequest[(x*7)+1]);
    MPI_Irecv(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &myrvRequest[(x*7)+2]);
    MPI_Irecv(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &myrvRequest[(x*7)+3]);
    MPI_Irecv(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &myrvRequest[(x*7)+4]);
    MPI_Irecv(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &myrvRequest[(x*7)+5]);
    MPI_Irecv(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &myrvRequest[(x*7)+6]);
  }
  
  for(int x=0;x<num_export_unique;x++)
  {
    int i = export_unique_procs[x];
  
    printf("RANK[%d]: send to rank %d\n", myrank, i);
    //Asychronous Communication
    MPI_Isend(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myRequest[(x*7)+1]);
    MPI_Isend(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &myRequest[(x*7)+2]);
    MPI_Isend(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &myRequest[(x*7)+3]);  
    MPI_Isend(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &myRequest[(x*7)+4]);
    MPI_Isend(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &myRequest[(x*7)+5]);
    MPI_Isend(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &myRequest[(x*7)+6]);
  }
  
  for(int x=0;x<num_import_unique;x++)
  {
    MPI_Wait(&myrvRequest[(x*7)+1], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+2], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+3], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+4], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+5], MPI_STATUS_IGNORE); 
    MPI_Wait(&myrvRequest[(x*7)+6], MPI_STATUS_IGNORE);
    int i = import_unique_procs[x];

    for(unsigned int j=0;j<pivot[i];j++)
    {
      for(int k=0;k<3;k++)
      {
        t[0][k][receive_idx] = tbuffer[0][(i*size*3)+(j*3)+(k)];        
        t[1][k][receive_idx] = tbuffer[1][(i*size*3)+(j*3)+(k)]; 
        t[2][k][receive_idx] = tbuffer[2][(i*size*3)+(j*3)+(k)]; 
        
        v[k][receive_idx] = vbuffer[(i*size*3)+(j*3)+(k)];
        p[k][receive_idx] = pbuffer[(i*size*3)+(j*3)+(k)];
        q[k][receive_idx] = qbuffer[(i*size*3)+(j*3)+(k)];
      }
      receive_idx++;
    }
  }
  
  for(int x=0;x<num_export_unique;x++)
  {
    MPI_Wait(&myRequest[(x*7)+1], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*7)+2], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*7)+3], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*7)+4], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*7)+5], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*7)+6], MPI_STATUS_IGNORE); 
  }

  *nt = *nt + (num_import-num_export);
  
  for(int i=0; i<3;i++)
  {//free memory
    free(tbuffer[i]);
  }
    free(pivot);
    free(vbuffer);
    free(send_idx); 
    free(myRequest);
    free(myrvRequest);
}

int main (int argc, char **argv)
{
  iREAL *t[3][3]; /* triangles */
  iREAL *v[3]; /* velocities */
  iREAL *distance; /*distance */
  iREAL *p[3],*q[3];
  unsigned int *tid; /* triangle identifiers */
  iREAL lo[3] = {-255, -255, -255}; /* lower corner */
  iREAL hi[3] = {255, 255, 255}; /* upper corner */
  iREAL mylo[3], myhi[3];
  unsigned int nt; /* number of triangles */
  int *rank; /* migration ranks */
  unsigned int *pid; /*particle identifier */
  unsigned int size; /* ranks size */
  int nprocs;
  int nNeighbors=0;
  int myrank;

  /* init */ 
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  
  int *neighborhood = (int *) malloc(sizeof(int[nprocs]));
  
  if (myrank == 0)
  {
    /* set nt */
    if (argc > 1) nt = atoi (argv[1]);
    else nt = 1000000;

    /* buffers */
    size = nt;

    for (int i = 0; i < 3; i ++)
    { 
      ERRMEM (t[0][i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (t[1][i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (t[2][i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (v[i] = (iREAL *) malloc (sizeof(iREAL[size])));
    
      ERRMEM (p[i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (q[i] = (iREAL *) malloc (sizeof(iREAL[size])));
    }
    ERRMEM (distance = (iREAL *) malloc (sizeof(iREAL[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX;
    
    //generate_triangles_and_velocities (lo, hi, nt, t, v, tid, pid);
    
    iREAL mint, maxt;
    nt = load_pointsVTK(t, tid, &mint, &maxt); 
    //normalize(nt, t, mint, maxt); 
   
    gen_velocities(lo, hi, nt, v);
  }
  else
  {
    /* set nt */
    nt = 0;

    /* buffers */
    if (argc > 1) size = atoi (argv[1])*4;
    else size = 1000000;

    for (int i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (t[1][i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (t[2][i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (v[i] = (iREAL *) malloc (sizeof(iREAL[size])));
      
      ERRMEM (p[i] = (iREAL *) malloc (sizeof(iREAL[size])));
      ERRMEM (q[i] = (iREAL *) malloc (sizeof(iREAL[size])));
    }
    ERRMEM (distance = (iREAL*) malloc (sizeof(iREAL[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
      
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX;
  }
  
  int num_import, *import_procs, num_export, *export_procs;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids, export_local_ids;
  
  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  iREAL step = 1E-3, time; unsigned int timesteps=0;
  
  //for (time = 0.0; time < 1.0; time += step)
  for(time = 0; time < 2; time++)
  {
    if(myrank == 0){printf("\n\n\n\n\nTIMESTEP: %i\n\n\n\n", timesteps);}

    loba_balance (lb, nt, t[0], tid, 1.1,
                  &num_import, &import_procs, 
                  &num_export, &export_procs, 
                  &import_global_ids, &import_local_ids, 
                  &export_global_ids, &export_local_ids);
    printf("passed load balance\n");
    
    if(time>0)
    for(int i = 0; i < nt+1; i++)
    {
      printf("Before - RANK[%i]: tid = %i\nt[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f,\nt[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f,\nt[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n\n", myrank, tid[i], t[0][0][i], t[0][1][i], t[0][2][i], t[1][0][i], t[1][1][i], t[1][2][i], t[2][0][i], t[2][1][i], t[2][2][i]);
    }

    migrate_triangles (size, &nt, t, v, p, q, tid, pid, 
                        num_import, import_procs, 
                        num_export, export_procs, 
                        import_global_ids, import_local_ids, 
                        export_global_ids, export_local_ids);
    if(time>0)
    for(int i = 0; i < nt; i++)
    {
      printf("After - RANK[%i]: tid = %i\nt[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f,\nt[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f,\nt[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n\n", myrank, tid[i], t[0][0][i], t[0][1][i], t[0][2][i],t[1][0][i], t[1][1][i], t[1][2][i], t[2][0][i], t[2][1][i], t[2][2][i]);
    }

    //printf("passed migration\n");
    //loba_getAdjacent(lb, myrank, neighborhood, &nNeighbors);
    //loba_migrateGhosts(lb, myrank, neighborhood, nNeighbors, size, &nt, t, v, p, q, distance, tid, pid);
    
    //printf("passed ghosts migration\n"); 
    //integrate (step, lo, hi, nt, t, v);
    loba_getbox (lb, myrank, mylo, myhi);//get local subdomain boundary box
    
    write_pointsVTK(myrank, nt, t, v, mylo, myhi, timesteps);
      
    timesteps++;
  }
  
  printf("finished\n");
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
