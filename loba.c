#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <float.h>
#include "error.h"
#include "loba.h"

struct zoltan_args
{
  unsigned int n;
  iREAL *p[3];
  unsigned int *id;
};

/* number of objects for balacing */
static int obj_count (struct zoltan_args *args, int *ierr)
{
  *ierr = ZOLTAN_OK;

#if 0
  return args->n > 0 ? args->n : 1;
#else
  return args->n;
#endif
}

/* list of object identifiers for load balancing */
static void obj_list (struct zoltan_args *args, int num_gid_entries, int num_lid_entries,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float *obj_wgts, int *ierr)
{
  int i;
  
  for (i = 0; i < args->n; i ++)
  {
    global_ids [i * num_gid_entries] = args->id[i];
    local_ids [i * num_lid_entries] = i;
    obj_wgts [i * wgt_dim] = 1.0;
  }

#if 0
  if (i == 0) /* XXX: Zoltan workaround */
  {
    global_ids [0] = UINT_MAX;
    obj_wgts [0] = 1.0;
  }
#endif

  *ierr = ZOLTAN_OK;
}

/* number of spatial dimensions */
static int dimensions (struct zoltan_args *args, int *ierr)
{
  *ierr = ZOLTAN_OK;
  return 3;
}

/* list of object points exploited during load balancing */
static void obj_points (struct zoltan_args *args, int num_gid_entries, int num_lid_entries, int num_obj,
  ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int num_dim, double *geom_vec, int *ierr)
{
  double *v;
  int i, j;

#if 0
  if (num_obj == 1 && global_ids [0] == UINT_MAX) /* XXX: Zoltan workaround */
  {
    geom_vec[0] = 0;
    geom_vec[1] = 0;
    geom_vec[2] = 0;
  }
  else
#endif
  for (i = 0; i < num_obj; i ++)
  {
    j = local_ids [i * num_lid_entries];
    v = &geom_vec [i * num_dim];

    v[0] = args->p[0][j];
    v[1] = args->p[1][j];
    v[2] = args->p[2][j];
  }

  *ierr = ZOLTAN_OK;
}

/* create load balancer */
struct loba* loba_create (enum algo al)
{
  struct loba *lb;

  ERRMEM (lb = (struct loba*) malloc (sizeof (struct loba)));

  switch (al)
  {
  case ZOLTAN_RCB:
  {
    /* create Zoltan object */
    ASSERT (lb->zoltan = Zoltan_Create (MPI_COMM_WORLD), "Zoltan initialisation failed");

    /* set general parameters */
    Zoltan_Set_Param (lb->zoltan, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param (lb->zoltan, "DEBUG_MEMORY", "0");
    Zoltan_Set_Param (lb->zoltan, "NUM_GID_ENTRIES", "1");
    Zoltan_Set_Param (lb->zoltan, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param (lb->zoltan, "OBJ_WEIGHT_DIM", "1");

    /* load balancing parameters */
    Zoltan_Set_Param (lb->zoltan, "LB_METHOD", "RCB");
    Zoltan_Set_Param (lb->zoltan, "IMBALANCE_TOL", "1.3");
    Zoltan_Set_Param (lb->zoltan, "AUTO_MIGRATE", "FALSE");
    Zoltan_Set_Param (lb->zoltan, "RETURN_LISTS", "IMPORT AND EXPORT");
    //Zoltan_Set_Param (lb->zoltan, "RETURN_LISTS", "EXPORT");
    
    /* RCB parameters */
    Zoltan_Set_Param (lb->zoltan, "RCB_OVERALLOC", "1.3");
    Zoltan_Set_Param (lb->zoltan, "RCB_REUSE", "TRUE");
    Zoltan_Set_Param (lb->zoltan, "RCB_OUTPUT_LEVEL", "0");
    Zoltan_Set_Param (lb->zoltan, "CHECK_GEOM", "1");
    Zoltan_Set_Param (lb->zoltan, "KEEP_CUTS", "TRUE");
    Zoltan_Set_Param (lb->zoltan, "REDUCE_DIMENSIONS", "0");
  }
  break;
  case ZOLTAN_RIB:
  {
    /* TODO */
  }
break;
  }

  lb->al = al;

  return lb;
}

/* balance points up to tolerance; output migration ranks */
void loba_balance (struct loba *lb, unsigned int n, iREAL *p[3], unsigned int *id, iREAL tol,
                    int *num_import, int **import_procs, int *num_export, int **export_procs, 
                    ZOLTAN_ID_PTR *import_global_ids, ZOLTAN_ID_PTR *import_local_ids, 
                    ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids) 
{
  switch (lb->al)
  {
  case ZOLTAN_RCB:
  {
    struct zoltan_args args = {n, {p[0], p[1], p[2]}, id};
    
    /* callbacks */
    Zoltan_Set_Fn (lb->zoltan, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) obj_count, &args);
    Zoltan_Set_Fn (lb->zoltan, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) obj_list, &args);
    Zoltan_Set_Fn (lb->zoltan, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) dimensions, &args);
    Zoltan_Set_Fn (lb->zoltan, ZOLTAN_GEOM_MULTI_FN_TYPE, (void (*)()) obj_points, &args);

    /* update imbalance */
    char str[128];
    snprintf (str, 128, "%g", tol);
    Zoltan_Set_Param (lb->zoltan, "IMBALANCE_TOL", str);
    
    int changes, num_gid_entries, num_lid_entries; /* TODO: do we need this outside? */
    
    /* update partitioning */
    ASSERT (Zoltan_LB_Balance (lb->zoltan, &changes, &num_gid_entries, &num_lid_entries,
	    num_import, import_global_ids, import_local_ids, import_procs,
	    num_export, export_global_ids, export_local_ids, export_procs) == ZOLTAN_OK, "Zoltan load balancing failed");
  }
  break;
  case ZOLTAN_RIB:
  {
  }
  }
}

/* find ranks overlapped by the [lo,hi] box */
void loba_query (struct loba *lb, int node, iREAL lo[3], iREAL hi[3], int *ranks, int *nranks)
{
  switch (lb->al)
  {
    case ZOLTAN_RCB:
    {
      Zoltan_LB_Box_Assign (lb->zoltan, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], ranks, nranks);
      break;
    }
    case ZOLTAN_RIB:
    {
    }
  }
}

void loba_getAdjacent(struct loba *lb, int myrank, int *ranks, int *nprocs)
{
  iREAL mylo[3], myhi[3], lo[3], hi[3];

  loba_getbox(lb, myrank, mylo, myhi); 
  
  iREAL mypoint[8][3];
  iREAL point[8][3];
  
  int isNeighbor;
  int counter = 0;

  mypoint[0][0] = mylo[0];
  mypoint[0][1] = mylo[1];
  mypoint[0][2] = mylo[2];

  mypoint[1][0] = mylo[0];
  mypoint[1][1] = myhi[1];
  mypoint[1][2] = mylo[2];

  mypoint[2][0] = mylo[0];
  mypoint[2][1] = myhi[1];
  mypoint[2][2] = myhi[2];

  mypoint[3][0] = mylo[0];
  mypoint[3][1] = mylo[1];
  mypoint[3][2] = myhi[2]; 

  mypoint[4][0] = myhi[0];
  mypoint[4][1] = myhi[1];
  mypoint[4][2] = myhi[2];

  mypoint[5][0] = myhi[0];
  mypoint[5][1] = mylo[1];
  mypoint[5][2] = myhi[2];

  mypoint[6][0] = myhi[0];
  mypoint[6][1] = mylo[1];
  mypoint[6][2] = mylo[2];

  mypoint[7][0] = myhi[0];
  mypoint[7][1] = myhi[1];
  mypoint[7][2] = mylo[2]; 
    
  int nranks;
  MPI_Comm_size(MPI_COMM_WORLD, &nranks);
  for(int i=0; i<nranks; i++)
  {
    if(i == myrank) continue;
    loba_getbox(lb, i, lo, hi); 
    
    point[0][0] = lo[0];
    point[0][1] = lo[1];
    point[0][2] = lo[2];

    point[1][0] = lo[0];
    point[1][1] = hi[1];
    point[1][2] = lo[2];

    point[2][0] = lo[0];
    point[2][1] = hi[1];
    point[2][2] = hi[2];

    point[3][0] = lo[0];
    point[3][1] = lo[1];
    point[3][2] = hi[2]; 

    point[4][0] = hi[0];
    point[4][1] = hi[1];
    point[4][2] = hi[2];

    point[5][0] = hi[0];
    point[5][1] = lo[1];
    point[5][2] = hi[2];

    point[6][0] = hi[0];
    point[6][1] = lo[1];
    point[6][2] = lo[2];

    point[7][0] = hi[0];
    point[7][1] = hi[1];
    point[7][2] = lo[2]; 
    
    for(int j=0; j<8; j++)
    {
      for(int jj=0; j<8; j++)
      {
        for(int z=0; z<3; z++)
        {
          for(int zz=0; zz<3; zz++)
          {
            if(point[j][z] == mypoint[jj][zz])
            {
              isNeighbor = 1;
            }
          }
        }
      }
    }
    if(isNeighbor == 1)
    {
      ranks[counter] = i;
      counter++;
    }
  }
  *nprocs = counter;
}
 
void loba_getbox (struct loba *lb, int part, iREAL lo[3], iREAL hi[3])
{
  switch (lb->al)
  {
    case ZOLTAN_RCB:
    {  
      int ndim;
      Zoltan_RCB_Box(lb->zoltan, part, &ndim, &lo[0], &lo[1], &lo[2], &hi[0], &hi[1], &hi[2]);  
  
      for(int j = 0; j < 3; j++)
      {
        if(lo[j] < -FLT_MAX)
        {
          lo[j] = -255.;
        } 
        
        if(hi[j] > FLT_MAX)
        {
          hi[j] = 255.;
        } 
      }
      
      break;
    }
    case ZOLTAN_RIB:
    {
    }
  }
}

void loba_migrateGhosts(struct loba *lb, int  myrank, int *neighborhood, int nNeighbors, unsigned int size, unsigned int *nt, iREAL *t[3][3], iREAL *v[3], iREAL *p[3], iREAL *q[3], unsigned int *tid, unsigned int *pid)
{
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  //allocate memory for tmp buffers
  int **send_idx, *pivot, *rcvpivot, **tid_buffer; 
  iREAL *tbuffer[3], *vbuffer, *pbuffer, *qbuffer;
  tbuffer[0] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  tbuffer[1] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  tbuffer[2] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL)); 
  vbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  //initially there is nothing in p,q buffers
  pbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  qbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));

  send_idx = (int **) malloc(nproc*sizeof(int*));
  tid_buffer = (int **) malloc(nproc*sizeof(int*));
  pivot = (int *) malloc(nproc*sizeof(int));
  rcvpivot = (int *) malloc(nproc*sizeof(int));

  unsigned int n = *nt;
  //prepare export buffers
  for (unsigned int i = 0; i < nNeighbors; i++) 
  {
    int proc = neighborhood[i];
     
    pivot[proc] = 0;
    send_idx[proc] = (int *) malloc(nproc*size*sizeof(int));
    tid_buffer[proc] = (int *) malloc(nproc*size*sizeof(int));
    
    for(unsigned int j = 0; j < n; j++)
    {
      //set send indices and pivots for buffers
      send_idx[proc][j] = tid[j];
      pivot[proc]++;
    }
  }

  printf("passed prepare for export buffers\n");

  //assign values to tmp export buffers
  for(int i=0;i<nNeighbors;i++)//n processes to prepare buffers for
  {
    int proc = neighborhood[i];
    for(unsigned int j=0;j<pivot[proc];j++)//pivot gives n number of ids to loop through
    {
      for(int k=0;k<3;k++)//loop through the xyz axis
      {
        tbuffer[0][(proc*size*3)+(j*3)+k] = t[0][k][send_idx[proc][j]]; //point 0        
        tbuffer[1][(proc*size*3)+(j*3)+k] = t[1][k][send_idx[proc][j]]; //point 1
        tbuffer[2][(proc*size*3)+(j*3)+k] = t[2][k][send_idx[proc][j]]; //point 2

        vbuffer[(proc*size*3)+(j*3)+(k)] = v[k][send_idx[proc][j]];
        pbuffer[(proc*size*3)+(j*3)+(k)] = p[k][send_idx[proc][j]];
        qbuffer[(proc*size*3)+(j*3)+(k)] = q[k][send_idx[proc][j]];
      }
    }
  }
  MPI_Request *myRequest = (MPI_Request*) malloc(nNeighbors*7*sizeof(MPI_Request));//7 sends
  MPI_Request *myrvRequest = (MPI_Request*) malloc(nNeighbors*7*sizeof(MPI_Request));//7 sends 
  
  //blocking communication
  for(int i=0; i<nNeighbors; i++)
  {
    int j = neighborhood[i];
    
    MPI_Send(&pivot[j], 1, MPI_INT, j, 1, MPI_COMM_WORLD); 
    MPI_Recv(&rcvpivot[j], 1, MPI_INT, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
  }
  
  for(int i=0; i<nNeighbors; i++)
  {
    int j = neighborhood[i];
    MPI_Irecv(&tid_buffer[j][0], rcvpivot[j], MPI_INT, j, 2, MPI_COMM_WORLD, &myrvRequest[(j*7)]);  
    MPI_Isend(&send_idx[j][0], pivot[j], MPI_INT, j, 2, MPI_COMM_WORLD, &myRequest[(j*7)]);
    
    MPI_Wait(&myrvRequest[(j*7)], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(j*7)], MPI_STATUS_IGNORE);
  }
  
  //prepare import buffers
  unsigned int receive_idx = n; //set to last id
  for(unsigned int i=0; i < nNeighbors; i++)
  {
    int x = neighborhood[i];
    for(unsigned int j=0; j < rcvpivot[x]; j++)
    {
      tid[receive_idx] = tid_buffer[x][j]; //tids to import add them to tid
      receive_idx++;
    }
  }
  receive_idx = *nt;

 
  for(int x=0;x<nNeighbors;x++)
  { 
    if(rcvpivot[neighborhood[x]] > 0)
    {//safe check     
      int i = neighborhood[x];
      MPI_Irecv(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myrvRequest[(x*7)+1]);
      MPI_Irecv(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &myrvRequest[(x*7)+2]);
      MPI_Irecv(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &myrvRequest[(x*7)+3]);
      
      MPI_Irecv(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &myrvRequest[(x*7)+4]);
      MPI_Irecv(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &myrvRequest[(x*7)+5]);
      MPI_Irecv(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &myrvRequest[(x*7)+6]);
    }
  }

  for(int x=0;x<nNeighbors;x++)
  {
    if(pivot[neighborhood[x]] > 0)
    {//safe check
      int i = neighborhood[x];
      MPI_Isend(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myRequest[(x*7)+1]);
      MPI_Isend(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &myRequest[(x*7)+2]);
      MPI_Isend(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &myRequest[(x*7)+3]);
      
      MPI_Isend(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &myRequest[(x*7)+4]);
      MPI_Isend(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &myRequest[(x*7)+5]);
      MPI_Isend(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &myRequest[(x*7)+6]);
    }
  }
  
  for(int x=0;x<nNeighbors;x++)
  {
    MPI_Wait(&myrvRequest[(x*7)+1], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+2], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+3], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+4], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*7)+5], MPI_STATUS_IGNORE); 
    MPI_Wait(&myrvRequest[(x*7)+6], MPI_STATUS_IGNORE);
    int i = neighborhood[x];
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

  for(int x=0;x<nNeighbors;x++)
  {
    if(pivot[neighborhood[x]] > 0)
    {//safe check
      MPI_Wait(&myRequest[(x*7)+1], MPI_STATUS_IGNORE);
      MPI_Wait(&myRequest[(x*7)+2], MPI_STATUS_IGNORE);
      MPI_Wait(&myRequest[(x*7)+3], MPI_STATUS_IGNORE);
      MPI_Wait(&myRequest[(x*7)+4], MPI_STATUS_IGNORE);
      MPI_Wait(&myRequest[(x*7)+5], MPI_STATUS_IGNORE);
      MPI_Wait(&myRequest[(x*7)+6], MPI_STATUS_IGNORE); 
    }
  }
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

/* free load balancer */
void loba_destroy (struct loba *lb)
{
  if (lb->zoltan) Zoltan_Destroy (&lb->zoltan);
  free (lb);
}
