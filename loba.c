#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include "error.h"
#include "loba.h"

struct zoltan_args
{
  unsigned int n;
  REAL *p[3];
  unsigned int *id;
};

/* number of objects for balacing */
static int obj_count (struct zoltan_args *args, int *ierr)
{
  *ierr = ZOLTAN_OK;
  
  return args->n > 0 ? args->n : 1;
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

  if (i == 0) /* XXX: Zoltan workaround */
  {
    global_ids [0] = UINT_MAX;
    obj_wgts [0] = 1.0;
  }

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

  if (num_obj == 1 && global_ids [0] == UINT_MAX) /* XXX: Zoltan workaround */
  {
    geom_vec[0] = 0;
    geom_vec[1] = 0;
    geom_vec[2] = 0;
  }
  else for (i = 0; i < num_obj; i ++)
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

  ERRMEM (lb = malloc (sizeof (struct loba)));

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
    //Zoltan_Set_Param (lb->zoltan, "RETURN_LISTS", "IMPORT");

    /* RCB parameters */
    Zoltan_Set_Param (lb->zoltan, "RCB_OVERALLOC", "1.3");
    Zoltan_Set_Param (lb->zoltan, "RCB_REUSE", "1");
    Zoltan_Set_Param (lb->zoltan, "RCB_OUTPUT_LEVEL", "0");
    Zoltan_Set_Param (lb->zoltan, "CHECK_GEOM", "1");
    Zoltan_Set_Param (lb->zoltan, "KEEP_CUTS", "1");
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
void loba_balance (struct loba *lb, unsigned int n, REAL *p[3], unsigned int *id, REAL tol, int *rank, 
                    int *num_gid_entries, int *num_lid_entries, int *num_import, int *import_procs, int *num_export, int *export_procs, 
                    unsigned int *export_local_id) 
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
    
    int changes; //num_gid_entries, num_lid_entries, num_import, *import_procs, num_export, *export_procs;
    ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids, export_local_ids;
    
    /* update partitioning */
    ASSERT (Zoltan_LB_Balance (lb->zoltan, &changes, num_gid_entries, num_lid_entries,
	    num_import, &import_global_ids, &import_local_ids, &import_procs,
	    num_export, &export_global_ids, &export_local_ids, &export_procs) == ZOLTAN_OK, "Zoltan load balancing failed");

    unsigned int myrank;
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
    
    for (int i = 0; i < n; i ++) rank[i] = myrank;
    for (int i = 0; i < *num_export; i++) 
    {
      rank[export_local_ids[i]] = export_procs[i];
      //load to transfer
      export_local_id[i] = export_local_ids[i];
    }
   
    printf("RANK[%i]:num_import:%d\n", myrank, *num_import);
    printf("RANK[%i]:num_export:%d\n", myrank, *num_export);
    
    //Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs, &export_global_ids, &export_local_ids, &export_procs);
  }
  break;
  case ZOLTAN_RIB:
  {
  }
  }
}

/* find ranks overlapped by the [lo,hi] box */
void loba_query (struct loba *lb, int node, REAL lo[3], REAL hi[3], int *nranks, int *ranks)
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

/* free load balancer */
void loba_destroy (struct loba *lb)
{
  if (lb->zoltan) Zoltan_Destroy (&lb->zoltan);

  free (lb);
}
