#ifndef __rcb__
#define __rcb__

#include <zoltan.h>

enum algo
{
  ZOLTAN_RCB,
  ZOLTAN_RIB
};

struct loba
{
  struct Zoltan_Struct *zoltan; /* Zoltan load balancer */

  enum algo al;
};

/* create load balancer */
struct loba* loba_create (enum algo al);

/* balance points up to tolerance; output migration ranks */
//void loba_balance (struct loba *lb, unsigned int n, REAL *p[3], unsigned int *id, REAL tol, int *rank);

void loba_balance (struct loba *lb, unsigned int n, REAL *p[3], unsigned int *id, REAL tol, int *rank, 
                    int *num_gid_entries, int *num_lid_entries, int *num_import, int *import_procs, int *num_export, int *export_procs,
                    unsigned int *export_local_id, ZOLTAN_ID_PTR *import_global_ids, ZOLTAN_ID_PTR *import_local_ids, 
                    ZOLTAN_ID_PTR *export_global_ids, ZOLTAN_ID_PTR *export_local_ids); 

/* find ranks overlapped by the [lo,hi] box */
void loba_query (struct loba *lb, int node, REAL lo[3], REAL hi[3], int *nranks, int *ranks);

/* free load balancer */
void loba_destroy (struct loba *lb);

#endif
