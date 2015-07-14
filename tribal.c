#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <limits.h>
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
                
                tid[i] = i;
            }
        }
    } while (ch != EOF);
    return nt;
}

void write_pointsVTK(unsigned int nt, REAL *t[3][3], REAL *v[3], unsigned int timesteps)
{
    for(unsigned int i=0;i<timesteps;i++)
    {
      char iter[15];
      sprintf(iter, "%u", i);
      char filename[50] = "output/output"; //care or buffer overflow
      strcat(filename, iter);
      strcat(filename, ".vtk");
      printf("%s\n", filename);
        printf("%i\n", i);
      
      FILE *fp = fopen(filename, "w+");
      
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
}

void migrate_status(int myrank, unsigned int nt, REAL *t[3][3], unsigned int timesteps, int num_import, int *import_procs, int num_export, int *export_procs, unsigned int *export_local_id) 
{  
  /*
  if(myrank != 0)
  {
    //send import and export to rank0   
    MPI_Send(&num_import, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    MPI_Send(&num_export, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 
    MPI_Send(&nt, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
  } else if(myrank == 0)
  {
    int **nparticle, **import, **export;
    ERRMEM (nparticle = (int **) malloc(timesteps * sizeof(int *)));
    ERRMEM (import = (int **) malloc(timesteps * sizeof(int *)));
    ERRMEM (export = (int **) malloc(timesteps * sizeof(int *)));
        
    for(unsigned int i=0;i<timesteps;i++)
    {
      ERRMEM (nparticle[i] = (int *) malloc(num_lid_entries * sizeof(REAL))); 
      ERRMEM (import[i] = (int *) malloc(num_lid_entries * sizeof(REAL)));
      ERRMEM (export[i] = (int *) malloc(num_lid_entries * sizeof(REAL)));
    }
    
    int nnodes;
    MPI_Comm_size(MPI_COMM_WORLD, &nnodes);
    for(unsigned int i=0;i<timesteps;i++)
    {
      for(int j=0;j<nnodes;j++)
      {
        if(j != 0)
        {
          int rank_import, rank_export, nn;
          //receive import and export;
          MPI_Recv(&rank_import, 1, MPI_INT, j, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          MPI_Recv(&rank_export, 1, MPI_INT, j, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
          MPI_Recv(&nn, 1, MPI_INT, j, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          nparticle[i][j] = nn;
          import[i][j] = rank_import;
          export[i][j] = rank_export;
        }
        else
        {
          nparticle[i][j] = nt;
          import[i][j] = num_import;
          export[i][j] = num_export;
        }
      }
    }
  
    FILE *fp = fopen("mpi.csv", "w+");
    
    for(int i=0;i<nnodes;i++)
    {
      fprintf(fp,"ID, nt[%u], imp[%u], exp[%u]", i, i, i);
      if(i!=nnodes-1)
      {
        fprintf(fp,", ");
      }
    }
    fprintf(fp,"\n");
    
    //replace with time/step for all timesteps
    for(unsigned int i=0;i<timesteps;i++)
    { 
      for(int j=0;j<nnodes;j++)
      {
        if(j == nnodes-1)
        { 
          fprintf(fp,"%u, %i, %i, %i\n", i, nparticle[i][j], import[i][j], export[i][j]);
        }
        else
        {
          fprintf(fp,"%u, %i, %i, %i, ", i, nparticle[i][j], import[i][j], export[i][j]);
        }
      }
    } 
    fclose(fp); 
  }*/
}

/* migrate triangles "in-place" to new ranks */
static void migrate_triangles (unsigned int size, unsigned int *nt, REAL *t[3][3], REAL *v[3],  
                              REAL *p[3], REAL *q[3], unsigned int *tid, unsigned int *pid,  
                              int num_import, int *import_procs, int num_export, int *export_procs, 
                              ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
                              ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids)
{
  int nproc, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  
  //allocate memory for tmp buffers
  int **send_idx, *pivot, *tid_buffer; 
  REAL *tbuffer[3], *vbuffer, *pbuffer, *qbuffer;
  tbuffer[0] = (REAL *) malloc(nproc*size*3*sizeof(REAL));
  tbuffer[1] = (REAL *) malloc(nproc*size*3*sizeof(REAL));
  tbuffer[2] = (REAL *) malloc(nproc*size*3*sizeof(REAL)); 
  vbuffer = (REAL *) malloc(nproc*size*3*sizeof(REAL));
  //initially there is nothing in p,q buffers
  pbuffer = (REAL *) malloc(nproc*size*3*sizeof(REAL));
  qbuffer = (REAL *) malloc(nproc*size*3*sizeof(REAL));

  send_idx = (int **) malloc(nproc*sizeof(int*));
  pivot = (int *) malloc(nproc*sizeof(int));

  for(int i=0;i<nproc;i++)
  {
    pivot[i] = 0;
    send_idx[i] = (int *) malloc(size*sizeof(int));
  } 
  //prepare export buffers
  int *export_unique_procs = (int*) malloc(nproc*sizeof(int));
  int num_export_unique=0;
  int idx=0;
  for (unsigned int i = 0; i < num_export; i++) 
  {
    int proc = export_procs[i];
    int exists = -1;
    for(unsigned int j = 0; j < idx;j++)
    {
      if(proc == export_unique_procs[j]) 
      {
        exists = 1;
        break;
      } else if (j == idx-1)
      {
        exists = 0;
      }
    }
 
    if(exists == 0 || i == 0)
    {
      export_unique_procs[idx] = proc;
      num_export_unique++;
      idx++;
      if(num_export_unique > 2)
      { //sort
        if(export_unique_procs[idx-1] > export_unique_procs[idx-2])
        {//swap
          int tmp = export_unique_procs[idx-2];
          export_unique_procs[idx-2] = export_unique_procs[idx-1];
          export_unique_procs[idx-1] = tmp;
        }
      }
    }
    //set send indices and pivots for buffers
    send_idx[proc][pivot[proc]] = export_local_ids[i];
    pivot[proc]++;
    
    //mark exported tid
    tid[export_local_ids[i]] = UINT_MAX;
  }

  //assign values to tmp export buffers
  for(int i=0;i<nproc;i++)
  {
    for(unsigned int j=0;j<pivot[i];j++)
    {
      for(int k=0;k<3;k++)
      {
        tbuffer[0][(i*nproc*pivot[i]*3)+j+(k+1)] = t[0][k][send_idx[i][j]];        
        tbuffer[1][(i*nproc*pivot[i]*3)+j+(k+1)] = t[1][k][send_idx[i][j]]; 
        tbuffer[2][(i*nproc*pivot[i]*3)+j+(k+1)] = t[2][k][send_idx[i][j]]; 
        
        vbuffer[(i*nproc*pivot[i]*3)+j+(k+1)] = v[k][send_idx[i][j]];
        pbuffer[(i*nproc*pivot[i]*3)+j+(k+1)] = p[k][send_idx[i][j]];
        qbuffer[(i*nproc*pivot[i]*3)+j+(k+1)] = q[k][send_idx[i][j]];
      }
    }
  }
  
  ///////////////////////////////////////////////
  //refine local arrays and ids (memory gaps)
  unsigned int pv=*nt-1;
  for(unsigned int i=0;i<num_export;i++)
  {//be cautious bug may be hidden here;
  
    unsigned int j=0;
    for(j=pv;j>export_local_ids[i];j--)
    {
      if(tid[j] != UINT_MAX)
      {
        tid[export_local_ids[i]] = tid[j];
        tid[j] = UINT_MAX;
        for(int k=0;k<3;k++)
        {
          t[0][k][export_local_ids[i]] = t[0][k][j];
          t[1][k][export_local_ids[i]] = t[1][k][j];
          t[2][k][export_local_ids[i]] = t[2][k][j];
          v[k][export_local_ids[i]] = v[k][j]; 
          p[k][export_local_ids[i]] = p[k][j]; 
          q[k][export_local_ids[i]] = q[k][j];
        }
        break;
      }
    }
  }
  
  //unsigned int ntexport = *nt - num_export;
  //printf("AFTER EXPORT Rank[%i], n: %i\n", myrank, ntexport);
 
  //////////////////////////////////////////////////////////////////////////
  //prepare import buffers
  int *import_unique_procs = (int*) malloc(nproc*sizeof(int));
  int num_import_unique=0;
  idx=0;
  unsigned int receive_idx = *nt - num_export; //set to last id
  for(unsigned int i=0; i < num_import; i++)
  {
    int proc = import_procs[i];
    int exists = -1;
    for(unsigned int j = 0; j < idx; j++)
    {
      if(proc == import_unique_procs[j])
      {
        exists = 1;
        break;
      } else if(j == idx-1)
      {
        exists = 0;
      }
    }

    if(exists == 0 || i == 0) 
    {
      import_unique_procs[idx] = proc;
      num_import_unique++;
      idx++;
      if(num_import_unique > 2)
      { //sort
        if(import_unique_procs[idx-1] > import_unique_procs[idx-2])
        {//swap
          int tmp = import_unique_procs[idx-2];
          import_unique_procs[idx-2] = import_unique_procs[idx-1];
          import_unique_procs[idx-1] = tmp;
        }
      }
    }
    tid[receive_idx] = import_global_ids[i];
    receive_idx++;
  }
  if(num_import > 0) *nt = *nt + num_import;//set new nt
  if(num_export > 0) *nt = *nt - num_export;
 
  //for(int i=0;i < num_import_unique;i++)
    //printf("%d:RANK[%d] - import_ID:%d, local_id:%d\n", i, myrank, import_unique_procs[i], import_local_ids[i]);
  //for(int i=0;i < num_export_unique;i++)
    //printf("%d:RANK[%d] - export_ID:%d, local_id:%d\n", i, myrank, export_unique_procs[i], export_local_ids[i]);

  //printf("num_import_unique: %d\n", num_import_unique);
  //printf("num_export_unique: %d\n", num_export_unique);
  
  for(int x=0;x<num_export_unique;x++)
  {
    if(pivot[export_unique_procs[x]] > 0)
    {//safe check
      int i = export_unique_procs[x];
      MPI_Send(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&tbuffer[0][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      MPI_Send(&tbuffer[1][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
      MPI_Send(&tbuffer[2][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
      //for(int x=0;x<pivot[i];x++) printf("sent-tbuffer: %f\n", tbuffer[0][(i*nproc*pivot[i]*3)+x+(0+1)]);

      MPI_Send(&vbuffer[(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
      MPI_Send(&pbuffer[(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
      MPI_Send(&qbuffer[(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
    }
  }
  receive_idx = *nt-receive_idx; // set to last id
  
  for(int x=0;x<num_import_unique;x++)
  {
    int i = import_unique_procs[x];
    //printf("RANK[%d]: receive from rank %d\n", myrank, i);
    MPI_Recv(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tbuffer[0][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tbuffer[1][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tbuffer[2][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //for(int x=0;x<pivot[i];x++) printf("received-tbuffer: %f\n", tbuffer[0][(i*nproc*pivot[i]*3)+x+(0+1)]);
    
    MPI_Recv(&vbuffer[(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pbuffer[(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&qbuffer[(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for(unsigned int j=0;j<pivot[i];j++)
    {
      for(int k=0;k<3;k++)
      {
        t[0][k][receive_idx] = tbuffer[0][(i*nproc*pivot[i]*3)+j+(k+1)];        
        t[1][k][receive_idx] = tbuffer[1][(i*nproc*pivot[i]*3)+j+(k+1)]; 
        t[2][k][receive_idx] = tbuffer[2][(i*nproc*pivot[i]*3)+j+(k+1)]; 
        
        v[k][receive_idx] = vbuffer[(i*nproc*pivot[i]*3)+j+(k+1)];
        p[k][receive_idx] = pbuffer[(i*nproc*pivot[i]*3)+j+(k+1)];
        q[k][receive_idx] = qbuffer[(i*nproc*pivot[i]*3)+j+(k+1)];
      }
      receive_idx++;
    }
  }

  for(int i=0; i<3;i++)
  {//free memory
    free(tbuffer[i]);
  }
    free(pivot);
    free(vbuffer);
    free(send_idx); 
}

int main (int argc, char **argv)
{
  REAL *t[3][3]; /* triangles */
  REAL *v[3]; /* velocities */
  REAL *d; /*distance */
  REAL *p[3],*q[3];
  unsigned int *tid; /* triangle identifiers */
  REAL lo[3] = {0, 0, 0}; /* lower corner */
  REAL hi[3] = {1, 1, 1}; /* upper corner */
  unsigned int nt; /* number of triangles */
  int *rank; /* migration ranks */
  unsigned int *pid; /*particle identifier */
  unsigned int size; /* ranks size */
  int myrank;

  /* init */ 
  MPI_Init (&argc, &argv);

  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  
  if (myrank == 0)
  {
    /* set nt */
    if (argc > 1) nt = atoi (argv[1]);
    else nt = 1000000;

    /* buffers */
    size = 10*nt;

    for (int i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = (REAL *) malloc (sizeof(REAL[size])));
    
      ERRMEM (p[i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (q[i] = (REAL *) malloc (sizeof(REAL[size])));
    }
    ERRMEM (d = (REAL *) malloc (nt * nt * sizeof(REAL)));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX;
    
    nt = load_pointsVTK(t, tid);
    /* generate triangles and velocities */
    //generate_triangles_and_velocities (lo, hi, nt, t, v, tid, pid);
  }
  else
  {
    /* set nt */
    nt = 0;

    /* buffers */
    if (argc > 1) size = atoi (argv[1])*4;
    else size = 1000000*4;

    for (int i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = (REAL *) malloc (sizeof(REAL[size])));
      
      ERRMEM (p[i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (q[i] = (REAL *) malloc (sizeof(REAL[size])));
    }
    ERRMEM (d = (REAL *) malloc (size * size * sizeof(REAL)));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
      
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX;
  }
  
  int  num_import, *import_procs, num_export, *export_procs;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids, export_local_ids;
  
  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  REAL step = 1E-3, time; unsigned int timesteps=0;
  
  //for (time = 0.0; time < 1.0; time += step)
  for(time = 0; time < 2; time++)
  {
    loba_balance (lb, nt, t[0], tid, 1.1, 
                  &num_import, &import_procs, 
                  &num_export, &export_procs, 
                  &import_global_ids, &import_local_ids, 
                  &export_global_ids, &export_local_ids);
    
    //printf("RANK[%i]:num_import:%d\n", myrank, num_import);
    //printf("RANK[%i]:num_export:%d\n", myrank, num_export);
   
    for(int i = 0; i < nt; i++)
    {
      printf("Before - RANK[%i]: tid = %i, t[0][0][i] = %f, p[0] = %f\n", myrank, tid[i], t[0][0][i], p[0][i]);
    }
   
    migrate_triangles (size, &nt, t, v, p, q, tid, pid, num_import, import_procs, num_export, export_procs, import_global_ids, import_local_ids, export_global_ids, export_local_ids);
    
    printf("RANK[%i]:NT:%i\n\n\n\n", myrank, nt);
    
    //contact_distance(nt, t, p, q, d);
      
    for(int i = 0; i < nt; i++)
    {
      printf("After - RANK[%i]: tid = %i, t[0][0][i] = %f, p[0] = %f\n", myrank, tid[i], t[0][0][i], p[0][i]);
    }
    
    //integrate_triangles (step, lo, hi, nt, t, v);
    
    timesteps++;
  }
  
  if(myrank == 0)
  write_pointsVTK(nt, t, v, timesteps);

  /* finalise */
  loba_destroy (lb);

  for (int i = 0; i < 2; i ++)
  {
    free (t[0][i]);
    free (t[1][i]);
    free (t[2][i]);
    free (v[i]);
    free (p[i]);
    free (q[i]);
  }

  Zoltan_LB_Free_Data (&import_global_ids, &import_local_ids, &import_procs, &export_global_ids, &export_local_ids, &export_procs);
  
  MPI_Finalize ();

  return 0;
}
