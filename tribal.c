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

void migrate_status(int myrank, unsigned int nt, REAL *t[3][3], unsigned int timesteps, int num_gid_entries, int num_lid_entries, int num_import, int *import_procs, int num_export, int *export_procs, unsigned int *export_local_id) 
{   
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
  }
}

/* migrate triangles "in-place" to new ranks */
static void migrate_triangles (int myrank, unsigned int size, unsigned int nt, REAL *t[3][3], REAL *v[3], int *rank, int num_gid_entries, int num_lid_entries, 
                                int num_import, int *import_procs, int num_export, int *export_procs, unsigned int *export_local_id) 
{
  /*
  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if(nt>0)
  for (unsigned int i = 0; i < num_export; i++) 
  {
    printf("export_id:%i\n", export_local_id[i]);   
  }
  
  int **send_idx, *pivot; 
  REAL **tbuffer[3][3], **vbuffer[3];
  for(int i=0;i<3;i++)
  {
    tbuffer[0][i] = (REAL **) malloc(nproc*sizeof(REAL*));
    tbuffer[1][i] = (REAL **) malloc(nproc*sizeof(REAL*));
    tbuffer[2][i] = (REAL **) malloc(nproc*sizeof(REAL*));
    vbuffer[i] = (REAL **) malloc(nproc*sizeof(REAL*));
  }

  send_idx = (int **) malloc(nproc*sizeof(int*));
  pivot = (int *) malloc(nproc*sizeof(int));
  for(int i=0;i<nproc;i++)
  {
    for(int j=0;j<3;j++)
    {
      tbuffer[j][0][i] = (REAL *) malloc(size*sizeof(REAL));
      tbuffer[j][1][i] = (REAL *) malloc(size*sizeof(REAL));
      tbuffer[j][2][i] = (REAL *) malloc(size*sizeof(REAL));
      vbuffer[j][i] = (REAL *) malloc(size*sizeof(REAL));
    }
    pivot[i] = 0;
    send_idx[i] = (int *) malloc(size*sizeof(int));
  }
  
  if(nt>0)
  for (unsigned int i = 0; i < num_export; i++) 
  {
    send_idx[rank[export_local_id[i]]][pivot[rank[export_local_id[i]]]] = export_local_id[i];
    pivot[rank[export_local_id[i]]]++;
  }
  
  int *pivot_buffer = (int *) malloc(nproc*sizeof(int));
  
  if(nt>0)
  for(int i=0;i<nproc;i++)
  {
    pivot_buffer[i] = pivot[i];
    for(unsigned int j=0;j<pivot[i];j++)
    {
      for(int k=0;k<3;k++)
      {
        tbuffer[0][k][i][j] = t[0][k][send_idx[i][j]];  
        
        tbuffer[1][k][i][j] = t[1][k][send_idx[i][j]]; 
        
        tbuffer[2][k][i][j] = t[2][k][send_idx[i][j]]; 
        
        vbuffer[k][i][j] = v[k][send_idx[i][j]];
      }
    }
  }

  if(myrank==0)
  {
    printf("t0:%f\n", tbuffer[0][0][1][0]);
  }

  //send buffers loop through procs is not correct, should loop through procs that need to send something/same for receive.(send has the pivot check*)
  for(int i=0;i<nproc;i++)
  {
    if(myrank != i && pivot[i] > 0)
    {  
      MPI_Send(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&tbuffer[i], pivot[i], MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      //MPI_Send(&vbuffer[i], pivot[i], MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
    }     
  }

  for(int i=0;i<nproc;i++)
  {
    if(myrank != i && myrank != 0)
    {
      MPI_Recv(&pivot_buffer[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&tbuffer[0][0][i][0], pivot_buffer[i], MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int x=0;x<pivot_buffer[i];x++)
      {
        printf("received-tbuffer: %f\n", tbuffer[0][0][0][x]);
      }
        //MPI_Recv(&vbuffer[i], pivot_buffer[i], MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
   // MPI_Sendrecv(&pivot_buffer[0], 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &pivot_buffer[0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); 
   // MPI_Sendrecv(&tbuffer[0], pivot[0], MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &tbuffer[0], pivot_buffer[0], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
   // MPI_Sendrecv(&vbuffer[0], pivot[0], MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &vbuffer[0], pivot_buffer[0], MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);  
*/

  int nproc;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  if(nt>0)
  for (unsigned int i = 0; i < num_export; i++) 
  {
    printf("export_id:%i\n", export_local_id[i]);   
  }
  
  int **send_idx, *pivot; 
  REAL *tbuffer[3], *vbuffer;
  tbuffer[0] = (REAL *) malloc(nproc*size*3*sizeof(REAL));
  tbuffer[1] = (REAL *) malloc(nproc*size*3*sizeof(REAL));
  tbuffer[2] = (REAL *) malloc(nproc*size*3*sizeof(REAL)); 
  vbuffer = (REAL *) malloc(nproc*size*3*sizeof(REAL)); 
  
  send_idx = (int **) malloc(nproc*sizeof(int*));
  pivot = (int *) malloc(nproc*sizeof(int));
  for(int i=0;i<nproc;i++)
  {
    pivot[i] = 0;
    send_idx[i] = (int *) malloc(size*sizeof(int));
  }

  
  if(nt>0)
  for (unsigned int i = 0; i < num_export; i++) 
  {
    send_idx[rank[export_local_id[i]]][pivot[rank[export_local_id[i]]]] = export_local_id[i];
    pivot[rank[export_local_id[i]]]++;
  }
  
  int *pivot_buffer = (int *) malloc(nproc*sizeof(int));
  
  if(nt>0)
  for(int i=0;i<nproc;i++)
  {
    pivot_buffer[i] = pivot[i];
    for(unsigned int j=0;j<pivot[i];j++)
    {
      for(int k=0;k<3;k++)
      {
        tbuffer[0][(i*nproc*pivot[i]*3)+j+(k+1)] = t[0][k][send_idx[i][j]];  
        
        tbuffer[1][(i*nproc*pivot[i]*3)+j+(k+1)] = t[1][k][send_idx[i][j]]; 
        
        tbuffer[2][(i*nproc*pivot[i]*3)+j+(k+1)] = t[2][k][send_idx[i][j]]; 
        
        vbuffer[(i*nproc*pivot[i]*3)+j+(k+1)] = v[k][send_idx[i][j]];
      }
    }
  }
 
  //test  
  if(myrank==0)
  {
  //printf("velocity1Origin:%f\n", v[0][send_idx[1][0]]);
  //printf("velocity1:%f\n", vbuffer[(1*nproc*pivot[1]*3)+0+(0+1)]);
  //printf("t0Origin:%f\n", t[0][0][send_idx[1][0]]);  
  //printf("t0:%f\n", tbuffer[0][(1*nproc*pivot[1]*3)+0+(0+1)]);
  }

  //send buffers loop through procs is not correct, should loop through procs that need to send something/same for receive.(send has the pivot check*)
  for(int i=0;i<nproc;i++)
  {
    if(myrank != i && pivot[i] > 0)
    {  
      MPI_Send(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send(&tbuffer[0][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      MPI_Send(&tbuffer[1][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
      MPI_Send(&tbuffer[2][(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
      for(int x=0;x<pivot[i];x++)
      {
        printf("sent-tbuffer: %f\n", tbuffer[0][(i*nproc*pivot[i]*3)+x+(0+1)]);
      }
      MPI_Send(&vbuffer[(i*nproc*pivot[i]*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
    }     
  }

  for(int i=0;i<nproc;i++)
  {
    if(myrank != i && myrank != 0)
    {
      MPI_Recv(&pivot_buffer[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&tbuffer[0][(i*nproc*pivot_buffer[i]*3)], pivot_buffer[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&tbuffer[1][(i*nproc*pivot_buffer[i]*3)], pivot_buffer[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      MPI_Recv(&tbuffer[2][(i*nproc*pivot_buffer[i]*3)], pivot_buffer[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for(int x=0;x<pivot_buffer[i];x++)
      {
        printf("received-tbuffer: %f\n", tbuffer[0][(i*nproc*pivot_buffer[i]*3)+x+(0+1)]);
      }
      MPI_Recv(&vbuffer[(i*nproc*pivot_buffer[i]*3)], pivot_buffer[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }
   // MPI_Sendrecv(&pivot_buffer[0], 1, MPI_INT, 1, 0, MPI_COMM_WORLD, &pivot_buffer[0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status); 
   // MPI_Sendrecv(&tbuffer[0], pivot[0], MPI_DOUBLE, 1, 1, MPI_COMM_WORLD, &tbuffer[0], pivot_buffer[0], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
   // MPI_Sendrecv(&vbuffer[0], pivot[0], MPI_DOUBLE, 1, 2, MPI_COMM_WORLD, &vbuffer[0], pivot_buffer[0], MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);  

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
  unsigned int size; /* ranks size */
  int myrank;
  unsigned int *export_local_id;

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
    ERRMEM (export_local_id = (unsigned int *) malloc (sizeof(unsigned int[size])));
    
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
    ERRMEM (rank = (int *) malloc (sizeof(int[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (export_local_id = (unsigned int *) malloc (sizeof(unsigned int[size])));
  }
  
  int num_gid_entries, num_lid_entries, num_import, *import_procs, num_export, *export_procs;
  
  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  REAL step = 1E-3, time; unsigned int timesteps=0;
  
  //for (time = 0.0; time < 1.0; time += step)
  {
    
    loba_balance (lb, nt, t[0], tid, 1.1, rank, &num_gid_entries, &num_lid_entries, &num_import, 
                  import_procs, &num_export, export_procs, export_local_id); 
    
    for (unsigned int j = 0; j < nt; j++) printf ("rank[%u] = %d\n", j, rank[j]);
    
    migrate_triangles (myrank, size, nt, t, v, rank, num_gid_entries, num_lid_entries, 
                       num_import, import_procs, num_export, export_procs, export_local_id); 

    //contact_distance(nt, t, p, q, d); 

    //integrate_triangles (step, lo, hi, nt, t, v);
    timesteps++;
  }
  
  if(myrank == 0)
  {
   // write_pointsVTK(nt, t, v);
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
