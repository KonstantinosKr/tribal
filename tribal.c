#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include "loba.h"
#include "error.h"
//#include "tribal_ispc.h"
#include "contact.h"


unsigned int load_pointsVTK(REAL *t[3][3], unsigned int tid[], REAL *mint, REAL *maxt)
{
    FILE *fp1 = fopen("input.vtk", "r");
    unsigned int nt;
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
            printf("found!\n");
            ch = fscanf(fp1,"%s",word);
            unsigned int n = atoi(word);
            //get points
            ch = fscanf(fp1,"%s",word);
            
            ERRMEM (point[0] = (double *)malloc (sizeof(double[n])));
            ERRMEM (point[1] = (double *)malloc (sizeof(double[n])));
            ERRMEM (point[2] = (double *)malloc (sizeof(double[n])));
            
            *mint = DBL_MAX;
            *maxt = DBL_MIN;
            for(unsigned int i=0;i<n;i++)
            {
                ch = fscanf(fp1, "%s", word);
                point[0][i] = atof(word);
                ch = fscanf(fp1, "%s", word);
                point[1][i] = atof(word);
                ch = fscanf(fp1, "%s", word);
                point[2][i] = atof(word);
                printf("POINT[0] = %f | POINT[1] = %f | POINT[2] = %f\n", point[0][i], point[1][i], point[2][i]);
                
                if(point[0][i] < *mint) 
                {
                  *mint = point[0][i];
                }

                if(point[1][i] < *mint) 
                {
                  *mint = point[1][i];
                }

                if(point[2][i] < *mint)
                {
                  *mint = point[2][i];
                }

                /////////////////////
                
                if(point[0][i] > *maxt) 
                {
                  *maxt = point[0][i];
                }

                if(point[1][i] > *maxt) 
                {
                  *maxt = point[1][i];
                }

                if(point[2][i] > *maxt)
                {
                  *maxt = point[2][i];
                }
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
                
                //printf("idx:%s T[0][0] = %f | T[0][1] = %f | T[0][2] = %f\n", word, t[0][0][i], t[0][1][i], t[0][2][i]);
                
                ch = fscanf(fp1,"%s",word);
                index = atoi(word);
                t[1][0][i] = point[0][index];
                t[1][1][i] = point[1][index];
                t[1][2][i] = point[2][index];
                
                //printf("idx:%s T[1][0] = %f | T[1][1] = %f | T[1][2] = %f\n", word, t[1][0][i], t[1][1][i], t[1][2][i]);
                
                ch = fscanf(fp1,"%s",word);
                index = atoi(word);
                t[2][0][i] = point[0][index];
                t[2][1][i] = point[1][index];
                t[2][2][i] = point[2][index];
                
                //printf("idx:%s T[2][0] = %f | T[2][1] = %f | T[2][2] = %f\n", word, t[2][0][i], t[2][1][i], t[2][2][i]);
                
                tid[i] = i;
            }
        }
    } while (ch != EOF);
    return nt;
}

void write_pointsVTK(int myrank, unsigned int nt, REAL *t[3][3], REAL *v[3], REAL lo[3], REAL hi[3], unsigned int timesteps)
{
    char iter[15];
    sprintf(iter, "%u_%i.vtk", timesteps, myrank);
    char filename[50] = "output/output"; //care or buffer overflow
    strcat(filename, iter);
    printf("%s\n", filename);
      
    FILE *fp = fopen(filename, "w+");
      
    fprintf(fp,"# vtk DataFile Version 2.0\nOutput vtk file\nASCII\n\nDATASET UNSTRUCTURED_GRID\nPOINTS %i float\n", (nt*3)+8);
      
    unsigned i;
    for(i = 0; i < nt; i++)
    {
      fprintf(fp,"%.5f %.5f %.5f\n%.5f %.5f %.5f\n%.5f %.5f %.5f\n", t[0][0][i], t[0][1][i], t[0][2][i], t[1][0][i], t[1][1][i], t[1][2][i], t[2][0][i], t[2][1][i], t[2][2][i]);
    }
 
    fprintf(fp, "%.5f %.5f %.5f\n%.5f %.5f %.5f\n%.5f %.5f %.5f\n%.5f %.5f %.5f\n", lo[0], lo[1], lo[2], lo[0], hi[1], lo[2], lo[0], hi[1], hi[2], lo[0], lo[1], hi[2]); 

    fprintf(fp, "%.5f %.5f %.5f\n%.5f %.5f %.5f\n%.5f %.5f %.5f\n%.5f %.5f %.5f\n", hi[0], hi[1], hi[2], hi[0], lo[1], hi[2], hi[0], lo[1], lo[2], hi[0], hi[1], lo[2]); 
    
    fprintf(fp,"\nCELLS %i %i\n", nt+12, nt+nt*3+(12*3));
    for(i = 0; i < nt*3; i=i+3)
    {
        fprintf(fp,"3 %i %i %i\n", i, i+1, i+2);
    }
    
    //left
    //1. lo,lo,lo
    //2. lo,hi,lo
    //3. lo,hi,hi
    //4. lo,lo,hi
    
    //link:
    //1. -> 2.
    //2. -> 3.
    //3. -> 4.
  
    //right
    //5. hi,hi,hi
    //6. hi,lo,hi
    //7. hi,lo,lo
    //8. hi,hi,lo
    
    //link:
    //5. -> 6.
    //6. -> 7.
    //7. -> 8.
   
    //interlink:
    //1. -> 6.
    //2. -> 5.
    //3. -> 7.
    //4. -> 8.
    //
    unsigned int ii = i;
    for(int j = 0; j < 3; j++)
    {
      fprintf(fp, "2 %i %i\n", i, i+1);
      fprintf(fp, "2 %i %i\n", i+4, i+4+1);
      i = i+1;
    }
    fprintf(fp, "2 %i %i\n", i, ii);
    fprintf(fp, "2 %i %i\n", i+4, ii+4);
    
    fprintf(fp, "2 %i %i\n", ii, ii+6);
    fprintf(fp, "2 %i %i\n", ii+1, ii+7);
    fprintf(fp, "2 %i %i\n", ii+2, ii+4);
    fprintf(fp, "2 %i %i\n", ii+3, ii+5);
    
    fprintf(fp,"\nCELL_TYPES %i\n", nt+12);
    for(i = 0; i<nt; i++)
    {
        fprintf(fp,"5\n");
    }
    for(int j = 0; j < 12; j++)
    {
      fprintf(fp, "3\n");
    }

    fclose(fp);
}


void normalize(unsigned int nt, REAL *t[3][3], REAL mint, REAL maxt) 
{  
  //range -255 to 255
  REAL inv_range = 510.0/(maxt-mint);

  for(unsigned int i=0; i < nt; i++)
  {
    for(int j=0; j < 3; j++)
    {
      for(int k=0; k < 3; k++)
      {
        t[j][k][i] = -255.0 + (t[j][k][i] - mint) * inv_range;
      }
    }
  }
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
  for(int i=0;i<nproc;i++)//n processes to prepare buffers for
  {
    for(unsigned int j=0;j<pivot[i];j++)//pivot gives n number of ids to loop through
    {
      for(int k=0;k<3;k++)//loop through the xyz axis
      {
        tbuffer[0][(i*size*3)+(j*3)+k] = t[0][k][send_idx[i][j]]; //point 0        
        tbuffer[1][(i*size*3)+(j*3)+k] = t[1][k][send_idx[i][j]]; //point 1
        tbuffer[2][(i*size*3)+(j*3)+k] = t[2][k][send_idx[i][j]]; //point 2

        ///printf("POSITION:%i\n\n", (j*3)+k);
        vbuffer[(i*size*3)+(j*3)+(k)] = v[k][send_idx[i][j]];
        pbuffer[(i*size*3)+(j*3)+(k)] = p[k][send_idx[i][j]];
        qbuffer[(i*size*3)+(j*3)+(k)] = q[k][send_idx[i][j]];
      }
        //printf("rank[%i] = tsend[0][0]:%f, tsend[0][1]:%f, tsend[0][2]:%f\n", i, t[0][0][send_idx[i][j]], t[0][1][send_idx[i][j]], t[0][2][send_idx[i][j]]);
        //printf("rank[%i] = tbuff[0][0]:%f, tbuff[0][1]:%f, tbuff[0][2]:%f\n", i, tbuffer[0][(i*size*3)+(j*3)+0], tbuffer[0][(i*size*3)+(j*3)+(1)], tbuffer[0][(i*size*3)+(j*3)+(2)]);
    }
  }

  /*for(int i=0;i<pivot[1];i++)
  printf("Processed - RANK[%i]: tid = %i\n t[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f\n t[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f\n t[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n", myrank, tid[i], 
              t[0][0][send_idx[1][i]], t[0][1][send_idx[1][i]], t[0][2][send_idx[1][i]], 
              t[1][0][send_idx[1][i]], t[1][1][send_idx[1][i]], t[1][2][send_idx[1][i]], 
              t[2][0][send_idx[1][i]], t[2][1][send_idx[1][i]], t[2][2][send_idx[1][i]]);

  for(int i=0;i<pivot[1];i++)
  {
      printf("ProcessedBUFFER - RANK[%i]: tid = %i\n t[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f\n t[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f\n t[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n", myrank, tid[i],
        tbuffer[0][(1*size*3)+(i*3)+(0)], tbuffer[0][(1*size*3)+(i*3)+(1)], tbuffer[0][(1*size*3)+(i*3)+(2)],
        tbuffer[1][(1*size*3)+(i*3)+(0)], tbuffer[1][(1*size*3)+(i*3)+(1)], tbuffer[1][(1*size*3)+(i*3)+(2)], 
        tbuffer[2][(1*size*3)+(i*3)+(0)], tbuffer[2][(1*size*3)+(i*3)+(1)], tbuffer[2][(1*size*3)+(i*3)+(2)]);
  } 
*/
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
 
  MPI_Request *myRequest;
  for(int x=0;x<num_export_unique;x++)
  {
    if(pivot[export_unique_procs[x]] > 0)
    {//safe check
      int i = export_unique_procs[x];
  
      //Asychronous Communication
      myRequest = malloc(sizeof(MPI_Request)); 
      MPI_Isend(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, myRequest);
      MPI_Wait(myRequest, MPI_STATUS_IGNORE);
      free(myRequest);
      myRequest = malloc(sizeof(MPI_Request)); 
      MPI_Isend(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, myRequest);
      MPI_Wait(myRequest, MPI_STATUS_IGNORE);
      free(myRequest);
      myRequest = malloc(sizeof(MPI_Request)); 
      MPI_Isend(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, myRequest);
      MPI_Wait(myRequest, MPI_STATUS_IGNORE);
      free(myRequest);
      myRequest = malloc(sizeof(MPI_Request)); 
      MPI_Isend(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, myRequest);  
      MPI_Wait(myRequest, MPI_STATUS_IGNORE);
      free(myRequest);
      myRequest = malloc(sizeof(MPI_Request)); 
      MPI_Isend(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, myRequest);
      MPI_Wait(myRequest, MPI_STATUS_IGNORE);
      free(myRequest);
      myRequest = malloc(sizeof(MPI_Request)); 
      MPI_Isend(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, myRequest);
      MPI_Wait(myRequest, MPI_STATUS_IGNORE);
      free(myRequest);
      myRequest = malloc(sizeof(MPI_Request)); 
      MPI_Isend(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, myRequest);
      MPI_Wait(myRequest, MPI_STATUS_IGNORE);
      free(myRequest);

    //BLOCKING COMM
/*      printf("RANK[%d]: send to rank %d\n", myrank, i);
      MPI_Send(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
      printf("send pivot: %i\n", pivot[i]);

      printf("pivot[%i]*3:%i\n", i, pivot[i]*3);
      printf("tbuffer[0][(i*size*3)]: %f\n", tbuffer[0][(i*size*3)]);
      MPI_Send(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      printf("send t1 points\n");
      MPI_Send(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
      printf("send t2 points\n");
      MPI_Send(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
      printf("send t3 points\n");
      
      MPI_Send(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
      printf("send v buffer\n");
      MPI_Send(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD);
      printf("send p buffer\n");
      MPI_Send(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD);
      printf("send q buffer\n");
  */  
      /*
      for(int j=0;j<pivot[i];j++) 
      {
        printf("SENT-BUFFER - RANK[%i]: tid = %i\n t[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f\n t[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f\n t[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n", myrank, tid[i],
        tbuffer[0][(i*size*3)+(j*3)+(0)], tbuffer[0][(i*size*3)+(j*3)+(1)], tbuffer[0][(i*size*3)+(j*3)+(2)],
        tbuffer[1][(i*size*3)+(j*3)+(0)], tbuffer[1][(i*size*3)+(j*3)+(1)], tbuffer[1][(i*size*3)+(j*3)+(2)], 
        tbuffer[2][(i*size*3)+(j*3)+(0)], tbuffer[2][(i*size*3)+(j*3)+(1)], tbuffer[2][(i*size*3)+(j*3)+(2)]);
      }
      */

    }
  }
  receive_idx = *nt-receive_idx; // set to last id
  
  for(int x=0;x<num_import_unique;x++)
  {
    int i = import_unique_procs[x];
    printf("RANK[%d]: receive from rank %d\n", myrank, i);
   
    myRequest = malloc(sizeof(MPI_Request)); 
    //non-blocking
    MPI_Irecv(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, myRequest);
    MPI_Wait(myRequest, MPI_STATUS_IGNORE);
    MPI_Irecv(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, myRequest);
    MPI_Wait(myRequest, MPI_STATUS_IGNORE);
    MPI_Irecv(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, myRequest);
    MPI_Wait(myRequest, MPI_STATUS_IGNORE);
    MPI_Irecv(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, myRequest);
    MPI_Wait(myRequest, MPI_STATUS_IGNORE);
    MPI_Irecv(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, myRequest);
    MPI_Wait(myRequest, MPI_STATUS_IGNORE);
    MPI_Irecv(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, myRequest);
    MPI_Wait(myRequest, MPI_STATUS_IGNORE);
    MPI_Irecv(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, myRequest);
    MPI_Wait(myRequest, MPI_STATUS_IGNORE);

    //BLOCKING COMM 
  /*  MPI_Recv(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&pbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&qbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    */
    /*for(int j=0;j<pivot[i];j++) 
    {
        printf("RECV-BUFFER - RANK[%i]: tid = %i\n t[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f\n t[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f\n t[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n", myrank, tid[i],
        tbuffer[0][(i*size*3)+(j*3)+(0)], tbuffer[0][(i*size*3)+(j*3)+(1)], tbuffer[0][(i*size*3)+(j*3)+(2)],
        tbuffer[1][(i*size*3)+(j*3)+(0)], tbuffer[1][(i*size*3)+(j*3)+(1)], tbuffer[1][(i*size*3)+(j*3)+(2)], 
        tbuffer[2][(i*size*3)+(j*3)+(0)], tbuffer[2][(i*size*3)+(j*3)+(1)], tbuffer[2][(i*size*3)+(j*3)+(2)]);
    }
*/

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

  for(int i=0; i<3;i++)
  {//free memory
    free(tbuffer[i]);
  }
    free(pivot);
    free(vbuffer);
    free(send_idx); 
}

void gen_velocities (REAL lo[3], REAL hi[3], unsigned int nt, REAL * v[3])
{
    for(unsigned int i = 0; i < nt; i++)
    {
      v[0][i] = hi[0] * drand48()*2;
      v[1][i] = hi[1] * drand48()*2;
      v[2][i] = hi[2] * drand48()*2;
    }
}

void integrate (REAL step, REAL lo[3], REAL hi[3], unsigned int nt, REAL * t[3][3], REAL * v[3])
{
    for(unsigned int i = 0; i < nt; i++)
    {
        t[0][0][i] += step*v[0][i];
        t[0][1][i] += step*v[1][i];
        t[0][2][i] += step*v[2][i];
        t[1][0][i] += step*v[0][i];
        t[1][1][i] += step*v[1][i];
        t[1][2][i] += step*v[2][i];
        t[2][0][i] += step*v[0][i];
        t[2][1][i] += step*v[1][i];
        t[2][2][i] += step*v[2][i];
        
        if (t[0][0][i] < lo[0]) v[0][i] *= -1;
        if (t[0][1][i] < lo[1]) v[1][i] *= -1;
        if (t[0][2][i] < lo[2]) v[2][i] *= -1;
        if (t[0][0][i] > hi[0]) v[0][i] *= -1;
        if (t[0][1][i] > hi[1]) v[1][i] *= -1;
        if (t[0][2][i] > hi[2]) v[2][i] *= -1;
        
        if (t[1][0][i] < lo[0]) v[0][i] *= -1;
        if (t[1][1][i] < lo[1]) v[1][i] *= -1;
        if (t[1][2][i] < lo[2]) v[2][i] *= -1;
        if (t[1][0][i] > hi[0]) v[0][i] *= -1;
        if (t[1][1][i] > hi[1]) v[1][i] *= -1;
        if (t[1][2][i] > hi[2]) v[2][i] *= -1;
        
        if (t[2][0][i] < lo[0]) v[0][i] *= -1;
        if (t[2][1][i] < lo[1]) v[1][i] *= -1;
        if (t[2][2][i] < lo[2]) v[2][i] *= -1;
        if (t[2][0][i] > hi[0]) v[0][i] *= -1;
        if (t[2][1][i] > hi[1]) v[1][i] *= -1;
        if (t[2][2][i] > hi[2]) v[2][i] *= -1;
    }
}

int main (int argc, char **argv)
{
  REAL *t[3][3]; /* triangles */
  REAL *v[3]; /* velocities */
  REAL *distance; /*distance */
  REAL *p[3],*q[3];
  unsigned int *tid; /* triangle identifiers */
  REAL lo[3] = {-255, -255, -255}; /* lower corner */
  REAL hi[3] = {255, 255, 255}; /* upper corner */
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
    else nt = 100000;

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
    ERRMEM (distance = (REAL *) malloc (sizeof(REAL[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX;
    
    //generate_triangles_and_velocities (lo, hi, nt, t, v, tid, pid);
    
    REAL mint, maxt;
    nt = load_pointsVTK(t, tid, &mint, &maxt); 
    normalize(nt, t, mint, maxt); 
   
    gen_velocities(lo, hi, nt, v);
  }
  else
  {
    /* set nt */
    nt = 0;

    /* buffers */
    if (argc > 1) size = atoi (argv[1])*4;
    else size = 100000*10;

    for (int i = 0; i < 3; i ++)
    {
      ERRMEM (t[0][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[1][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (t[2][i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (v[i] = (REAL *) malloc (sizeof(REAL[size])));
      
      ERRMEM (p[i] = (REAL *) malloc (sizeof(REAL[size])));
      ERRMEM (q[i] = (REAL *) malloc (sizeof(REAL[size])));
    }
    ERRMEM (distance = (REAL*) malloc (sizeof(REAL[size])));
    ERRMEM (tid = (unsigned int *) malloc (sizeof(unsigned int[size])));
    ERRMEM (pid = (unsigned int *) malloc (sizeof(unsigned int[size])));
      
    for(unsigned int i=0;i<size;i++) tid[i] = UINT_MAX;
  }
  
  int num_import, *import_procs, num_export, *export_procs;
  ZOLTAN_ID_PTR import_global_ids, import_local_ids, export_global_ids, export_local_ids;
  
  /* create load balancer */
  struct loba *lb = loba_create (ZOLTAN_RCB);

  /* perform time stepping */
  REAL step = 1E-3, time; unsigned int timesteps=0;
  
  //for (time = 0.0; time < 1.0; time += step)
  for(time = 0; time < 100; time++)
  {
    if(myrank == 0){printf("timestep: %i\n", timesteps);}

    loba_balance (lb, nt, t[0], tid, 1.1,
                  &num_import, &import_procs, 
                  &num_export, &export_procs, 
                  &import_global_ids, &import_local_ids, 
                  &export_global_ids, &export_local_ids);
      
    //printf("RANK[%i]:num_import:%d\n", myrank, num_import);
    //printf("RANK[%i]:num_export:%d\n", myrank, num_export);
   
    /*for(int i = 0; i < nt; i++)
      printf("Before - RANK[%i]: tid = %i\n t[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f\n t[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f\n t[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n", myrank, tid[i], 
              t[0][0][i], t[0][1][i], t[0][2][i], 
              t[1][0][i], t[1][1][i], t[1][2][i], 
              t[2][0][i], t[2][1][i], t[2][2][i]);*/
    
    migrate_triangles (size, &nt, t, v, p, q, tid, pid, 
                        num_import, import_procs, 
                        num_export, export_procs, 
                        import_global_ids, import_local_ids, 
                        export_global_ids, export_local_ids);
     
  /*for(int i = 0; i < nt; i++)
    printf("After - RANK[%i]: tid = %i\n t[0][0][i] = %f, t[0][1][i] = %f, t[0][2][i] = %f\n t[1][0][i] = %f, t[1][1][i] = %f, t[1][2][i] = %f\n t[2][0][i] = %f, t[2][1][i] = %f, t[2][2][i] = %f\n", myrank, tid[i], 
              t[0][0][i], t[0][1][i], t[0][2][i], 
              t[1][0][i], t[1][1][i], t[1][2][i], 
              t[2][0][i], t[2][1][i], t[2][2][i]);*/

/*  lo[0] = -100.0;
    lo[1] = -100.0;
    lo[2] = -100.0;

    hi[0] = 1000.0;
    hi[1] = 1000.0;
    hi[2] = 1000.0;

    loba_query(lb, 0, lo, hi, rank, &nprocs);*/
    
    //loba_getAdjacent(lb, myrank, neighborhood, &nNeighbors);
    //printf("RANK[%i]: AdjacentRanks:%i\n", myrank, nNeighbors);
    
    //for(int i=0;i<nNeighbors;i++){printf("RANK[%i]: Neighborhood[%i]:%i\n", myrank, i, neighborhood[i]);}
    
    //loba_migrateGhosts(lb, myrank, neighborhood, nNeighbors, size, &nt, t, v, p, q, tid, pid);

    //contact_distance(nt, t, p, q, distance); 
    
    integrate (step, lo, hi, nt, t, v);
    
    REAL mylo[3], myhi[3];
    loba_getbox (lb, myrank, mylo, myhi);
    
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
