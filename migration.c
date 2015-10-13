#include "migration.h"

/* migrate triangles "in-place" to new ranks */
void migrate_triangles (unsigned long long int size, unsigned int *nt, iREAL *t[3][3], iREAL *v[3],
                              unsigned int *tid, unsigned int *pid,  
                              int num_import, int *import_procs, int num_export, int *export_procs, 
                              ZOLTAN_ID_PTR import_global_ids, ZOLTAN_ID_PTR import_local_ids,
                              ZOLTAN_ID_PTR export_global_ids, ZOLTAN_ID_PTR export_local_ids)
{
  int nproc, myrank;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank (MPI_COMM_WORLD, &myrank);
  
  //allocate memory for tmp buffers
  int **send_idx, *pivot, *pid_buffer, *rcvpid_buffer; 
  iREAL *tbuffer[3], *vbuffer;
  tbuffer[0] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  tbuffer[1] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  tbuffer[2] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL)); 
  vbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  pid_buffer = (int *) malloc(nproc*size*3*sizeof(int));

  int *rcvpivot;
  iREAL *trvbuffer[3], *vrvbuffer;
  trvbuffer[0] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  trvbuffer[1] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  trvbuffer[2] = (iREAL *) malloc(nproc*size*3*sizeof(iREAL)); 
  vrvbuffer = (iREAL *) malloc(nproc*size*3*sizeof(iREAL));
  rcvpid_buffer = (int *) malloc(nproc*size*3*sizeof(int));
  
  send_idx = (int **) malloc(nproc*sizeof(int*));
  pivot = (int *) malloc(nproc*sizeof(int));

  rcvpivot = (int *) malloc(nproc*sizeof(int));

  for(int i=0;i<nproc;i++)
  {
    rcvpivot[i] = 0;
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
      }
        pid_buffer[(i*size*3)+(j*3)] = pid[send_idx[i][j]];
    }
  }

  ///////////////////////////////////////////////
  //refine local arrays and ids (memory gaps)
  unsigned int pv=*nt-1;
  for(unsigned int i=0;i<num_export;i++)
  {//be cautious bug may be hidden here;
  
    for(unsigned int j=pv;j>export_local_ids[i];j--)//from last towards first but only until gap of exported
    {
      if(tid[j] != UINT_MAX)//if not marked as to be exported switch fill gaps
      {
        tid[export_local_ids[i]] = tid[j]; //send from 'last to first' the tids to 'first to last' in tid array
        tid[j] = UINT_MAX; //mark moved tid
        pid[export_local_ids[i]] = pid[j];
        pid[j] = UINT_MAX;
        for(int k=0;k<3;k++)
        {
          t[0][k][export_local_ids[i]] = t[0][k][j];
          t[1][k][export_local_ids[i]] = t[1][k][j];
          t[2][k][export_local_ids[i]] = t[2][k][j];
          
          v[k][export_local_ids[i]] = v[k][j];
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

  if(*nt > 0 && num_export > 0) 
  {
    receive_idx = *nt - num_export; //set to last id
  } else if(*nt >= 0 && num_export <= 0){
    receive_idx = *nt;
  }
  MPI_Request *myRequest = (MPI_Request*) malloc(num_export_unique*4*sizeof(MPI_Request));//4 sends
  MPI_Request *myrvRequest = (MPI_Request*) malloc(num_import_unique*4*sizeof(MPI_Request));//4 sends 
  
  for(int x=0;x<num_export_unique;x++)
  {
    int i = export_unique_procs[x];
    MPI_Send(&pivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD);
  }
  
  for(int x=0;x<num_import_unique;x++)
  {
    int i = import_unique_procs[x];  
    MPI_Recv(&rcvpivot[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  for(int x=0;x<num_import_unique;x++)
  {
    int i = import_unique_procs[x];  
    //printf("RANK[%d]: receive from rank %d\n", myrank, i);
   
    MPI_Irecv(&trvbuffer[0][(i*size*3)], rcvpivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myrvRequest[(x*4)+0]);
    MPI_Irecv(&trvbuffer[1][(i*size*3)], rcvpivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &myrvRequest[(x*4)+1]);
    MPI_Irecv(&trvbuffer[2][(i*size*3)], rcvpivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &myrvRequest[(x*4)+2]);
    MPI_Irecv(&vrvbuffer[(i*size*3)], rcvpivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &myrvRequest[(x*4)+3]);
  }

  
  for(int x=0;x<num_export_unique;x++)
  {
    int i = export_unique_procs[x];
  
    //printf("RANK[%d]: send to rank %d\n", myrank, i);
    MPI_Isend(&tbuffer[0][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &myRequest[(x*4)+0]);
    MPI_Isend(&tbuffer[1][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &myRequest[(x*4)+1]);
    MPI_Isend(&tbuffer[2][(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &myRequest[(x*4)+2]);  
    MPI_Isend(&vbuffer[(i*size*3)], pivot[i]*3, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &myRequest[(x*4)+3]);
  }
  
  for(int x=0;x<num_import_unique;x++)
  {
    MPI_Wait(&myrvRequest[(x*4)], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*4)+1], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*4)+2], MPI_STATUS_IGNORE);
    MPI_Wait(&myrvRequest[(x*4)+3], MPI_STATUS_IGNORE);
    int i = import_unique_procs[x];

    for(unsigned int j=0;j<rcvpivot[i];j++)
    {
      for(int k=0;k<3;k++)
      {
        t[0][k][receive_idx] = trvbuffer[0][(i*size*3)+(j*3)+(k)];        
        t[1][k][receive_idx] = trvbuffer[1][(i*size*3)+(j*3)+(k)]; 
        t[2][k][receive_idx] = trvbuffer[2][(i*size*3)+(j*3)+(k)]; 
        
        v[k][receive_idx] = vrvbuffer[(i*size*3)+(j*3)+(k)];
      }
      receive_idx++;
    }
  }
  
  for(int x=0;x<num_export_unique;x++)
  {
    MPI_Wait(&myRequest[(x*4)], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*4)+1], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*4)+2], MPI_STATUS_IGNORE);
    MPI_Wait(&myRequest[(x*4)+3], MPI_STATUS_IGNORE);
  }

  *nt = *nt + (num_import-num_export);
  
  for(int i=0; i<3;i++)
  {//free memory
    free(tbuffer[i]);
    free(trvbuffer[i]);
  }
    free(pivot);
    free(vbuffer);
    free(vrvbuffer);
    free(send_idx); 
    free(myRequest);
    free(myrvRequest);
}
