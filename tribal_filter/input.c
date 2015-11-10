#include "input.h" 

void load_enviroment(int ptype[], unsigned int *nt, unsigned int nParticles, iREAL *t[3][3], unsigned int tid[], unsigned int pid[], iREAL *mint, iREAL *maxt)
{
  unsigned int n = 0;
  *nt = 0;
  for(unsigned int i = 0; i < nParticles; i++)
  {
    load_points(ptype[i], &n, i, *nt, t, tid, pid, mint, maxt);
    *nt = n + *nt;
    n = 0;
  }
}

/*
void save_enviroment()
{


}

void resize_enviroment()
{

}

void translate_enviroment()
{
  
}

void init_enviroment()
{


}
*/

void load_points(int ptype, unsigned int *nt, unsigned int bodyID, unsigned int startIDX, iREAL *t[3][3], unsigned int tid[], unsigned int pid[], iREAL *mint, iREAL *maxt)
{
  //////////VTK format////////////

  //Input Type
  //0: Triangulated Mesh
  //1: Triangle
  //2: Sphere
  //3: Square
  //4: Hexahedron
  
  iREAL min = DBL_MAX;
  iREAL max = DBL_MIN;
  FILE *fp1;
  if(ptype == 0)
  {
    char filename[100] = "input/mesh";
    char strtmp[100];
    sprintf(strtmp, "%i.vtk", bodyID);
    strcat(filename, strtmp);
    fp1 = fopen(filename, "r+");
    printf("%s\n", filename);
  } else if(ptype == 1)
  {
    fp1 = fopen("input/shapes/triangle.vtk", "r");
  } else if(ptype == 2)
  {
    fp1 = fopen("input/shapes/sphere.vtk", "r");
  } else if(ptype == 3)
  {
    fp1 = fopen("input/shapes/square.vtk", "r");
  } else if(ptype == 4)
  {
    fp1 = fopen("input/shapes/hexahedron.vtk", "r");
  }
  if( fp1 == NULL )
  {
      perror("Error while opening the file.\n");
      exit(EXIT_FAILURE);
  }
  
  char ch, word[500];
  iREAL *point[3];
  
  do {
      ch = fscanf(fp1,"%s",word);
      if(strcmp(word, "POINTS")==0)
      {
          printf("found!\n");
          ch = fscanf(fp1,"%s",word);
          unsigned int n = atol(word);
          //get points
          ch = fscanf(fp1,"%s",word);
          //printf("will read: %llu\n",n); 
          point[0] = (iREAL *)malloc (n*sizeof(iREAL));
          point[1] = (iREAL *)malloc (n*sizeof(iREAL));
          point[2] = (iREAL *)malloc (n*sizeof(iREAL));
          
          for(unsigned int i=0;i<n;i++)
          {
              fscanf(fp1, "%lf", &point[0][i]);
              fscanf(fp1, "%lf", &point[1][i]);
              fscanf(fp1, "%lf", &point[2][i]);
              //printf("POINT[0] = %f | POINT[1] = %f | POINT[2] = %f\n", point[0][i], point[1][i], point[2][i]);
              
              if(point[0][i] < min) 
              {
                min = point[0][i];
              }

              if(point[1][i] < min) 
              {
                min = point[1][i];
              }

              if(point[2][i] < min)
              {
                min = point[2][i];
              }

              /////////////////////
              
              if(point[0][i] > max) 
              {
                max = point[0][i];
              }

              if(point[1][i] > max) 
              {
                max = point[1][i];
              }

              if(point[2][i] > max)
              {
                max = point[2][i];
              }
          }
      }
      if(strcmp(word, "CELLS")==0)
      { 
          ch = fscanf(fp1,"%s",word);
          unsigned int n = atol(word);
          *nt = n;
          ch = fscanf(fp1,"%s",word);
          printf(":::%u::\n",n);
          for(unsigned int i=startIDX;i<startIDX+n;i++)
          {
              ch = fscanf(fp1,"%s",word);
              ch = fscanf(fp1,"%s",word);
              
              unsigned int index = atol(word);
              t[0][0][i] = point[0][index];
              t[0][1][i] = point[1][index];
              t[0][2][i] = point[2][index];
              
              //printf("idx:%s T[0][0] = %f | T[0][1] = %f | T[0][2] = %f\n", word, t[0][0][i], t[0][1][i], t[0][2][i]);
              
              ch = fscanf(fp1,"%s",word);
              index = atol(word);
              t[1][0][i] = point[0][index];
              t[1][1][i] = point[1][index];
              t[1][2][i] = point[2][index];
              
              //printf("idx:%s T[1][0] = %f | T[1][1] = %f | T[1][2] = %f\n", word, t[1][0][i], t[1][1][i], t[1][2][i]);
              
              ch = fscanf(fp1,"%s",word);
              index = atol(word);
              t[2][0][i] = point[0][index];
              t[2][1][i] = point[1][index];
              t[2][2][i] = point[2][index];
              
              //printf("idx:%s T[2][0] = %f | T[2][1] = %f | T[2][2] = %f\n", word, t[2][0][i], t[2][1][i], t[2][2][i]);
              
              tid[i] = i;
              pid[i] = bodyID;
          }
      }
  } while (ch != EOF);
  *mint = min;
  *maxt = max;
}

void normalize(unsigned int nt, iREAL *t[3][3], iREAL mint, iREAL maxt) 
{  
  //range -255 to 255
  iREAL inv_range = 510.0/(maxt-mint);

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