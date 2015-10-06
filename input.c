#include "input.h" 

unsigned int load_pointsVTK(iREAL *t[3][3], unsigned int tid[], iREAL *mint, iREAL *maxt)
{
    FILE *fp1 = fopen("input/input.vtk", "r");
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
