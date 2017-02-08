#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "header.h"

int get_seed ( );
void latin_random ( int dim_num, int point_num, int *seed, double x[] );
void latin_center ( int dim_num, int point_num, int *seed, double x[] );


int optimal_spacing()
{
  int i, j, ncubes=1,iseed;
  int DIM_NUM, POINT_NUM;
  FILE *fp;
  double xx[10*100];

  muh(0);

  if(ARGC>=3)
    {
      DIM_NUM = atoi(ARGV[2]);
      POINT_NUM = atoi(ARGV[3]);
    }
  //xx = calloc(sizeof(double),DIM_NUM*POINT_NUM);

  for(i=1;i<=ncubes;++i)
    {
      iseed = get_seed ( ); 
      muh(iseed);
      latin_random ( DIM_NUM, POINT_NUM, &iseed, xx );
      muh(1);
    }

  
  for(i=0;i<POINT_NUM;++i)
   {
      for(j=0;j<DIM_NUM;++j)
	printf("%f ",xx[i+j*DIM_NUM]);
      printf("\n");
    }
  
}
