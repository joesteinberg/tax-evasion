#ifndef __DIFF_EVO_UTILS_C__
#define __DIFF_EVO_UTILS_C__

//////////////////////////////////////////////////////////////////////////
// comments

/*
This source file contains some utility functions for the differential
evolution algorithm. I put them here to keep the main program file 
easy to read.
*/

//////////////////////////////////////////////////////////////////////////
// includes and macros

#include "diff_evo_utils.h"

//////////////////////////////////////////////////////////////////////////
// function definitions

void help()
{
  printf("\nUsage: diff_evo_mpi <x> <y>\n");
  printf("\t\tx=0: No evasion\n");
  printf("\t\tx=1: Evasion\n");
  printf("\t\ty: sensitivity analysis (0,1,2,3,4)\n");
}

int fcmp (const double x1, const double x2, const double epsilon)
{
  int exponent;
  double delta, difference;

  {
    double max = (fabs (x1) > fabs (x2)) ? x1 : x2;

    frexp (max, &exponent);
  }

  delta = ldexp (epsilon, exponent);
  difference = x1 - x2;

  if (difference > delta)       /* x1 > x2 */
    {
      return 1;
    }
  else if (difference < -delta) /* x1 < x2 */
    {
      return -1;
    }
  else                          /* -delta <= difference <= delta */
    {
      return 0;                 /* x1 ~=~ x2 */
    }
}

void assignd(int dim, double dest[], const double src[])
{
  int i;
  for(i=0; i<dim; i++)
    {
      dest[i] = src[i];
    }
}

void seed_rng(unsigned int seed)
{
  time_t now = time(NULL);
  unsigned char *p = (unsigned char *)&now;
  size_t i;

  for (i = 0; i < sizeof now; i++)
    seed = seed * (UCHAR_MAX + 2U) + p[i];

  srand(seed);
}

double rand_uniform(int bigrand)
{
	return bigrand * (1.0 / (RAND_MAX + 1.0) );
}

int rand_ind(const int lo, const int hi )
{
	return lo + rand_uniform( rand() ) * (hi - lo);
}

double rand_double(const double lo, const double hi)
{
	return lo + rand_uniform( rand() ) * (hi - lo);
}

// read diff evo parameters
int read_input_file
(
 const char *fname,
 int *strategy,
 int *genmax,
 int *refresh,
 int *dim,
 int *npop,
 double *F,
 double *CR,
 double *ubound,
 double *lbound
)
{
  FILE * infile = fopen(fname,"r");
  if(infile==NULL)
    {
      printf("\nCannot open input file!\n");
      return 1;
    }
  
  if(fscanf(infile,"%d",strategy) != 1)
    {
      printf("\nError reading strategy!\n");
      fclose(infile);
      return 1;
    }
  if(fscanf(infile,"%d",genmax) != 1)
    {
       printf("\nError reading genmax!\n");
      fclose(infile);
      return 1;
    }
  if(fscanf(infile,"%d",refresh) != 1)
    {
    printf("\nError reading refresh rate!\n");
      fclose(infile);
      return 1;
    }
  if(fscanf(infile,"%d",dim) != 1)
    {
        printf("\nError reading dim!\n");
      fclose(infile);
      return 1;
    }
  if(fscanf(infile,"%d",npop) != 1)
    {
    printf("\nError reading npop!\n");
      fclose(infile);
      return 1;
    }
  if(fscanf(infile,"%lf",F) != 1)
    {
    printf("\nError reading F!\n");
      fclose(infile);
      return 1;
    }
  if(fscanf(infile,"%lf",CR) != 1)
    {
      printf("\nError CR!\n");
      fclose(infile);
      return 1;
    }

  // check problem dimension before reading bounds
  if(*dim>MAXDIM)
    {
      printf("\nError in input file! dim = %d > MAXDIM = %d\n",*dim,MAXDIM);
      fclose(infile);
      return 1;
    }
  if(*dim<=0)
    {
      printf("\nError in input file! dim = %d, should be positive\n",*dim);
      fclose(infile);
      return 1;
    }


  // read bounds
  int i;
  for(i=0; i<*dim; i++)
    {
      if(fscanf(infile, "%lf %lf",lbound + i,ubound+i) != 2)
	{
    printf("\nError reading bounds!\n");
	  fclose(infile);
	  return 1;
	}

    }

  // close file and do the rest of the checking
  fclose(infile);

  // check population size
  if(*npop > MAXPOP)
    {
      printf("\nError in input file! npop = %d > MAXPOP = %d\n",*npop,MAXPOP);
      return 1;
    }
  if(*npop<=0)
    {
      printf("\nError in input file! npop = %d, should be positive\n",*npop);
      return 1;
    }

  // check genmax
  if(*genmax <= 0)
    {
      printf("\nError in input file! genmax = %d, should be positive\n",*genmax);
      return 1;
    }
  
  // check crossover rate
  if(*CR < 0 || *CR>1.0)
    {
      printf("\nError in input file! CR = %f, should be in [0,1]\n",*CR);
      return 1;
    }

  // check refresh rate
  if(*refresh <= 0)
    {
      printf("\nError in input file! refresh = %d, should be positive\n",*refresh);
      return 1;
    }

  // check strategy
  if(*strategy <= 0 || *strategy > 10)
    {
      printf("\nError in input file! strategy = %d, should be in {1,2,..,10}\n",*strategy);
      return 1;
    }

  // check bounds
  for(i=0; i<*dim; i++)
    {
      if(ubound[i]<lbound[i])
	{
	  printf("\nError in input file! ubound[%d] = %f < %f = lbound[%d]\n",i,ubound[i],lbound[i],i);
	  return 1;
	}
    }
  
  return 0;
  
}

#endif
