#ifndef __DIFF_EVO_UTILS_H__
#define __DIFF_EVO_UTILS_H__

//////////////////////////////////////////////////////////////////////////
// comments

/*
This header file declares some utility functions for the differential
evolution algorithm. I put them here to keep the main program file 
easy to read.
*/

//////////////////////////////////////////////////////////////////////////
// includes and macros

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include <memory.h>

#ifdef MPI
#include <mpi.h>
#endif

#define MAXPOP 500
//#define RESTART_POP_FROM_FILE 0
#define MAXDIM 35
#define EPS 1.0e-8

int read_pop_from_file;

//////////////////////////////////////////////////////////////////////////
// function definitions

void help();
int fcmp (const double x1, const double x2, const double epsilon);
void assignd(int dim, double dest[], const double src[]);
void seed_rng(unsigned int seed);
double rand_uniform(int bigrand);
int rand_ind(const int lo, const int hi );
double rand_double(const double lo, const double hi);
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
 double ubound[],
 double lbound[]
 );

#endif
