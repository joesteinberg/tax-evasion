#ifndef __DIFF_EVO_MPI__
#define __DIFF_EVO_MPI__

//////////////////////////////////////////////////////////////////////////
// comments

/*
  This program implements an MPI-based parallel version of the 
  Differential Evolution optimization developed by Rainer Storn and
  Kenneth Price. More information on the algorithm can be found at 

  http://http.icsi.berkeley.edu/~storn/litera.html
*/

//////////////////////////////////////////////////////////////////////////
// includes and macros

#include "calibration.h"
#include "diff_evo_utils.h"
#include <omp.h>
#include <mpi.h>

double c[MAXPOP][MAXDIM], d[MAXPOP][MAXDIM];
double (*pold)[MAXPOP][MAXDIM], (*pnew)[MAXPOP][MAXDIM], (*pswap)[MAXPOP][MAXDIM];

//extern void load_eqm();
extern void setup(int evasion_flag, int evasion_elast_flag);
extern void cleanup();
double extern cal_obj_fun(int dim, double member[]); /* obj. funct. */
void evolve(int strategy, int npop, int dim, double F, double CR, const double bestit[], double tmp[][MAXDIM]);

void ensure_bounds
(
 const int dim,
 const double * lbnd,
 const double * hbnd,
 double member[]
 )
{
  unsigned int n;
  for(n=0; n<dim; n++)
    {
      if(member[n]<lbnd[n])
	{
	  member[n]=lbnd[n]+1.0e-10;
	}
      if(member[n]>hbnd[n])
	{
	  member[n]=hbnd[n]-1.0e-10;
	}
    }

  return;
}

//////////////////////////////////////////////////////////////////////////
// function source

// main function
int main(int argc, char *argv[])
{
  // non-MPI declarations
  char *strat[] = // strategy indicator
    {
      "",
      "DE/best/1/exp",
      "DE/rand/1/exp",
      "DE/rand-to-best",
      "DE/best/2/exp",
      "DE/rand/2/exp",
      "DE/best/1/bin",
      "DE/rand/1/bin",
      "DE/rand-to-best/1/bin",
      "DE/best/2/rand",
      "DE/rand/2/bin"
    };

  int read_pop_from_file=0;
  int i,j;
  int dim=0; // dimension of parameter vector
  int npop=0; // population size
  int strategy=0; // choice for screen output
  int refresh=1; // output refresh rate
  int genmax=0; // max number of generations
  //int seed; // seed value for rng
  int gen=0, imin=0; // gen number, index of pop member with best score

  long nfeval; /* number of function evaluations */

  double F=0.0, CR=0.0; // diff evo control variables
  double cvar=0.0; // cost variance
  double cmean=0.0; // cost mean
  double cmin=0.0; // ?
    
  double ubound[MAXDIM]; // upper bounds for parameters 
  double lbound[MAXDIM]; // lower bounds for parameters
  double tmp[MAXPOP][MAXDIM]; // population member array
  double best[MAXDIM]; // all-time best member
  double bestit[MAXDIM]; // best member of this iteration
  double cost[MAXPOP]; // members' objective function values
  double trial_cost[MAXPOP]; // buffer ?

  FILE *outfile;
  
  // MPI declarations
  int err; // error code for MPI function calls
  int id, ntasks, nslaves; // process id, number of processes, number of slave processes
  int source_id, dest_id; // message-sending process id numbers
  MPI_Status status;
  const int tag = 42;
  double cost_mpi; // message array
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int provided;
  int namelen;
  int nt=1, it=0;
  int evasion_flag=0, evasion_elast_flag=0;

  // initialize MPI ----------------------------------------------------
  err = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  if(err != MPI_SUCCESS)
    {
      printf("diff_evo_mpi.c: MPI initialization failed!\n");
      return 1;
    }

  // get number of tasks
  err = MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  if(err != MPI_SUCCESS)
    {
      printf("diff_evo_mpi.c: failed to get number of MPI tasks!\n");
      MPI_Finalize();
      return 1;
    }
  else if(ntasks<2)
    {
      printf("diff_evo_mpi.c: at least 2 processors required!\n");
      MPI_Finalize();
      return 1;
    }

  // get processor name
  err = MPI_Get_processor_name(processor_name, &namelen);
  if(err != MPI_SUCCESS)
    {
      printf("diff_evo_mpi.c: failed to get process name!\n");
      MPI_Finalize();
      return 1;
    }

  // get process id
  err = MPI_Comm_rank(MPI_COMM_WORLD, &id);
  if(err != MPI_SUCCESS)
    {
      printf("diff_evo_mpi.c: failed to get process id!\n");
      MPI_Finalize();
      return 1;
    }
  printf("MPI process %d of %d initialized.\n",id,ntasks);

  // check arguments
  if(argc != 3)
    {
      if(id==0)
	{
	  help();
	}
      MPI_Finalize();
      return 1;
    }

  if(strcmp(argv[1],"0")==0)
    {
      evasion_elast=0.0;
    }
  else if(strcmp(argv[1],"1")==0)
    {
      evasion_elast=0.5;
    }
  else if(strcmp(argv[1],"2")==0)
    {
      evasion_elast=2.5;
    }
  else if(strcmp(argv[1],"3")==0)
    {
      evasion_elast=7.5;
    }
  else
    {
      if(id==0)
	{
	  help();
	}
      help();
      MPI_Finalize();
      return 1;
    }

  if(strcmp(argv[2],"0")==0)
    {
      read_pop_from_file=0;
    }
  else if(strcmp(argv[2],"1")==0)
    {
      read_pop_from_file=1;
    }
  else
    {
      if(id==0)
	{
	  help();
	}
      MPI_Finalize();
      return 1;
    }

  if(id==0)
    {
      logfile = stdout;
    }
  else
    {
      open_log(id);
    }
 
  fprintf(logfile,"\n\n");
  fprintf(logfile,"Wealth tax optimization program (hybrid MPI-OpenMP parallelization)\n");
  fprintf(logfile,"Joseph Steinberg, University of Toronto\n\n");
  fprintf(logfile,"Evasion elasticity: %0.2f\n",evasion_elast);
  if(read_pop_from_file)
    {
      fprintf(logfile,"Restarting population from saved output from previous run\n");
    }
  else
    {
      fprintf(logfile,"Initializing fresh population\n");
    }
  fprintf(logfile,"\nMPI process number %d\n",id);  
  fprintf(logfile,"OpenMP parallel processing with %d threads\n",omp_get_num_threads());
  
  // all processes need to set up income process, etc.
  // read differential evolution parameters
  gsl_set_error_handler_off();
  seed_rng(id);

  if(read_input_file
     (
      "input/opt_params.txt",
      &strategy,
      &genmax,
      &refresh,
      &dim,
      &npop,
      &F,
      &CR,
      ubound,
      lbound
      ))
    {
      printf("\nFailure reading differential evolution parameters!\n");
      MPI_Finalize();
      return 1;
    }

  if(dim==1)
    {
      fprintf(logfile,"Optimizing over labor and wealth taxes. Zero capital income tax and no wealth tax threshold.\n");
    }
  else if(dim==2)
    {
      fprintf(logfile,"Optimizing over labor tax, wealth tax, and wealth tax threshold. Zero capital income tax.\n");

    }
    else if(dim==3)
    {
      fprintf(logfile,"Optimizing over labor tax, capital incomne tax, wealth tax, and wealth tax threshold.\n");
    }


  // number of slave processes
  nslaves = ntasks-1;
  if(ntasks != npop)
    {
      printf("You must use exactly %d processes\n",npop);
      MPI_Finalize();
      return 1;
    }
	  
  // slave process -----------------------------------------------------------------
  // basically wait for messages from the master thread...
  // when you receive one, evaluate the objective function using the new member data and report back
  /// then wait some more
  if(id>0)
    {
      fprintf(logfile,"Slave process waiting for tasks from master...\n\n");
      setup();
      while(1)
	{
	  i=0;
	  err = MPI_Recv(tmp[i], dim, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
	  if(isinf(tmp[0][0]))
	    {
	      fprintf(logfile,"Closing slave thread %d\n",id);
	      cleanup();
	      close_log();
	      MPI_Finalize();
	      return 0;
	    }
	  cost_mpi = cal_obj_fun(dim,tmp[0]);
	  dest_id = 0;
	  err = MPI_Send(&cost_mpi, 1, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
	}
    }

  // master process --------------------------------------------------------------
  
  // start timer
  fprintf(logfile,"Master process starting...\n\n");
  if(setup())
    {
      fprintf(logfile,"Failed to setup model economy!\n")
      for (j=0; j<nslaves; j++)
	{
	  i=0;
	  tmp[i][0] = HUGE_VAL;
	  dest_id=j+1;
	  err = MPI_Send(tmp[i], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
	}
      fprintf(logfile,"Closing master thread %d\n",id);
      close_log();
      MPI_Finalize();
      return 1;
    }

  time_t start, stop;
  time(&start);

  // open output files
  outfile = fopen("output/diff_evo_output.txt","w");
  if(outfile == NULL)
    {
      fprintf(logfile,"\nCannot open output file\n");
      for (j=0; j<nslaves; j++)
	{
	  i=0;
	  tmp[i][0] = HUGE_VAL;
	  dest_id=j+1;
	  err = MPI_Send(tmp[i], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
	}
      fprintf(logfile,"Closing master thread %d\n",id);
      close_log();
      MPI_Finalize();
      return 1;
    }
  
  nfeval = 0;

  // initialize the population
  if(read_pop_from_file)
    {
      FILE * popfile = fopen("output/popfile.txt","rb");
      if(popfile != NULL)
	{
	  int readcnt=0;
	  for(i=0; i<npop; i++)
	    {
	      for(j=0; j<dim; j++)
		{
		  readcnt += fscanf(popfile,"%lf ",&(c[i][j]));
		}
	      readcnt += fscanf(popfile,"%lf\n",&(cost[i]));
	    }
	  fclose(popfile);

	  if(readcnt != (npop)*(dim+1))
	    {
	      fprintf(logfile,"\nCannot read population file\n");
	      for (j=0; j<nslaves; j++)
		{
		  i=0;
		  tmp[i][0] = HUGE_VAL;
		  dest_id=j+1;
		  err = MPI_Send(tmp[i], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
		}
	      fprintf(logfile,"Closing master thread %d\n",id);
	      close_log();
	      MPI_Finalize();
	      return 1;
	    }
	}
      else
	{
	  fprintf(logfile,"\nCannot open population file\n");
	  for (j=0; j<nslaves; j++)
	    {
	      i=0;
	      tmp[i][0] = HUGE_VAL;
	      dest_id=j+1;
	      err = MPI_Send(tmp[i], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
	    }
	  fprintf(logfile,"Closing master thread %d\n",id);
	  close_log();
	  MPI_Finalize();
	  return 1;
	}
    }
  else
    {
      fprintf(logfile,"Randomly initializing population\n");
      for(i=0; i<npop; i++)
	{
	  for(j=0; j<dim; j++)
	    {
	      c[i][j] = lbound[j] + rand_uniform(rand())*(ubound[j]-lbound[j]);
	    }
	}
    }

  if(!read_pop_from_file)
    {
      // do the first round of work
      fprintf(logfile,"Calculating objective function values for initial population\n");

      // first send all the slaves their new member data and ask them to compute energies
      //for(i=0; i< (int)((npop/nslaves)+1); i++)
      for(j=0; j<nslaves; j++)
	{
	  dest_id = j+1;
	  nfeval++;
	  err = MPI_Send(c[j+1],dim,MPI_DOUBLE,dest_id, tag, MPI_COMM_WORLD);
	}

      // evaluate the first member of the population using the master thread
      fprintf(logfile,"Master thread performing slave task while waiting for other slaves...\n\n");
      cost[0] = cal_obj_fun(dim,c[0]);
  
      // now receive their responses
      for(j=0; j<nslaves; j++)
	{
	  err = MPI_Recv(&cost_mpi, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
	  source_id = (status.MPI_SOURCE);
	  cost[source_id]=cost_mpi;
	}
    }

  // figure out which member is the best and save it
  cmin = cost[0];
  imin = 0;
  for(i=0; i<npop; i++)
    {
      if(isnan(cost[i]))
	{
	  fprintf(logfile, "NaN value received from process %d!\n",i);
	  for (j=0; j<nslaves; j++)
	    {
	      i=0;
	      tmp[i][0] = HUGE_VAL;
	      dest_id=j+1; /* Destination address */
	      err = MPI_Send(tmp[i], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
	    }
	  fprintf(logfile,"Closing master thread %d\n",id);
	  close_log();
	  MPI_Finalize();
	  return 1;
	}
      if(cost[i]<cmin)
	{
	  cmin=cost[i];
	  imin=i;
	}
    }

  assignd(dim,best,c[imin]);
  assignd(dim,bestit,c[imin]); 
    

  pold = &c;
  pnew = &d;

  // main iteration loop
  gen = 0;
  while(gen<genmax)
    {
      fprintf(logfile,"Generation %d starting\n",gen);
      fflush(logfile);
      gen++;
      imin=0;

      // first, evolve!
      evolve(strategy, npop, dim, F, CR, bestit, tmp);
      for(i=0; i<npop; i++)
	{
	  ensure_bounds(dim,lbound,ubound,tmp[i]);
	}

      // send out the new member data
      for (j=0; j<nslaves; j++)
	{
	  dest_id=j+1;
	  nfeval++;
	  err = MPI_Send(tmp[j+1], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
	}

      // evaluate the first member of the population using the master thread
      fprintf(logfile,"Master thread performing slave task while waiting for other slaves...\n\n");
      trial_cost[0] = cal_obj_fun(dim,tmp[0]);

      // receive the responses
      for (j=0; j<nslaves; j++)
	{
	  err = MPI_Recv(&cost_mpi, 1, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
	  source_id = (status.MPI_SOURCE);
	  trial_cost[source_id]=cost_mpi;
	}
      
      // evaluate the responses
      for (i=0; i<npop; i++)
	{
	  if(isnan(trial_cost[i]))
	    {
	      fprintf(logfile, "NaN value received from process %d!\n",i);
	      for (j=0; j<nslaves; j++)
		{
		  i=0;
		  tmp[i][0] = HUGE_VAL;
		  dest_id=j+1; /* Destination address */
		  err = MPI_Send(tmp[i], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
		}
	      fprintf(logfile,"Closing master thread %d\n",id);
	      close_log();
	      MPI_Finalize();
	      return 1;
	    }
	  if (trial_cost[i] <= cost[i])
	    {
	      cost[i]=trial_cost[i];
	      assignd(dim,(*pnew)[i],tmp[i]);
	      if (trial_cost[i]<cmin)
		{                    
		  cmin=trial_cost[i];
		  imin=i;
		  assignd(dim,best,tmp[i]);
		}
	    }
	  else
	    {
	      assignd(dim,(*pnew)[i],(*pold)[i]);
	    }
	}

      assignd(dim,bestit,best);

      pswap = pold;
      pold  = pnew;
      pnew  = pswap;

      cmean = 0.;
      for (j=0; j<npop; j++)
	{
	  cmean += cost[j];
	}
      cmean = cmean/npop;

      cvar = 0.;
      for (j=0; j<npop; j++)
	{
	  cvar += (cost[j] - cmean)*(cost[j] - cmean);
	}
      cvar = cvar/(npop-1);

      // output
      //if (gen>1 && gen%refresh==1)
      if(gen>=1)
	{
	  fprintf(logfile,"\n\n\tBest-so-far cost funct. value=%-15.10g\n",cmin);

	  for (j=0;j<dim;j++)
	    {
	      fprintf(logfile,"\n\tbest[%d]=%-15.10g",j,best[j]);
	    }
	  fprintf(logfile,"\n\n\tGeneration=%d  NFEs=%ld   Strategy: %s    ",gen,nfeval,strat[strategy]);
	  fprintf(logfile,"\n\tNP=%d    F=%-4.2g    CR=%-4.2g   cost-variance=%-10.5g\n",
		  npop,F,CR,cvar);
	}

      fprintf(outfile,"%ld   %-15.10g\n",nfeval,cmin);
      fflush(outfile);

      FILE * popfile = fopen("output/popfile.txt","wb");
      for(i=0; i<npop; i++)
	{
	  for(j=0; j<dim; j++)
	    {
	      fprintf(popfile,"%0.16f ",(*pnew)[i][j]);
	    }
	  fprintf(popfile,"%0.16f\n",cost[i]);
	}

      fclose(popfile);
    }

  // ---------------------------------------------------------------
  // wrap everything up

  // log output
  fprintf(logfile,"\n\n\tTerminal best obj. funct. value = %-15.10g\n",cmin);
  for (j=0;j<dim;j++)
    {
      fprintf(logfile,"\n\tbest[%d]=%-15.10g",j,best[j]);
    }
  fprintf(logfile,"\n\tGeneration=%d\tNFEs=%ld\tStrategy: %s    ",gen,nfeval,strat[strategy]);
  fprintf(logfile,"\n\tNP=%d\tF=%-4.2g\tCR=%-4.2g\tCost-variance=%-10.5g\n",npop,F,CR,cvar);

  // file output
  fprintf(outfile,"\n\n\nTerminal best obj. funct. value = %-15.10g\n",cmin);
  for (j=0;j<dim;j++)
    {
      fprintf(outfile,"\nbest[%d]=%-15.10g",j,best[j]);
    }
  fprintf(outfile,"\n\nGeneration=%d\tNFEs=%ld\tStrategy: %s    ",gen,nfeval,strat[strategy]);
  fprintf(outfile,"\n NP=%d\tF=%-4.2g\tCR=%-4.2g\tCost-variance=%-10.5g\n",npop,F,CR,cvar);
  fclose(outfile);

  time(&stop);
  fprintf(logfile,"\nTime: %0.2f seconds\n\n",difftime(stop,start));

  // tell the slave threads to exit
  for (j=0; j<nslaves; j++)
    {
      i=0;
      tmp[i][0] = HUGE_VAL;
      dest_id=j+1; /* Destination address */
      err = MPI_Send(tmp[i], dim, MPI_DOUBLE, dest_id, tag, MPI_COMM_WORLD);
    }
  
  fprintf(logfile,"Closing master thread %d\n",id);
  close_log();
  MPI_Finalize();
  
  return 0;
  
}

//////////////////////////////////////////////////////////////////////////
// differential evolution

void evolve(int strategy, int npop, int dim, double F, double CR, const double bestit[], double tmp[][MAXDIM])
{
  int i, r1, r2, r3, r4, r5, n, L;

  // loop through the population...
  for(i=0; i<npop; i++)
    {
      // first find 5 other random members, then do the mutation
      do
	{
	  r1=(int)(rand_uniform(rand())*npop);
	}
      while(r1==i);
      do
	{
	  r2=(int)(rand_uniform(rand())*npop);
	}
      while(r2==i || r2==r1);
      do
	{
	  r3=(int)(rand_uniform(rand())*npop);
	}
      while(r3==i || r3==r1 || r3==r2);
      do
	{
	  r4=(int)(rand_uniform(rand())*npop);
	}
      while(r4==i || r4==r1 || r4==r2 || r4==r3);
      do
	{
	  r5=(int)(rand_uniform(rand())*npop);
	}
      while(r5==i || r5==r1 || r5==r2 || r5==r3 || r5==r4);

      // DE/best/1/exp
      if (strategy == 1)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  L = 0;
	  do
	    {
	      tmp[i][n] = bestit[n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
	      n = (n+1)%dim;
	      L++;
	    }
	  while((rand_uniform(rand()) < CR) && (L < dim));
	}

      // DE/rand/1/exp
      else if (strategy == 2)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  L = 0;
	  do
	    {
	      tmp[i][n] = (*pold)[r1][n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
	      n = (n+1)%dim;
	      L++;
	    }while((rand_uniform(rand()) < CR) && (L < dim));
	}

      // DE/rand-to-best/1/exp
      else if (strategy == 3)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  L = 0;
	  do
	    {
	      tmp[i][n] = tmp[i][n] + F*(bestit[n] - tmp[i][n]) + F*((*pold)[r1][n]-(*pold)[r2][n]);
	      n = (n+1)%dim;
	      L++;
	    }while((rand_uniform(rand()) < CR) && (L < dim));
	}

      // DE/best/2/exp
      else if (strategy == 4)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  L = 0;
	  do
	    {
	      tmp[i][n] = bestit[n] +
		((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
	      n = (n+1)%dim;
	      L++;
	    }while((rand_uniform(rand()) < CR) && (L < dim));
	}

      // DE/rand/2/exp
      else if (strategy == 5)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  L = 0;
	  do
	    {
	      tmp[i][n] = (*pold)[r5][n] +
		((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
	      n = (n+1)%dim;
	      L++;
	    }while((rand_uniform(rand()) < CR) && (L < dim));
	}
      
      // DE/best/1/bin
      else if (strategy == 6)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  for (L=0; L<dim; L++)
	    {
	      if ((rand_uniform(rand()) < CR) || L == (dim-1))
		{
		  tmp[i][n] = bestit[n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
		}
	      n = (n+1)%dim;
	    }
	}


      // DE/rand/1/bin
      else if (strategy == 7)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  for (L=0; L<dim; L++)
	    {
	      if ((rand_uniform(rand()) < CR) || L == (dim-1))
		{
		  tmp[i][n] = (*pold)[r1][n] + F*((*pold)[r2][n]-(*pold)[r3][n]);
		}
	      n = (n+1)%dim;
	    }
	}
      
      // DE/rand-to-best/1/bin
      else if (strategy == 8)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  for (L=0; L<dim; L++)
	    {
	      if ((rand_uniform(rand()) < CR) || L == (dim-1))
		{
		  tmp[i][n] = tmp[i][n] + F*(bestit[n] - tmp[i][n]) + F*((*pold)[r1][n]-(*pold)[r2][n]);
		}
	      n = (n+1)%dim;
	    }
	}
      
      // DE/best/2/bin
      else if (strategy == 9)
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  for (L=0; L<dim; L++)
	    {
	      if ((rand_uniform(rand()) < CR) || L == (dim-1))
		{
		  tmp[i][n] = bestit[n] +
		    ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
		}
	      n = (n+1)%dim;
	    }
	}
      
      // DE/rand/2/bin
      else
	{
	  assignd(dim,tmp[i],(*pold)[i]);
	  n = (int)(rand_uniform(rand())*dim);
	  for (L=0; L<dim; L++)
	    {
	      if ((rand_uniform(rand()) < CR) || L == (dim-1))
		{
		  tmp[i][n] = (*pold)[r5][n] +
		    ((*pold)[r1][n]+(*pold)[r2][n]-(*pold)[r3][n]-(*pold)[r4][n])*F;
		}
	      n = (n+1)%dim;
	    }
	}
    } 
}

#endif
