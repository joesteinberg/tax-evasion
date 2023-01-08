#ifndef __EXTERNS_C__
#define __EXTERNS_C__

#include "calibration.h"
#include "eqm.h"

opt_flag=0;

extern int setup()
{
  set_params();
  alloc_spline_mem();
  ss0.tax_flag=1;
  init_eqm_t(&ss0);
  if(load_eqm_txt_t(&ss0,"input/ss0"))
    {
      fprintf(logfile,"Failed to read benchmark equilibrium from file!\n");
      return 1;
    }
  else
    {
      return 0;
    }
}

extern void cleanup()
{
  free_spline_mem();
}

/*
objective function for calibration program 
1. double tauk: capital income tax
2. double taul: labor income tax
3. double abar: wealth tax threshold
*/

extern double cal_obj_fun(int D, double targeted_params[])
{
  fprintf(logfile,"\n----------------------------------------------------------------------------------------------\n");
  fprintf(logfile,"----------------------------------------------------------------------------------------------\n");
  fprintf(logfile,"----------------------------------------------------------------------------------------------\n");

  if(D==1)
    {
      tauk=0.0;
      tauv=0.0;
      abar=0.0;
      taul=targeted_params[0];
      fprintf(logfile,"\nSolving for equilibrium given taul = %0.6f\n\n",taul);
      fflush(logfile);
      snprintf(pref,128,"output/ss_opt1_%0.2f_%0.6f_%0.6f",
	       evasion_elast,taul,taua);
      snprintf(csv,128,"output/ss_opt1_%0.2f_%0.6f_%0.6f.csv",
	       evasion_elast,taul,taua);

    }
  else if(D==2)
    {
      tauk=0.0;
      tauv=0.0;
      taul=targeted_params[0];
      abar=targeted_params[1];
      fprintf(logfile,"\nSolving for equilibrium given (taul, abar) = (%0.6f,%0.6f)\n\n",taul,abar);
      fflush(logfile);
            snprintf(pref,128,"output/ss_opt2_%0.2f_%0.6f_%0.6f_%0.6f",
	       evasion_elast,taul,taua,abar);
	    snprintf(csv,128,"output/ss_opt2_%0.2f_%0.6f_%0.6f_%0.6f.csv",
	       evasion_elast,taul,taua,abar);

    }
    else if(D==3)
    {
      taul=targeted_params[0];
      abar=targeted_params[1];
      tauk= targeted_params[2];
      tauv = tauk;
      fprintf(logfile,"\nSolving for equilibrium given (tauk,taul, abar) = (%0.6f,%0.6f,%0.6f)\n\n",taul,abar);
      fflush(logfile);
      snprintf(pref,128,"output/ss_opt3_%0.2f_%0.6f_%0.6f_%0.6f_%0.6f",
	       evasion_elast,tauk,taul,taua,abar);
      snprintf(csv,128,"output/ss_opt3_%0.2f_%0.6f_%0.6f_%0.6f_%0.6f.csv",
	       evasion_elast,tauk,taul,taua,abar);
    }
  else
    {
      return GSL_NAN;
    }
  
  abar = M50_threshold*abar;
  
  if(test_dist(&ss0))
    {
      fprintf(logfile,"\nInvalid initial distribution! Sending NaN to master thread...\n");
      fflush(logfile);
      return GSL_NAN;
      //return +HUGE_VAL;
    }
  
  if(solve_ss(&ss0,0))
    {
      fprintf(logfile,"\nFailed to solve for equilibrium! Sending NaN to master thread...\n");
      fflush(logfile);
      return GSL_NAN;
      //return +HUGE_VAL;
    }
  
  fprintf(logfile,"Welfare: %0.16f\n", ss0.welfare);
  
  fprintf(logfile,"Storing equilibrium data...\n");
    
  write_eqm_t_csv(&ss0,csv);
  store_txt_eqm_t(&ss0,pref);

  fprintf(logfile,"\n----------------------------------------------------------------------------------------------\n");
  fprintf(logfile,"----------------------------------------------------------------------------------------------\n");
  fprintf(logfile,"----------------------------------------------------------------------------------------------\n");
      
  fflush(logfile);
  
  return -ss0.welfare_newborn;
}

#endif
