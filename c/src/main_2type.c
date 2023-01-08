#define EXTERN
#define EXTERN2

#ifndef __MAIN_C__
#define __MAIN_C__

#include "calibration.h"
#include "eqm_2type.h"

void help()
{
  fprintf(logfile,"Bad user input!\n");

  fprintf(logfile,"Usage:\n");
  
  fprintf(logfile,"model_2type 0: steady state with current US policies\n");
  fprintf(logfile,"model_2type 1: short-run partial equilibrium with new tax system\n");
  fprintf(logfile,"model_2type 2: Warren wealrh tax\n");  
  fprintf(logfile,"model_2type 3 <c> <d> <e>: series of capital income tax reforms\n");
  fprintf(logfile,"\td: max tax hike\n");
  fprintf(logfile,"\te: number of reforms in each direction\n");
  fprintf(logfile,"\tf: short run partial equilibrium\n");
  fprintf(logfile,"model_2type 4 <c> <d>: series of wealth tax reforms\n");
  fprintf(logfile,"\tc: max wealth tax\n");
  fprintf(logfile,"\td: number of reforms\n");
  fprintf(logfile,"\te: short run partial equilibrium\n");
 }

int main(int argc, char * argv[])
{
  time_t start, stop;
  time(&start);
  write_binary=1;
  evasion_type=1;
  logfile = stdout;
  linebreak();
  fprintf(logfile,"Tax Evasion and Capital Taxation\nShahar Rotberg and Joseph Steinberg\n\n");
  fprintf(logfile,"Model with 2 Agent Types (Evaders and Non-Evaders)\n\n");

#ifdef _OPENMP
  omp_set_num_threads(NTH);
  fprintf(logfile,"Parallel processing with %d threads\n",NTH);
#else
  fprintf(logfile,"No parallelization\n");
#endif
  
  alloc_mem();

  //ss0->frac_v = 0.625;
  //ss1->frac_v = 0.625;

  ss0->frac_v = 0.35;
  ss1->frac_v = 0.35;

  //ss0->frac_v = 0.070;
  //ss0->frac_v = 0.070;
  
  fprintf(logfile,"Evader fraction = %0.3f\n\n",ss0->frac_v);

  constraint_flag=0;
  detection_revenue_flag=0;
  benchmark_eqm = 1;
  taua_clear_gbc=0;
  taul_clear_gbc=0;
  wealth_tax_type=0;
  lump_sum=0;
  public_goods=0;
  taul_flag=0;
  
  
  //---------------------------------------------------------------------------
  // Initial steady state
  if(argc==2 && strcmp(argv[1],"0")==0)
    {
      fprintf(logfile,"Solving for initial steady state with current US tax code\n");
      set_params(1);     
      dist_tol = 5.0e-10;
      wage_tol = 3e-4;
      init_eqm_t(ss0);
      init_eqm_t(ss1);
	    
      snprintf(fname0,128,"output/ss0_2type_%0.3f",ss0->frac_v);
      snprintf(csv,128,"output/ss0_2type_%0.3f.csv",ss0->frac_v);

      //load_benchmark_eqm_t(ss0);
      load_bin_eqm_t(ss0,"output/ss0_2type_0.075.bin");
      verbose=1;
      benchmark_eqm=1;
      solve_ss(ss0,0);
      store_bin_eqm_t(ss0,fname0);
      copy_eqm_t(ss1,ss0);
      
      free_mem();
    }
  
  //---------------------------------------------------------------------------
  // short run elasticity to 1.5% tax
  else if(argc==2 && strcmp(argv[1],"1")==0)
    {
      fprintf(logfile,"Solving for short-run partial equilibrium with 1pct flat wealth tax\n");

      set_params(1);
      dist_tol = 5.0e-10;
      wage_tol = 3e-4;
      init_eqm_t(ss0);
      init_eqm_t(ss1);

      snprintf(fname0,128,"output/ss0_2type_%0.3f.bin",ss0->frac_v);
      snprintf(csv,128,"output/ss1_srpe_2type_%0.3f.csv",ss0->frac_v);
      load_bin_eqm_t(ss0,fname0);
      
      load_bin_eqm_t(ss1,fname0);
      taua = 0.015;
      abar=0;
      verbose=1;
      dist_max_iter=1;
      wage_max_iter=1;
      solve_ss(ss1,1);
      calc_wealth_tax_evasion_elast(1,ss1);
      
      free_mem();
    }
  //---------------------------------------------------------------------------
  // Warren wealth tax
  else if(argc==2 && strcmp(argv[1],"2")==0)
    {
      
      fprintf(logfile,"Solving for Warren tax on wealth above 50M\n");

      set_params(1);
      //dist_tol = 5.0e-10;
      //wage_tol = 3e-4;
      init_eqm_t(ss0);
      init_eqm_t(ss1);

      snprintf(fname0,128,"output/ss0_2type_%0.3f.bin",ss0->frac_v);
      load_bin_eqm_t(ss0,fname0);
      load_bin_eqm_t(ss1,fname0);

      taua_clear_gbc=0;
      wealth_tax_type=1;
      lump_sum=1;
      public_goods=0;
      set_wealth_tax_params(ss0,wealth_tax_type);
      benchmark_eqm=0;
      verbose=1;
      snprintf(csv,128,"output/ss1_warren_2type_%0.3f.csv",ss0->frac_v);
      snprintf(pref1,128,"output/ss1_warren_2type_%0.3f",ss0->frac_v);
      solve_ss(ss1,0);

      calc_wealth_tax_evasion_elast(1,ss1);

      free_mem();
    }

  //---------------------------------------------------------------------------
  // fprintf(logfile,"model 3 <d> <e> <f> series of capital income tax reforms\n");
  // fprintf(logfile,"\td: max tax hike\n");
  // fprintf(logfile,"\te: number of reforms\n");
  // fprintf(logfile,"\tf: short-run partial equilibrium\n");
  else if(argc==5 && strcmp(argv[1],"3")==0)
    {
      fprintf(logfile,"Solving for series of steady states with different capital income tax rates\n");

      set_params(1);

      // set the file to load initial eqm
      dist_tol = 5.0e-10;
      wage_tol = 3e-4;

      init_eqm_t(ss0);
      init_eqm_t(ss1);

      snprintf(fname0,128,"output/ss0_2type_%0.3f.bin",ss0->frac_v);
      load_bin_eqm_t(ss0,fname0);
      load_bin_eqm_t(ss1,fname0);

      // set up the exercise
      taua_clear_gbc=0;
      lump_sum=1;
      public_goods=0;
      benchmark_eqm=0;
      verbose=1;

      double tauk0=tauk;
      double tauk_hi=tauk0 + strtof(argv[2],NULL);
      int n_tau = atoi(argv[3]);
      double d = (tauk_hi-tauk0)/(n_tau-1);

      int pe = 0;
      if(strcmp(argv[4],"1")==0)
	{
	  pe=1;
	  dist_max_iter=1;
	  wage_max_iter=1;
	  fprintf(logfile,"Short-run partial equilibrium\n");
	}

      
      for(int n=0; n<n_tau; n++)
	{
	  tauk=tauk0+(n+1)*d;
	  
	  linebreak();
	  fprintf(logfile,"Solving for tauk = %0.6f\n",tauk);

	  if(pe==0)
	    {
	      snprintf(csv,128,"output/ss1_tauk_%0.6f_2type_%0.3f.csv",tauk,ss0->frac_v);
	      snprintf(pref1,128,"output/ss1_tauk_%0.6f_2type_%0.3f",tauk,ss0->frac_v);
	    }
	  else
	    {
	      snprintf(csv,128,"output/ss1_tauk_%0.6f_2type_%0.3f_pe.csv",tauk,ss0->frac_v);
	      snprintf(pref1,128,"output/ss1_tauk_%0.6f_2type_%0.3f_pe",tauk,ss0->frac_v);
	    }
	  solve_ss(ss1,pe);
	}
      
      free_mem();
    }
  //---------------------------------------------------------------------------
  // fprintf(logfile,"model 4 <c> <d> <e>: series of wealth tax reforms\n");
  // fprintf(logfile,"\tc: max wealth tax\n");
  // fprintf(logfile,"\td: number of reforms\n");
  // fprintf(logfile,"\te: short-run partial equilibrium\n");
  else if(argc==5 && strcmp(argv[1],"4")==0)
    {
      fprintf(logfile,"Solving for series of steady states with different flat wealth tax rates\n");
      
      set_params(1);

      dist_tol = 5.0e-10;
      wage_tol = 3e-4;
      
      // parse the additional arguments
      int pe=0;

      // set the file to load initial eqm
      init_eqm_t(ss0);
      init_eqm_t(ss1);

      snprintf(fname0,128,"output/ss0_2type_%0.3f.bin",ss0->frac_v);
      load_bin_eqm_t(ss0,fname0);
      load_bin_eqm_t(ss1,fname0);

      // set up the exercise
      verbose=1;
      taua_clear_gbc=0;
      taul_clear_gbc=0;
      lump_sum=1;
      public_goods=0;
      benchmark_eqm=0;
      
      double taua0=taua;
      double taua_hi=strtof(argv[2],NULL);
      int n_tau = atoi(argv[3]);

      if(strcmp(argv[4],"1")==0)
	{
	  pe=1;
	  dist_max_iter=1;
	  wage_max_iter=1;
	  fprintf(logfile,"Short-run partial equilibrium\n");
	}

      double d = (taua_hi-taua0)/(n_tau-1);      
      for(int n=0; n<n_tau; n++)
	{
	  if(pe)
	    load_bin_eqm_t(ss1,fname0);
	  
	  taua=taua0+(n+1)*d;
	  
	  linebreak();
	  fprintf(logfile,"Solving for taua = %0.6f\n",taua);

	  if(pe==0)
	    {
	      snprintf(csv,128,"output/ss1_taua_%0.6f_2type_%0.3f.csv",taua,ss0->frac_v);
	      snprintf(pref1,128,"output/ss1_taua_%0.6f_2type_%0.3f",taua,ss0->frac_v);
	    }
	  else
	    {
	      snprintf(csv,128,"output/ss1_taua_%0.6f_2type_%0.3f_pe.csv",taua,ss0->frac_v);
	      snprintf(pref1,128,"output/ss1_taua_%0.6f_2type_%0.3f_pe",taua,ss0->frac_v);
	    }
	  solve_ss(ss1,pe);
	  if(pe)
	    calc_wealth_tax_evasion_elast(1,ss1);

	}
      
      free_mem();
    }

  else
    {
      help();
      return 1;
    }
  
  time(&stop);
  fprintf(logfile,"\nTotal time: %0.2f seconds\n",difftime(stop,start));

  linebreak();

  
  return 0;
}

#endif
