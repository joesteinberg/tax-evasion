#ifndef __EXTERNS_C__
#define __EXTERNS_C__

#include "calibration.h"
#include "eqm.h"


extern void setup()
{
  alloc_mem();

  constraint_flag=0;
  detection_revenue_flag=0;
  benchmark_eqm = 0;
  taua_clear_gbc=0;
  taul_clear_gbc=0;
  wealth_tax_type=0;
  lump_sum=1;
  public_goods=0;

  set_params(evasion_type);
  init_eqm_t(ss0);
  init_eqm_t(ss1);

  if(evasion_type==0)
    {
      load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");
      load_bin_eqm_t(ss1,"output/ss0_no_evasion.bin");
    }
  else
    {
      if(sens_type==0)
	{
	  load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	  load_bin_eqm_t(ss1,"output/ss0_evasion.bin");
	}
      else if(sens_type==4)
	{
	  load_bin_eqm_t(ss0,"output/ss0_evasion_nodrev.bin");
	  load_bin_eqm_t(ss1,"output/ss0_evasion_nodrev.bin");
	  detection_revenue_flag=1;
	}
      else if(sens_type==1)
	{
	  load_bin_eqm_t(ss0,"output/ss0_evasion_hi_ec1.bin");
	  load_bin_eqm_t(ss1,"output/ss0_evasion_hi_ec1.bin");
	  eta=1.5*eta;
	}
      else if(sens_type==2)
	{
	  load_bin_eqm_t(ss0,"output/ss0_evasion_hi_ec2.bin");
	  load_bin_eqm_t(ss1,"output/ss0_evasion_hi_ec2.bin");
	  p2e = p2e * 2.0;
	  p2n = p2n * 2.0;
	      
	}
      else if(sens_type==3)
	{
	  load_bin_eqm_t(ss0,"output/ss0_evasion_hi_ec3.bin");
	  load_bin_eqm_t(ss1,"output/ss0_evasion_hi_ec3.bin");
	  penalty_frac_stock = penalty_frac_stock * 1.5;
	  penalty_frac_tauk = penalty_frac_tauk * 1.5;
	  penalty_frac_taua = penalty_frac_taua * 1.5;

	}

    }

}

extern void cleanup()
{
  free_mem();
}

/*
objective function for calibration program 
1. double taua: tax rate
2. double abar: threshold where wealth tax starts
*/

extern double cal_obj_fun(int D, double targeted_params[])
{
  if(D!=2)
    {
      return GSL_NAN;
    }
  else
    {
      double taua_ = targeted_params[0];
      double abar_ = targeted_params[1];

      fprintf(logfile,"\n----------------------------------------------------------------------------------------------\n");
      fprintf(logfile,"----------------------------------------------------------------------------------------------\n");
      fprintf(logfile,"----------------------------------------------------------------------------------------------\n");
      fprintf(logfile,"\nSolving for equilibrium given (taua, abar) = (%0.6f,%0.6f)\n\n",taua_,abar_);
      fflush(logfile);

      /*
      if(evasion_type==0)
	{
	  load_bin_eqm_t(ss1,"output/ss0_no_evasion.bin");
	}
      else
	{
	  if(sens_type==0)
	    load_bin_eqm_t(ss1,"output/ss0_evasion.bin");
	  else if(sens_type==4)
	    load_bin_eqm_t(ss1,"output/ss0_evasion_nodrev.bin");
	  else if(sens_type==1)
	    load_bin_eqm_t(ss1,"output/ss0_evasion_hi_ec1.bin");
	  else if(sens_type==2)
	    load_bin_eqm_t(ss1,"output/ss0_evasion_hi_ec2.bin");
	  else if(sens_type==3)
	    load_bin_eqm_t(ss1,"output/ss0_evasion_hi_ec3.bin");
	}
      */
      
      taua=taua_;
      abar = abar_/ss0->ylbar;

      if(evasion_type==0)
	{
	  snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_no_evasion",taua_,abar_);
	  snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_no_evasion.csv",taua_,abar_);
	}
      else
	{
	  if(sens_type==0)
	    {
	      snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_evasion",taua_,abar_);
	      snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_evasion.csv",taua_,abar_);
	    }
	  else if(sens_type==1)
	    {
	      snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_evasion_hi_ec1",taua_,abar_);
	      snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_evasion_hi_ec1.csv",taua_,abar_);
	    }
	  else if(sens_type==2)
	    {
	      snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_evasion_hi_ec2",taua_,abar_);
	      snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_evasion_hi_ec2.csv",taua_,abar_);
	    }
	  else if(sens_type==3)
	    {
	      snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_evasion_hi_ec3",taua_,abar_);
	      snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_evasion_hi_ec3.csv",taua_,abar_);
	    }
	  else if(sens_type==4)
	    {
	      snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_evasion_nodrev",taua_,abar_);
	      snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_evasion_nodrev.csv",taua_,abar_);
	    }

	}

      if(solve_ss(ss1,0))
	{
	  fprintf(logfile,"\nFailed to solve for equilibrium! Sending NaN to master thread...\n");
	  fflush(logfile);
	  return +99999;
	}
  
      fprintf(logfile,"Welfare: %0.16f\n", ss1->welfare_newborn);
      fprintf(logfile,"Storing equilibrium data...\n");
      
      //store_bin_eqm_t(ss1,pref1);
      //write_vf_dist(ss1,pref1);
      write_eqm_t_csv(ss1,csv);

      fprintf(logfile,"\n----------------------------------------------------------------------------------------------\n");
      fprintf(logfile,"----------------------------------------------------------------------------------------------\n");
      fprintf(logfile,"----------------------------------------------------------------------------------------------\n");
      
      fflush(logfile);

      return -ss1->welfare_newborn;
    }
}

#endif
