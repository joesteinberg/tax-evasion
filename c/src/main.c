#define EXTERN
#define EXTERN2

#ifndef __MAIN_C__
#define __MAIN_C__

#include "calibration.h"
#include "eqm.h"

void help()
{
  fprintf(logfile,"Bad user input!\n");

  fprintf(logfile,"Usage:\n");
  
  fprintf(logfile,"model 0 <a> <b>: steady state with current US policies\n");
  fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  fprintf(logfile,"\ta=1: Evasion (must compile with NV=12 and ND=2!!)\n");
  fprintf(logfile,"\tb=1: Higher transfer costs\n");
  fprintf(logfile,"\tb=2: Higher detection rate\n");
  fprintf(logfile,"\tb=3: Higher detection penalty\n");
  fprintf(logfile,"\tb=4: No detection revenenues");
  fprintf(logfile,"\tb=5: No collateral constraint\n");
  fprintf(logfile,"\tb=8: Hidden wealth cannot be collateralized\n");
  
  fprintf(logfile,"model 1 <a> <b> <c>: steady state with new tax system\n");
  fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  fprintf(logfile,"\tb=0: GKKOC\n");
  fprintf(logfile,"\tb=1: Warren\n");
  fprintf(logfile,"\tb=2: Sanders\n");
  fprintf(logfile,"\tb=3: Lower and upper tax cutoffs\n");
  fprintf(logfile,"\tc=1: Higher transfer costs\n");
  fprintf(logfile,"\tc=2: Higher detection rate\n");
  fprintf(logfile,"\tc=3: Higher detection penalty\n");
  fprintf(logfile,"\tc=4: No detection revenues\n");
  fprintf(logfile,"\tc=5: No collateral constraint\n");
  fprintf(logfile,"\tc=6: Short-run partial equilibrium\n");
  fprintf(logfile,"\tc=7: Small open economy\n");
  fprintf(logfile,"\tc=8: Hidden wealth cannot be collateralized\n");
 
  
  fprintf(logfile,"model 2 <a> <b>: transition to new tax system\n");
  fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  fprintf(logfile,"\ta=1: Evasion (must compile with NV=12 and ND=2!!)\n");
  fprintf(logfile,"\tb=1: Warren\n");
  fprintf(logfile,"\tb=2: Sanders\n");
  fprintf(logfile,"\tb=5: One-time COVID wealth tax\n");
  
  fprintf(logfile,"model 3 <a> <b> <c> <d> <e> <f>: series of capital income tax reforms\n");
  fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  fprintf(logfile,"\tb=1: Higher transfer costs\n");
  fprintf(logfile,"\tb=2: Higher detection rate\n");
  fprintf(logfile,"\tb=3: Higher detection penalty\n");
  fprintf(logfile,"\tb=4: No detection revenues\n");
  fprintf(logfile,"\tb=5: No collateral constraint\n");
  fprintf(logfile,"\tb=6: Short-run partial equilibrium\n");
  fprintf(logfile,"\tb=7: Small open economy\n");
  fprintf(logfile,"\tb=8: Hidden wealth cannot be collateralized\n");
  fprintf(logfile,"\tc: max tax cut\n");
  fprintf(logfile,"\td: max tax hike\n");
  fprintf(logfile,"\te: number of reforms in each direction\n");
  fprintf(logfile,"\tf: clear GBC using labor tax instead of lump-sum transfers\n");
  
  fprintf(logfile,"model 4 <a> <b> <c> <d> <e>: series of wealth tax reforms\n");
  fprintf(logfile,"\ta=0: No evasion (must compile with NV=1=ND=1!!)\n");
  fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  fprintf(logfile,"\tb=1: Higher transfer costs\n");
  fprintf(logfile,"\tb=2: Higher detection rate\n");
  fprintf(logfile,"\tb=3: Higher detection penalty\n");
  fprintf(logfile,"\tb=4: No detection revenues\n");
  fprintf(logfile,"\tb=5: No collateral constraint\n");
  fprintf(logfile,"\tb=6: Short-run partial equilibrium\n");
  fprintf(logfile,"\tb=7: Small open economy\n");
  fprintf(logfile,"\tb=8: Hidden wealth cannot be collateralized\n");
  fprintf(logfile,"\tb=16: Higher transfer costs + SR PE\n");
  fprintf(logfile,"\tb=26: Higher detection rate + SR PE\n");
  fprintf(logfile,"\tb=36: Higher detection penalty + SR PE\n");
  fprintf(logfile,"\tc: max wealth tax\n");
  fprintf(logfile,"\td: number of reforms\n");
  fprintf(logfile,"\te: clear GBC using labor tax instead of lump-sum transfers\n");
  
  fprintf(logfile,"model 6 <a> <b> <c>: Single progressive wealth tax\n");
  fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  fprintf(logfile,"\tb: Tax rate\n");
  fprintf(logfile,"\tc: Threshold (fraction of avg labor income)\n");
}

int main(int argc, char * argv[])
{
  time_t start, stop;
  time(&start);
  write_binary=1;
  logfile = stdout;
  linebreak();
  fprintf(logfile,"Tax Evasion and Capital Taxation\nShahar Rotberg and Joseph Steinberg\n\n");

#ifdef _OPENMP
  omp_set_num_threads(NTH);
  fprintf(logfile,"Parallel processing with %d threads\n",NTH);
#else
  fprintf(logfile,"No parallelization\n");
#endif
  
  alloc_mem();

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
  // fprintf(logfile,"model 0 <a> <b>: steady state with current US policies\n");
  // fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  // fprintf(logfile,"\ta=1: Evasion (must compile with NV=12 and ND=2!!)\n");
  // fprintf(logfile,"\tb=1: Higher transfer costs\n");
  // fprintf(logfile,"\tb=2: Higher detection rate\n");
  // fprintf(logfile,"\tb=3: Higher detection penalty\n");
  // fprintf(logfile,"\tb=4: No detection revenues\n");
  // fprintf(logfile,"\tb=5: Alternative constraint model\n");
  if(argc==4 && strcmp(argv[1],"0")==0)
    {
      fprintf(logfile,"Solving for initial steady state with current US tax code\n");

      if(strcmp(argv[3],"5")==0)
	{
	  constraint_flag=1;
	}

      if(NV==1)
	set_params(0);
      else
	set_params(1);     

      if(strcmp(argv[2],"0")==0 && NV==1 && ND==1)
	{
	  fprintf(logfile,"Scenario 0: no evasion\n");

	  if(strcmp(argv[3],"5")==0)
	    {
	      constraint_flag=1;
	      fprintf(logfile,"No collateral constraint\n");
	      snprintf(fname0,128,"output/ss0_no_evasion_noconst");
	      snprintf(csv,128,"output/ss0_no_evasion_noconst.csv");
	    }
	  else if(strcmp(argv[3],"0")==0)
	    {
	      snprintf(fname0,128,"output/ss0_no_evasion");
	      snprintf(csv,128,"output/ss0_no_evasion.csv");
	    }
	  else
	    {
	      help();
	      free_mem();
	      return 1;
	    }

	  init_eqm_t(ss0);
	  //load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");

	}     
      else if(strcmp(argv[2],"1")==0 && NV==12 && ND==2)
	{
	  fprintf(logfile,"Scenario 1: Evasion costs + detection penalty\n");

	  if(strcmp(argv[3],"1")==0)
	    {
	      eta=1.5*eta;
	      fprintf(logfile,"Higher transfer cost\n");
	      snprintf(fname0,128,"output/ss0_evasion_hi_ec1");
	      snprintf(csv,128,"output/ss0_evasion_hi_ec1.csv");
	    }
	  else if(strcmp(argv[3],"2")==0)
	    {
	      p2e = p2e * 2.0;
	      p2n = p2n * 2.0;
	      fprintf(logfile,"Higher detection rate\n");
	      snprintf(fname0,128,"output/ss0_evasion_hi_ec2");
	      snprintf(csv,128,"output/ss0_evasion_hi_ec2.csv");
	    }
	  else if(strcmp(argv[3],"3")==0)
	    {
	      penalty_frac_stock = penalty_frac_stock * 1.5;
	      penalty_frac_tauk = penalty_frac_tauk * 1.5;
	      penalty_frac_taua = penalty_frac_taua * 1.5;
	      fprintf(logfile,"Higher detection penalty\n");
	      snprintf(fname0,128,"output/ss0_evasion_hi_ec3");
	      snprintf(csv,128,"output/ss0_evasion_hi_ec3.csv");
	    }
	  else if(strcmp(argv[3],"4")==0)
	    {
	      detection_revenue_flag=1;
	      fprintf(logfile,"No detection revenues\n");
	      snprintf(fname0,128,"output/ss0_evasion_nodrev");
	      snprintf(csv,128,"output/ss0_evasion_nodrev.csv");
	    }
	  else if(strcmp(argv[3],"5")==0)
	    {
	      constraint_flag=1;
	      fprintf(logfile,"No collateral constraint\n");
	      snprintf(fname0,128,"output/ss0_evasion_noconstr");
	      snprintf(csv,128,"output/ss0_evasion_noconstr.csv");
	    }
	  else if(strcmp(argv[3],"8")==0)
	    {
	      chi=0.0;
	      fprintf(logfile,"Hidden wealth cannot be collateralized\n");
	      snprintf(fname0,128,"output/ss0_evasion_chi0");
	      snprintf(csv,128,"output/ss0_evasion_chi0.csv");
	    }

	  else if(strcmp(argv[3],"0")==0)
	    {
	      snprintf(fname0,128,"output/ss0_evasion");
	      snprintf(csv,128,"output/ss0_evasion.csv");
	    }
	  else
	    {
	      help();
	      free_mem();
	      return 1;
	    }

	  init_eqm_t(ss0);
	  //load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	  
	}
      else
	{
	  help();
	  free_mem();
	  return 1;
	}

      verbose=2;
      benchmark_eqm=1;

      /*
      if(strcmp(argv[3],"1")==0 || strcmp(argv[3],"2")==0 || strcmp(argv[3],"3")==0 || strcmp(argv[3],"4")==0)
	{
	benchmark_eqm=0;
	lump_sum=1;
	}
      */
      
      taua_clear_gbc=0;
      solve_ss(ss0,0);
      write_vf_dist(ss0,fname0);
      store_bin_eqm_t(ss0,fname0);
      
      if(strcmp(argv[2],"1")==0 && NV==12 && ND==2 && strcmp(argv[3],"0")==0)
	write_reg_data(ss0);
      
      free_mem();
    }
  //---------------------------------------------------------------------------
  // new steady state
  // fprintf(logfile,"model 1 <a> <b> <c>: steady state with new tax system\n");
  // fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  // fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  // fprintf(logfile,"\tb=0: GKKOC\n");
  // fprintf(logfile,"\tb=1: Warren\n");
  // fprintf(logfile,"\tb=2: Sanders\n");
  // fprintf(logfile,"\tb=3: Lower and upper tax cutoffs\n");
  // fprintf(logfile,"\tc=1: Higher transfer costs\n");
  // fprintf(logfile,"\tc=2: Higher detection rate\n");
  // fprintf(logfile,"\tc=3: Higher detection penalty\n");
  // fprintf(logfile,"\tc=4: No detection revenues\n");
  // fprintf(logfile,"\tc=5: Alternative constraint model\n");
  // fprintf(logfile,"\tc=6: Short-run partial equilibrium\n");
  // fprintf(logfile,"\tc=7: Small open economy\n");
  // fprintf(logfile,"\tc=8: Hidden wealth cannot be collateralized\n");
  else if(argc==5 && strcmp(argv[1],"1")==0)
    {
      fprintf(logfile,"Solving for steady state with new tax system\n");

      if(strcmp(argv[2],"0")==0 && NV==1 && ND==1)
	{
	  fprintf(logfile,"Scenario 0: no evasion\n");
	}
      else if(strcmp(argv[2],"1")==0 && NV==12 && ND==2)
	{
	  fprintf(logfile,"Scenario 1: Evasion costs + detection penalty\n");
	}
      else
	{
	  help();
	  free_mem();
	  return 1;
	}

      if(strcmp(argv[3],"5")==0)
	{
	  constraint_flag=1;
	}

      if(NV==1)
	set_params(0);
      else
	set_params(1);

      // parse the additional arguments
      int pe=0;
      char suff[24] = "";
      if(strcmp(argv[4],"6")==0)
	{
	  pe=1;
	  dist_max_iter=1;
	  wage_max_iter=1;
	  strcpy(suff,"_pe");
	  fprintf(logfile,"Short-run partial equilibrium\n");
	}
      else if(strcmp(argv[4],"7")==0)
	{
	  pe=2;
	  strcpy(suff,"_soe");
	  fprintf(logfile,"Small open economy\n");
	}

      else if(NV==12 && strcmp(argv[4],"1")==0)
	{
	  strcpy(suff,"_hi_ec1");
	  eta=1.5*eta;
	  fprintf(logfile,"Higher transfer costs\n");
	}
      else if(NV==12 && strcmp(argv[4],"2")==0)
	{
	  strcpy(suff,"_hi_ec2");
	  p2e = p2e*2.0;
	  p2n = p2n*2.0;
	  fprintf(logfile,"Higher detection rate\n");
	}
      else if(NV==12 && strcmp(argv[4],"3")==0)
	{
	  penalty_frac_stock = penalty_frac_stock*1.5;
	  penalty_frac_tauk = 1.5*penalty_frac_tauk;
	  penalty_frac_taua = 1.5*penalty_frac_taua;
	  strcpy(suff,"_hi_ec3");
	  fprintf(logfile,"Higher detection penalty\n");
	}
      else if(NV==12 && strcmp(argv[4],"4")==0)
	{
	  detection_revenue_flag=1;
	  strcpy(suff,"_nodrev");
	  fprintf(logfile,"No detection revenues\n");
	}
      else if(strcmp(argv[4],"5")==0)
	{
	  constraint_flag=1;
	  strcpy(suff,"_altconstr");
	  fprintf(logfile,"Alternative external financing setup\n");
	}
      else if(strcmp(argv[4],"8")==0)
	{
	  chi=0.0;
	  strcpy(suff,"_chi0");
	  fprintf(logfile,"Hidden wealth cannot be collateralized\n");
	}
      else if(strcmp(argv[4],"0")!=0)
	{
	  help();
	  free_mem();
	  return 1;
	}

      // set the file name to load initial eqm from      
      init_eqm_t(ss0);
      init_eqm_t(ss1);

      if(strcmp(argv[4],"6")==0 || strcmp(argv[4],"7")==0)
	{
	  if(NV==1)
	    snprintf(fname0,128,"output/ss0_no_evasion.bin");
	  else
	    snprintf(fname0,128,"output/ss0_evasion.bin");
	}
      else
	{
	  if(NV==1)
	    snprintf(fname0,128,"output/ss0_no_evasion%s.bin",suff);
	  else
	    snprintf(fname0,128,"output/ss0_evasion%s.bin",suff);
	}
      load_bin_eqm_t(ss0,fname0);
      load_bin_eqm_t(ss1,fname0);

      // set the output file names
      if(strcmp(argv[3],"0")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      snprintf(pref1,128,"output/ss1_gkkoc_no_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_gkkoc_no_evasion%s.csv",suff);
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      snprintf(pref1,128,"output/ss1_gkkoc_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_gkkoc_evasion%s.csv",suff);
	    }

	}
      else if(strcmp(argv[3],"1")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      snprintf(pref1,128,"output/ss1_warren_no_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_warren_no_evasion%s.csv",suff);
	      //load_bin_eqm_t(ss1,"output/ss1_warren_no_evasion.bin");
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      snprintf(pref1,128,"output/ss1_warren_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_warren_evasion%s.csv",suff);
	      //load_bin_eqm_t(ss1,"output/ss1_warren_evasion.bin");
	    }
	}
      else if(strcmp(argv[3],"2")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      snprintf(pref1,128,"output/ss1_sanders_no_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_sanders_no_evasion%s.csv",suff);
	      //load_bin_eqm_t(ss1,"output/ss1_sanders_no_evasion.bin");
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      snprintf(pref1,128,"output/ss1_sanders_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_sanders_evasion%s.csv",suff);
	      //load_bin_eqm_t(ss1,"output/ss1_sanders_evasion.bin");
	    }
	}
      else if(strcmp(argv[3],"3")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      snprintf(pref1,128,"output/ss1_prog_tauk_no_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_prog_tauk_no_evasion%s.csv",suff);
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      snprintf(pref1,128,"output/ss1_prog_tauk_evasion%s",suff);
	      snprintf(csv,128,"output/ss1_prog_tauk_evasion%s.csv",suff);
	    }
	}
     
      else
	{
	  help();
	  free_mem();
	  return 1;
	}  

      if(strcmp(argv[3],"0")==0)
	{
	  tauk=0.0;
	  taua_clear_gbc=1;
	  wealth_tax_type=0;
	  lump_sum=0;
	  public_goods=0;
	  abar=0.0;
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: Wealth tax instead of capital income tax\n");
	}
      else if(strcmp(argv[3],"1")==0)
	{
	  taua_clear_gbc=0;
	  wealth_tax_type=1;
	  lump_sum=1;
	  public_goods=0;
	  set_wealth_tax_params(ss0,wealth_tax_type);
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: Warren wealth tax\n");
	}
      else if(strcmp(argv[3],"2")==0)
	{
	  taua_clear_gbc=0;
	  wealth_tax_type=2;
	  lump_sum=1;
	  public_goods=0;
	  set_wealth_tax_params(ss0,wealth_tax_type);
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: Sanders wealth tax\n");
	}
      else if(strcmp(argv[3],"3")==0)
	{
	  taua_clear_gbc=0;
	  wealth_tax_type=0;
	  lump_sum=1;
	  public_goods=0;
	  set_wealth_tax_params(ss0,wealth_tax_type);
	  kbar = abar/2;
	  abar=0.0;
	  tauk2 = 0.01;
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: Progressive capital income tax (Biden plan)\n");
	}
      else
	{
	  help();
	  free_mem();
	  return 1;
	}  

      verbose=1;
      solve_ss(ss1,pe);
      if(strcmp(argv[4],"0")==0)
	{
	  write_vf_dist(ss1,pref1);
	  store_bin_eqm_t(ss1,pref1);
	}

      if(wealth_tax_type>=0 && abar>1.0e-8)
	calc_wealth_tax_evasion_elast(1,ss1);
      

      free_mem();
    }
  //---------------------------------------------------------------------------
  // transition
  // fprintf(logfile,"model 2 <a> <b>: transition to new tax system\n");
  // fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  // fprintf(logfile,"\ta=1: Evasion (must compile with NV=12 and ND=2!!)\n");
  // fprintf(logfile,"\tb=1: Warren\n");
  // fprintf(logfile,"\tb=2: Sanders\n");
  // fprintf(logfile,"\tb=5: One-time COVID wealth tax\n");
      
  else if(argc==4 && strcmp(argv[1],"2")==0)
    {
      fprintf(logfile,"Solving for transition to new tax system\n");

      if(strcmp(argv[2],"0")==0 && NV==1 && ND==1)
	{
	  fprintf(logfile,"Scenario 0: no evasion\n");
	  set_params(0);
	}
      else if(strcmp(argv[2],"1")==0 && NV==12 && ND==2)
	{
	  set_params(1);
	  fprintf(logfile,"Scenario 1: Evasion costs + detection penalty\n");
	}
      else
	{
	  help();
	  free_mem();
	  return 1;
	}

      alloc_trans_mem();

      init_eqm_t(ss0);
      init_eqm_t((eqm_trans[0]));
      init_eqm_t(ss1);
      init_eqm_t((eqm_trans[NT-1]));

      if(strcmp(argv[3],"1")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_no_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_warren_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_warren_no_evasion.bin");
	      snprintf(csv,128,"output/trans_warren_no_evasion.csv");
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_warren_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_warren_evasion.bin");
	      snprintf(csv,128,"output/trans_warren_evasion.csv");
	    }
	}
      else if(strcmp(argv[3],"2")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_no_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_sanders_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_sanders_no_evasion.bin");
	      snprintf(csv,128,"output/trans_sanders_no_evasion.csv");
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_sanders_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_sanders_evasion.bin");
	      snprintf(csv,128,"output/trans_sanders_evasion.csv");
	    }
	}
      else if(strcmp(argv[3],"5")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_no_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss0_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss0_no_evasion.bin");
	      snprintf(csv,128,"output/trans_covid_no_evasion.csv");
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss0_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss0_evasion.bin");
	      snprintf(csv,128,"output/trans_covid_evasion.csv");
	    }
	}
      else if(strcmp(argv[3],"k")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_no_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_tauk_0.344737_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_tauk_0.344737_no_evasion.bin");
	      snprintf(csv,128,"output/trans_tauk_no_evasion.csv");
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_tauk_0.344737_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_tauk_0.344737_evasion.bin");
	      snprintf(csv,128,"output/trans_tauk_evasion.csv");
	    }
	}
      else if(strcmp(argv[3],"a")==0)
	{	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_no_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_taua_0.021053_no_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_taua_0.021053_no_evasion.bin");
	      snprintf(csv,128,"output/trans_taua_no_evasion.csv");
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	      load_bin_eqm_t((eqm_trans[0]),"output/ss0_evasion.bin");
	      load_bin_eqm_t(ss1,"output/ss1_taua_0.020000_evasion.bin");
	      load_bin_eqm_t((eqm_trans[NT-1]),"output/ss1_taua_0.020000_evasion.bin");
	      snprintf(csv,128,"output/trans_taua_evasion.csv");
	    }
	}
      else
	{
	  help();
	  free_mem();
	  free_trans_mem();
	  return 1;
	}

         // copy initial steady state through the rest of the periods
      for(int t=1; t<NT-1; t++)
	{
	  copy_eqm_t((eqm_trans[t]),(eqm_trans[0]));
	}

      // copy over initial asset grid just to make sure they are all the same
      memcpy( (double *)(eqm_trans[NT-1]->a_grid),
	      (double *)(eqm_trans[0]->a_grid), NA*sizeof(double) );

      // construct initial guess
      double w0  = eqm_trans[0]->W;
      double w1  = eqm_trans[NT-1]->W;
      double wguess[NT];

      if(w0<w1)
	{
	  expspace(w0,w1,NT,2.5,wguess);
	}
      else
	{
	  expspace(w1,w0,NT,2.5,wguess);
	  reverse(wguess,NT);
	}

      double r0  = eqm_trans[0]->r;
      double r1  = eqm_trans[NT-1]->r;
      double rguess[NT];

      if(r0<r1)
	{
	  expspace(r0,r1,NT,2.5,rguess);
	}
      else
	{
	  expspace(r1,r0,NT,2.5,rguess);
	  reverse(rguess,NT);
	}
      	  
      double q0 = eqm_trans[0]->Q;
      double q1  = eqm_trans[NT-1]->Q;
      double qguess[NT];
      if(q0<q1)
	{
	  expspace(q0,q1,NT,2.5,qguess);
	}
      else
	{
	  expspace(q1,q0,NT,2.5,qguess);
	  reverse(qguess,NT);
	}

      double k0 = eqm_trans[0]->K;
      double k1  = eqm_trans[NT-1]->K;
      double kguess[NT];
      if(k0<k1)
	{
	  expspace(k0,k1,NT,2.5,kguess);
	}
      else
	{
	  expspace(k1,k0,NT,2.5,kguess);
	  reverse(kguess,NT);
	}
      
      for(int t=1; t<NT-1; t++)
	{
	  eqm_trans[t]->W = wguess[t];
	  eqm_trans[t]->r = rguess[t];
	  eqm_trans[t]->Q = qguess[t];
	  eqm_trans[t]->K = kguess[t];
	  eqm_trans[t]->P = alpha * pow(qguess[t],alpha-nu) *
	    pow(eqm_trans[t]->L,1.0-alpha-alpha2) * pow(kguess[t],alpha2);

	  if(gsl_isnan(eqm_trans[t]->P))
	    {
	      printf("%0.8f %0.8f %0.8f\n",qguess[t],kguess[t],eqm_trans[t]->L);
	      free_mem();
	      free_trans_mem();
	      return 1;
	    }
	}

      // set the tax params
      if(strcmp(argv[3],"1")==0)
	{
	  taua_clear_gbc=0;
	  wealth_tax_type=1;
	  lump_sum=1;
	  public_goods=0;
	  set_wealth_tax_params(ss0,wealth_tax_type);
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: Warren wealth tax\n");
	}
      else if(strcmp(argv[3],"2")==0)
	{
	  taua_clear_gbc=0;
	  wealth_tax_type=2;
	  lump_sum=1;
	  public_goods=0;
	  set_wealth_tax_params(ss0,wealth_tax_type);
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: Sanders wealth tax\n");
	}
      else if(strcmp(argv[3],"5")==0)
	{
	  taua_clear_gbc=0;
	  wealth_tax_type=5;
	  lump_sum=1;
	  public_goods=0;
	  set_wealth_tax_params(ss0,wealth_tax_type);
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: One-tiome COVID wealth tax\n");
	}
      else if(strcmp(argv[3],"k")==0)
	{
	  taua_clear_gbc=0;
	  lump_sum=1;
	  public_goods=0;
	  tauk = 0.344737;
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: 10pp increase in tauk\n");	  
	}
      else if(strcmp(argv[3],"a")==0)
	{
	  taua_clear_gbc=0;
	  lump_sum=1;
	  public_goods=0;
	  taua = 0.021053;
	  abar = 0.0;
	  benchmark_eqm=0;
	  fprintf(logfile,"Policy: 2pp increase in taua\n");	  
	}

      else
	{
	  help();
	  free_mem();
	  free_trans_mem();
	  return 1;
	}


      solve_trans();
      
      free_mem();
      free_trans_mem();
    }
      
  //---------------------------------------------------------------------------
  // fprintf(logfile,"model 3 <a> <b> <c> <d> <e> <f>: series of capital income tax reforms\n");
  // fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  // fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  // fprintf(logfile,"\tb=1: Higher transfer costs\n");
  // fprintf(logfile,"\tb=2: Higher detection rate\n");
  // fprintf(logfile,"\tb=3: Higher detection penalty\n");
  // fprintf(logfile,"\tb=4: No detection revenues\n");
  // fprintf(logfile,"\tb=5: Alternative constraint model\n");
  // fprintf(logfile,"\tb=6: Short-run partial equilibrium\n");
  // fprintf(logfile,"\tb=7: Small open economy\n");
  // fprintf(logfile,"\tb=8: Hidden wealth cannot be collateralized\n");
  // fprintf(logfile,"\tc: max tax cut\n");
  // fprintf(logfile,"\td: max tax hike\n");
  // fprintf(logfile,"\te: number of reforms in each direction\n");
  // fprintf(logfile,"\tf: clear GBC using labor tax instead of lump-sum transfers\n");
  else if(argc==8 && strcmp(argv[1],"3")==0)
    {
      fprintf(logfile,"Solving for series of steady states with different capital income tax rates\n");

      if(strcmp(argv[3],"5")==0)
	{
	  constraint_flag=1;
	}
      
      if(strcmp(argv[2],"0")==0 && NV==1 && ND==1)
	{
	  fprintf(logfile,"Scenario 0: no evasion\n");
	  set_params(0);
	}
      else if(strcmp(argv[2],"1")==0 && NV==12 && ND==2)
	{
	  set_params(1);
	  fprintf(logfile,"Scenario 1: Evasion costs + detection penalty\n");
	}
      else
	{
	  help();
	  free_mem();
	  return 1;
	}

      int pe=0;
      char suff[24] = "";
      
      if(strcmp(argv[3],"6")==0)
	{
	  pe=1;
	  dist_max_iter=0;
	  wage_max_iter=1;
	  strcpy(suff,"_pe");
	  fprintf(logfile,"Short-run partial equilibrium\n");
	}
      else if(strcmp(argv[3],"7")==0)
	{
	  pe=2;
	  strcpy(suff,"_soe");
	  fprintf(logfile,"Small open economy\n");
	}
      else if(NV==12 && strcmp(argv[3],"1")==0)
	{
	  strcpy(suff,"_hi_ec1");
	  eta = eta*1.5;
	  fprintf(logfile,"Higher transfer cost\n");
	}
      else if(NV==12 && strcmp(argv[3],"2")==0)
	{
	  strcpy(suff,"_hi_ec2");
	  p2e = p2e * 2.0;
	  p2n = p2n * 2.0;
	  fprintf(logfile,"Higher detection rate\n");
	}
      else if(NV==12 && strcmp(argv[3],"3")==0)
	{
	  strcpy(suff,"_hi_ec3");
	  penalty_frac_stock = penalty_frac_stock * 1.5;
	  penalty_frac_tauk = penalty_frac_tauk * 1.5;
	  penalty_frac_taua = penalty_frac_taua * 1.5;
	  fprintf(logfile,"Higher detection penalty\n");
	}
      else if(NV==12 && strcmp(argv[3],"4")==0)
	{
	  strcpy(suff,"_nodrev");
	  detection_revenue_flag=1;
	  fprintf(logfile,"No detection revenues\n");
	}
      else if(strcmp(argv[3],"5")==0)
	{
	  strcpy(suff,"_noconst");
	  constraint_flag=1;
	  fprintf(logfile,"No collateral constraint\n");
	}
      else if(strcmp(argv[3],"5")==0)
	{
	  strcpy(suff,"_chi0");
	  chi=0.0;
	  fprintf(logfile,"Hidden wealth cannot be collateralized\n");
	}


      else if(strcmp(argv[3],"0")!=0)
	{
	  help();
	  free_mem();
	  return 1;
	}

      // set the file to load initial eqm
      init_eqm_t(ss0);
      init_eqm_t(ss1);

      if(strcmp(argv[3],"6")==0 || strcmp(argv[3],"7")==0)
	{
	  if(NV==1)
	    snprintf(fname0,128,"output/ss0_no_evasion.bin");
	  else
	    snprintf(fname0,128,"output/ss0_evasion.bin");
	}
      else
	{
	  if(NV==1)
	    snprintf(fname0,128,"output/ss0_no_evasion%s.bin",suff);
	  else
	    snprintf(fname0,128,"output/ss0_evasion%s.bin",suff);
	}
      load_bin_eqm_t(ss0,fname0);
      load_bin_eqm_t(ss1,fname0);

      // set up the exercise
      taua_clear_gbc=0;
      lump_sum=1;
      public_goods=0;
      benchmark_eqm=0;
      verbose=1;

      double tauk0=tauk;
      double tauk_lo=tauk0 - strtof(argv[4],NULL);
      double tauk_hi=tauk0 + strtof(argv[5],NULL);
      int n_tau = atoi(argv[6]);

      taul_clear_gbc = atoi(argv[7]);
      if(taul_clear_gbc)
	{
	  fprintf(logfile,"Revenue used to reduce/increase labor income taxes, not lump-sum transfers\n");
	  lump_sum=0;
	  strcpy(suff,"_taul");
	}

      
      double d = (tauk0-tauk_lo)/(n_tau-1);

      for(int n=0; n<n_tau; n++)
	{
	  tauk=tauk0-(n+1)*d;

	  linebreak();
	  fprintf(logfile,"Solving for tauk = %0.6f\n",tauk);
	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      snprintf(csv,128,"output/ss1_tauk_%0.6f_no_evasion%s.csv",tauk,suff);
	      snprintf(pref1,128,"output/ss1_tauk_%0.6f_no_evasion%s",tauk,suff);
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      snprintf(csv,128,"output/ss1_tauk_%0.6f_evasion%s.csv",tauk,suff);
	      snprintf(pref1,128,"output/ss1_tauk_%0.6f_evasion%s",tauk,suff);
	    }

	  solve_ss(ss1,pe);
	  if(strcmp(argv[3],"0")==0)
	    {
	      store_bin_eqm_t(ss1,pref1);
	      write_vf_dist(ss1,pref1);
	    }
	}
      
      // now go upward
      //load_bin_eqm_t(ss1,fname0);
      d = (tauk_hi-tauk0)/(n_tau-1);
      
      for(int n=0; n<n_tau; n++)
	{
	  tauk=tauk0+(n+1)*d;
	  
	  linebreak();
	  fprintf(logfile,"Solving for tauk = %0.6f\n",tauk);
	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      snprintf(csv,128,"output/ss1_tauk_%0.6f_no_evasion%s.csv",tauk,suff);
	      snprintf(pref1,128,"output/ss1_tauk_%0.6f_no_evasion%s",tauk,suff);
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      snprintf(csv,128,"output/ss1_tauk_%0.6f_evasion%s.csv",tauk,suff);
	      snprintf(pref1,128,"output/ss1_tauk_%0.6f_evasion%s",tauk,suff);
	    }


	  solve_ss(ss1,pe);
	  if(strcmp(argv[3],"0")==0)
	    {
	      store_bin_eqm_t(ss1,pref1);
	      write_vf_dist(ss1,pref1);
	    }
	}
      
      free_mem();
    }
  //---------------------------------------------------------------------------
  // fprintf(logfile,"model 4 <a> <b> <c> <d> <e>: series of wealth tax reforms\n");
  // fprintf(logfile,"\ta=0: No evasion (must compile with NV=1=ND=1!!)\n");
  // fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  // fprintf(logfile,"\tb=1: Higher transfer costs\n");
  // fprintf(logfile,"\tb=2: Higher detection rate\n");
  // fprintf(logfile,"\tb=3: Higher detection penalty\n");
  // fprintf(logfile,"\tb=4: No detection revenues\n");
  // fprintf(logfile,"\tb=5: Alternative constraint model\n");
  // fprintf(logfile,"\tb=6: Short-run partial equilibrium\n");
  // fprintf(logfile,"\tb=7: Small open economy\n");
  // fprintf(logfile,"\tb=8: Hidden wealth cannot be collateralized"\n);
  // fprintf(logfile,"\tb=16: Higher transfer costs + SR PE\n");
  // fprintf(logfile,"\tb=26: Higher detection rate + SR PE\n");
  // fprintf(logfile,"\tb=36: Higher detection penalty + SR PE\n");
  // fprintf(logfile,"\tc: max wealth tax\n");
  // fprintf(logfile,"\td: number of reforms\n");
  // fprintf(logfile,"\te: clear GBC using labor tax instead of lump-sum transfers\n");

  else if(argc==7 && strcmp(argv[1],"4")==0)
    {
      fprintf(logfile,"Solving for series of steady states with different flat wealth tax rates\n");

      if(strcmp(argv[3],"5")==0)
	{
	  constraint_flag=1;
	}

      if(strcmp(argv[2],"0")==0 && NV==1 && ND==1)
	{
	  fprintf(logfile,"Scenario 0: no evasion\n");
	  set_params(0);
	}
      else if(strcmp(argv[2],"1")==0 && NV==12 && ND==2)
	{
	  set_params(1);
	  fprintf(logfile,"Scenario 1: Evasion costs + detection penalty\n");
	}
      else
	{
	  help();
	  free_mem();
	  return 1;
	}

      // parse the additional arguments
      int pe=0;
      char suff[24] = "";
      if(strcmp(argv[3],"6")==0)
	{
	  pe=1;
	  dist_max_iter=1;
	  wage_max_iter=1;
	  strcpy(suff,"_pe");
	  fprintf(logfile,"Short-run partial equilibrium\n");
	}
      else if(strcmp(argv[3],"7")==0)
	{
	  pe=2;
	  strcpy(suff,"_soe");
	  fprintf(logfile,"Small open economy\n");
	}

      else if(NV==12 && strcmp(argv[3],"1")==0)
	{
	  strcpy(suff,"_hi_ec1");
	  eta = eta*1.5;
	  fprintf(logfile,"Higher transfer cost\n");
	}
      else if(NV==12 && strcmp(argv[3],"2")==0)
	{
	  strcpy(suff,"_hi_ec2");
	  p2e = p2e * 2.0;
	  p2n = p2n * 2.0;
	  fprintf(logfile,"Higher detection rate\n");
	}
      else if(NV==12 && strcmp(argv[3],"3")==0)
	{
	  strcpy(suff,"_hi_ec3");
	  penalty_frac_stock = penalty_frac_stock * 1.5;
	  penalty_frac_tauk = penalty_frac_tauk * 1.5;
	  penalty_frac_taua = penalty_frac_taua * 1.5;
	  fprintf(logfile,"Higher detection penalty\n");
	}
      else if(NV==12 && strcmp(argv[3],"4")==0)
	{
	  strcpy(suff,"_nodrev");
	  detection_revenue_flag=1;
	  fprintf(logfile,"No detection revenues\n");
	}
      else if(strcmp(argv[3],"5")==0)
	{
	  strcpy(suff,"_noconst");
	  constraint_flag=1;
	  fprintf(logfile,"No collateral constraint\n");
	}
      else if(strcmp(argv[3],"8")==0)
	{
	  strcpy(suff,"_chi0");
	  fprintf(logfile,"Hidden wealth cannot be collateralized\n");
	}
      else if(NV==12 && strcmp(argv[3],"16")==0)
	{
	  pe=1;
	  dist_max_iter=1;
	  wage_max_iter=1;
	  strcpy(suff,"_hi_ec1");
	  eta = eta*1.5;
	  fprintf(logfile,"Higher transfer cost + SR PE\n");
	}
      else if(NV==12 && strcmp(argv[3],"26")==0)
	{
	  pe=1;
	  dist_max_iter=1;
	  wage_max_iter=1;
	  strcpy(suff,"_hi_ec2");
	  p2e = p2e * 2.0;
	  p2n = p2n * 2.0;
	  fprintf(logfile,"Higher detection rate + SR PE\n");
	}
      else if(NV==12 && strcmp(argv[3],"36")==0)
	{
	  pe=1;
	  dist_max_iter=1;
	  wage_max_iter=1;
	  strcpy(suff,"_hi_ec3");
	  penalty_frac_stock = penalty_frac_stock * 1.5;
	  penalty_frac_tauk = penalty_frac_tauk * 1.5;
	  penalty_frac_taua = penalty_frac_taua * 1.5;
	  fprintf(logfile,"Higher detection penalty + SR PE\n");
	}


      else if(strcmp(argv[3],"0")!=0)
	{
	  help();
	  free_mem();
	  return 1;
	}

      // set the file to load initial eqm
      init_eqm_t(ss0);
      init_eqm_t(ss1);

      if(strcmp(argv[3],"6")==0 || strcmp(argv[3],"7")==0)
	{
	  if(NV==1)
	    snprintf(fname0,128,"output/ss0_no_evasion.bin");
	  else
	    snprintf(fname0,128,"output/ss0_evasion.bin");
	}
      else
	{
	  if(NV==1)
	    snprintf(fname0,128,"output/ss0_no_evasion%s.bin",suff);
	  else
	    snprintf(fname0,128,"output/ss0_evasion%s.bin",suff);
	}
      load_bin_eqm_t(ss0,fname0);
      load_bin_eqm_t(ss1,fname0);

      if(strcmp(argv[3],"16")==0)
	{
	  strcpy(suff,"_hi_ec1_pe");
	}
      else if(strcmp(argv[3],"26")==0)
	{
	  strcpy(suff,"_hi_ec2_pe");
	}
      else if(strcmp(argv[3],"36")==0)
	{
	  strcpy(suff,"_hi_ec3_pe");
	}

      // set up the exercise
      verbose=0;
      taua_clear_gbc=0;
      taul_clear_gbc=0;
      lump_sum=1;
      public_goods=0;
      benchmark_eqm=0;
      
      double taua0=taua;
      double taua_hi=strtof(argv[4],NULL);
      int n_tau = atoi(argv[5]);

      taul_clear_gbc = atoi(argv[6]);
      if(taul_clear_gbc)
	{
	  fprintf(logfile,"Revenue used to reduce/increase labor income taxes, not lump-sum transfers\n");
	  lump_sum=0;
	  strcpy(suff,"_taul");
	}


      double d = (taua_hi-taua0)/(n_tau-1);      
      for(int n=0; n<n_tau; n++)
	{
	  if(pe)
	    load_bin_eqm_t(ss1,fname0);
	    
	  taua=taua0+(n+1)*d;
	  
	  linebreak();
	  fprintf(logfile,"Solving for taua = %0.6f\n",taua);
	  
	  if(strcmp(argv[2],"0")==0)
	    {
	      snprintf(csv,128,"output/ss1_taua_%0.6f_no_evasion%s.csv",taua,suff);
	      snprintf(pref1,128,"output/ss1_taua_%0.6f_no_evasion%s",taua,suff);
	    }
	  else if(strcmp(argv[2],"1")==0)
	    {
	      snprintf(csv,128,"output/ss1_taua_%0.6f_evasion%s.csv",taua,suff);
	      snprintf(pref1,128,"output/ss1_taua_%0.6f_evasion%s",taua,suff);
	    }

	  solve_ss(ss1,pe);
	  if(strcmp(argv[3],"0")==0)
	    {
	      write_vf_dist(ss1,pref1);
	      store_bin_eqm_t(ss1,pref1);
	    }

	}
      
      free_mem();
    }
  //---------------------------------------------------------------------------
  //fprintf(logfile,"model 6 <a> <b> <c>: Single progressive wealth tax\n");
  //fprintf(logfile,"\ta=0: No evasion (must compile with NV=ND=1!!)\n");
  //fprintf(logfile,"\ta=1: Evasion costs (must compile with NV=12 and ND=2!!)\n");
  //fprintf(logfile,"\tb: Tax rate\n");
  //fprintf(logfile,"\tc: Threshold (fraction of avg labor income)\n");
  else if(argc==5 && strcmp(argv[1],"6")==0)
    {
      double taua_=strtof(argv[3],NULL);
      double abar_=strtof(argv[4],NULL);

      fprintf(logfile,"Solving for single prog tax with wealth tax rate of %0.6f and threshold = %0.6f\n",taua_,abar_);
      
      if(strcmp(argv[2],"0")==0 && NV==1)
	{
	  fprintf(logfile,"Scenario 0: no evasion\n");
	  set_params(0);
	}
      else if(strcmp(argv[2],"1")==0 && NV==12)
	{
	  set_params(1);
	  fprintf(logfile,"Scenario 1: Evasion costs + detection penalty\n");
	}
      else
	{
	  help();
	  free_mem();
	  return 1;
	}

      // load the initial eqm and set the output file name
      init_eqm_t(ss0);
      init_eqm_t(ss1);
      if(strcmp(argv[2],"0")==0)
	{
	  load_bin_eqm_t(ss0,"output/ss0_no_evasion.bin");
	  load_bin_eqm_t(ss1,"output/ss0_no_evasion.bin");
	  snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_no_evasion.csv",taua_,abar_);
	  snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_no_evasion.bin",taua_,abar_);
	}
      else
	{
	  load_bin_eqm_t(ss0,"output/ss0_evasion.bin");
	  load_bin_eqm_t(ss1,"output/ss0_evasion.bin");
	  snprintf(csv,128,"output/ss_opt_%0.6f_%0.6f_evasion.csv",taua_,abar_);
	  snprintf(pref1,128,"output/ss_opt_%0.6f_%0.6f_evasion.bin",taua_,abar_);
	}
            
      // set up the exercise
      taua_clear_gbc=0;
      lump_sum=1;
      public_goods=0;
      benchmark_eqm=0;
      wealth_tax_type=0;     
      abar=abar_/ss0->ylbar;
      taua=taua_;
      verbose=2;
      solve_ss(ss1,0);

      write_vf_dist(ss1,pref1);
      
      if(wealth_tax_type>=0 && abar>1.0e-8)
	calc_wealth_tax_evasion_elast(1,ss1);

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
