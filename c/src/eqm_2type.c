#ifndef __EQM_2TYPE_C__
#define __EQM_2TYPE_C__

#include "eqm_2type.h"

/****************************************************************/
/*    Functions related to memory allocation                    */
/****************************************************************/

void alloc_mem()
{
  fprintf(logfile,"Allocating memory\n");

  ss0 = (eqm_t *)malloc(sizeof(eqm_t));
  ss1 = (eqm_t *)malloc(sizeof(eqm_t));
  
  int tn, j, ie, iz, ii, iv, id;
  for(tn=0; tn<NTH; tn++)
    {
      acc[tn] = gsl_interp_accel_alloc();
      spline_V[tn] = (gsl_spline *******)malloc(sizeof(gsl_spline ******)*J);
      
      for(j=0; j<J; j++)
	{
	  spline_V[tn][j] = (gsl_spline ******)malloc(sizeof(gsl_spline *****)*NE);
	  
	  for(ie=0; ie<NE; ie++)
	   {
	     spline_V[tn][j][ie] = (gsl_spline *****)malloc(sizeof(gsl_spline ****)*NZ);
	     
	      for(iz=0; iz<NZ; iz++)
		{
		  spline_V[tn][j][ie][iz] = (gsl_spline ****)malloc(sizeof(gsl_spline ***)*NI);

		  for(ii=0; ii<NI; ii++)
		    {
		      spline_V[tn][j][ie][iz][ii] = (gsl_spline ***)malloc(sizeof(gsl_spline **)*NV);
		      
		      for(iv=0; iv<NV; iv++)
			{
			  spline_V[tn][j][ie][iz][ii][iv] = (gsl_spline **)malloc(sizeof(gsl_spline *)*ND);
			  
			  for(id=0; id<ND; id++)
			    {
			      spline_V[tn][j][ie][iz][ii][iv][id] = gsl_spline_alloc(gsl_interp_linear,NA);
			    }
			}
		    }
		}
	    }
	}
    }
}

void free_mem()
{
  fprintf(logfile,"Freeing memory\n");

  free(ss0);
  free(ss1);
    
  int tn, j, ie, iz, ii, iv, id;
  for(tn=0; tn<NTH; tn++)
    {
      gsl_interp_accel_free(acc[tn]);

      for(j=0; j<J; j++)
	{
	  for(ie=0; ie<NE; ie++)
	    {
	      for(iz=0; iz<NZ; iz++)
		{
		  for(ii=0; ii<NI; ii++)
		    {
		      for(iv=0; iv<NV; iv++)
			{
			  for(id=0; id<ND; id++)
			    {
			      gsl_spline_free(spline_V[tn][j][ie][iz][ii][iv][id]);
			    }

			  free(spline_V[tn][j][ie][iz][ii][iv]);
			}
		      
		      free(spline_V[tn][j][ie][iz][ii]);
		    }

		  free(spline_V[tn][j][ie][iz]);
		}

	      free(spline_V[tn][j][ie]);
	    }
	  
	  free(spline_V[tn][j]);
	}
      
      free(spline_V[tn]);
    }
}

/****************************************************************/
/*    Functions related to initializing the equilibrium         */
/****************************************************************/

void copy_eqm_t(eqm_t * dest, const eqm_t * src)
{
  dest->L = src->L;
  dest->Q = src->Q;
  dest->K = src->K;
  dest->W = src->W;
  dest->r = src->r;
  dest->P = src->P;
  dest->ylbar = src->ylbar;
  dest->G = src->G;
  dest->lump_sum=src->lump_sum;
  
  memcpy( (double *)(dest->a_grid), (double *)(src->a_grid), NA*sizeof(double) );
  
  memcpy( (double *)(dest->ga), (double *)(src->ga), 2*J*NE*NZ*NI*NA*1*ND*sizeof(double) );
  memcpy( (int *)(dest->gv), (int *)(src->gv), 2*J*NE*NZ*NI*NA*1*ND*sizeof(int) );
  memcpy( (double *)(dest->V), (double *)(src->V), 2*J*NE*NZ*NI*NA*1*ND*sizeof(double) );
  memcpy( (double *)(dest->Psi), (double *)(src->Psi), 2*J*NE*NZ*NI*NA*1*ND*sizeof(double) );

}

void initialize_asset_grid(double W, double * a_grid)
{
  expspace(0.0,asset_grid_ub_mult*e_grid[NE-1]*W,NA,asset_grid_exp,a_grid);
}

void init_eqm_t(eqm_t * et)
{
  int j, ie, iz, ii, ia, iv, id;

  et->Q = 333.60737688359467;
  et->K = 65.750225513212570;
  et->L = 30.155518045063850;
  et->W = 1.4183938503265381;
  et->r = 2.9629606753587723E-002;
  et->lump_sum=0.0;
  et->ylbar = 1.0895181505953169;
  
  et->P = alpha * pow(et->Q,alpha-nu) * pow(et->L,1.0-alpha-alpha2) * pow(et->K,alpha2);
  et->Y = pow(et->Q,alpha)*pow(et->K,alpha2)*pow(et->L,1.0-alpha-alpha2);
  et->G = 0.0;
  et->welfare = -HUGE_VAL;
  et->welfare_newborn = -HUGE_VAL;

  if(benchmark_eqm)
    {
      initialize_asset_grid(et->W,et->a_grid);
    }

  SET_ALL_V(et->V,2*J*NE*NZ*NI*NA*NV*ND,0.0);
  
  init_dist(et);
  set_capital_policies(et);

  for(j=0; j<J; j++)
    {
      for(ie=0; ie<NE; ie++)
	{
	  for(iz=0; iz<NZ; iz++)
	    {
	      for(ii=0; ii<NI; ii++)
		{
		  for(iv=0; iv<NV; iv++)
		    {
		      for(id=0; id<ND; id++)
			{
			  for(ia=0; ia<NA; ia++)
			    {
			      et->ga[1][j][ie][iz][ii][iv][id][ia] = beta*et->a_grid[ia];
			      et->gv[1][j][ie][iz][ii][iv][id][ia] = 0;

			      if(iv==0 && id==0)
				et->ga[0][j][ie][iz][ii][iv][id][ia] = beta*et->a_grid[ia];
			      
			    }
			}
		    }
		}
	    }
	}
    }
}

void init_splines(eqm_t * et, int j, int v)
{
  int ie;
  for(ie=0; ie<NE; ie++)
    {
      int iz;
      for(iz=0; iz<NZ; iz++)
	{
	  int ii;
	  for(ii=0; ii<NI; ii++)
	    {
	      int iv;
	      int iv_ub = (evasion_type==1 && v==1 ? NV : 1);
	      for(iv=0; iv<iv_ub; iv++)
		{
		  int id;
		  int id_ub = (evasion_type==1 && iv>0 ? ND: 1);
		  for(id=0; id<id_ub; id++)
		    {
		      for(int tn=0; tn<NTH; tn++)
			{
			  gsl_spline_init(spline_V[tn][j][ie][iz][ii][iv][id],
					  et->a_grid,
					  et->V[v][j][ie][iz][ii][iv][id],
					  NA);
			  
			}
		      
		    }
		}
	    }
	}
    }
}

/****************************************************************/
/*    Functions related to constructing the income policy       */
/****************************************************************/

void set_wealth_tax_params(const eqm_t * et, int wealth_tax_type)
{
  if(wealth_tax_type==1)
    {
      taua=0.02;
    }
  else if(wealth_tax_type==2)
    {
      taua=0.01;
    }
  else if(wealth_tax_type==3)
    {
      taua=0.02;
    }
  else if(wealth_tax_type==5)
    {
      taua=0.02;
    }
  else if(wealth_tax_type==6)
    {
      taua=0.015;
    }
  
  double total_mass = SUM(et->Psi,2*J*NE*NZ*NI*NV*NA*ND);
  
  double mass=0.0;
  
  for(int ia=NA-1; ia>=0; ia--)
    {
      for(int iv=NV-1; iv>=0; iv--)
	{
	  for(int id=0; id<ND; id++)
	    {
	      for(int j=0; j<J; j++)
		{
		  for(int ie=0; ie<NE; ie++)
		    {
		      for(int iz=0; iz<NZ; iz++)
			{
			  for(int ii=0; ii<NI; ii++)
			    {
			      mass += et->Psi[0][j][ie][iz][ii][iv][id][ia];
			      mass += et->Psi[1][j][ie][iz][ii][iv][id][ia];
			      
			      if((wealth_tax_type <= 2 || wealth_tax_type==6) && mass/total_mass>=0.000583)
				{
				  abar=et->a_grid[ia];
				  printf("%0.8f\n",abar);
				  return;
				}
			      else if(wealth_tax_type == 3)
				{
				  if(abar2<1.0e-8 && mass/total_mass >= 0.001)
				    {
				      abar2=et->a_grid[ia];
				    }
				  else if(abar<1.0e-8 && mass/total_mass >= 0.01)
				    {
				      abar=et->a_grid[ia];
				      return;
				    }
				}
			      else if(wealth_tax_type == 5 && mass/total_mass>=0.000583)
				{
				  abar=et->a_grid[ia];
				  printf("%0.8f\n",abar);
				  return;
				}
			    }
			}
		    }
		}
	    }
	}
    }
}



void set_capital_policies(eqm_t * et)
{
  int iz;
  for(iz=0; iz<NZ; iz++)
    {
      int ii;
      for(ii=0; ii<NI; ii++)
	{
	  double z=z_grid[iz];
	  if(ii==1)
	    {
	      z = 0.0;
	    }

	  int iv;
	  for(iv=0; iv<NV; iv++)
	    {
	      int ia;
	      for(ia=0; ia<NA; ia++)
		{
		  double a = et->a_grid[ia];
		  et->gks[iz][ii][iv][ia] = a;
		  double at = a*(1-evasion_grid[iv]);
		  double av = a*evasion_grid[iv];
		  
		  if(iz==0 && ii==0)
		    {
		      if(wealth_tax_type==0 || wealth_tax_type==6)
			{
			  et->gta[iv][ia] = taua*fmax(0.0,at-abar);
			  et->gta_evaded[iv][ia] = taua*fmax(0.0,a-abar) - et->gta[iv][ia];
			}
		      else if(wealth_tax_type==1 || wealth_tax_type==5)
			{
			  et->gta[iv][ia] = taua * fmax(0.0,at-abar) + (taua>1.0e-10)*(taua+0.01)*fmax(0.0,at-abar*20.);
			  et->gta_evaded[iv][ia] = taua*fmax(0.0,a-abar) + (taua>1.0e-10)*(taua+0.01)*fmax(0.0,a-abar*20.)
			    - et->gta[iv][ia];	      
			}
		      else if(wealth_tax_type==2)
			{
			  et->gta[iv][ia] = 0.0;
			  double gta_owed=0.0;
			  double abar_[8] = {abar*32./50.,abar,abar*5.,abar*10.,abar*20.,abar*50.,abar*100.,abar*200.};
			  double taua_[8] = {taua,taua+0.01,taua+0.02,taua+0.03,taua+0.04,taua+0.05,taua+0.06,taua+0.07};
			  for(int ibar=0; ibar<8; ibar++)
			    {
			      if(at>abar_[ibar])
				{
				  double tmpt=0.0;
				  if(ibar<8-1)
				    {
				      tmpt = (fmin(at,abar_[ibar+1])-abar_[ibar]);
				    }
				  else
				    {
				      tmpt = (at-abar_[ibar]);
				    }
				  et->gta[iv][ia] += taua_[ibar] * tmpt;
				}

			      if(a>abar_[ibar])
				{
				  double tmp=0.0;
				  if(ibar<8-1)
				    {
				      tmp = (fmin(a,abar_[ibar+1])-abar_[ibar]);
				    }
				  else
				    {
				      tmp = (a-abar_[ibar]);
				    }
				  gta_owed += taua_[ibar] * tmp;
				}

			      if(et->gta[iv][ia]<0.0 && ibar<(8-1))
				{
				  printf("Error!\n");
				  printf("%0.8f,%0.8f,%0.8f,%0.8f\n",a,at,abar_[ibar],abar_[ibar+1]);
				}
			    }
			  et->gta_evaded[iv][ia] = gta_owed - et->gta[iv][ia];

			}
		      else if(wealth_tax_type==3)
			{
			  et->gta[iv][ia] = taua * fmax(0.0,fmin(abar2,at)-abar);
			  et->gta_evaded[iv][ia] = taua*fmax(0.0,fmin(a,abar2)-abar) - et->gta[iv][ia];	      
			}
		    }
		  
		  // if you got unlucky and lost your idea, you get no capital income, just passive interest income
		  if(ii==1)
		    {
		      et->gk[iz][ii][iv][ia] = 0.0;
		      
		      et->gtk[iz][ii][iv][ia] = tauk*at*fmax(et->r,0.0);
		      et->gtk[iz][ii][iv][ia] += tauk2* fmax(0.0, at*fmax(et->r,0.0) - kbar);
		      
		      et->gyk[iz][ii][iv][ia] = (a*et->r) - et->gtk[iz][ii][iv][ia];
		    }
		  // otherwise you operate your business
		  else
		    {
		      // levered profit
		      double kstar = pow( (nu*et->P*pow(z,nu))/(et->r+delta), 1.0/(1.0-nu) );
		      double kconst=0.0;
		      
		      if(constraint_flag==1)  
			{
			  // k = a_r + d, d<=max{lambda*(a_r+a_h),lambda2*pi}
			  // k - a_r <= max{lambda*(a_r+a_h),lambda2*pi}
			  // k <= max{lambda*(a_r+a_h),lambda2*pi} + a_r

			  kconst = fmax(lambda*(at + chi*av) , kstar);
			  double kpi=et->P*pow(z*kconst,nu);
			  //kconst = fmax(lambda*a,lambda2*pi) + at;

			  double diff=HUGE_VAL;
			  int iter=0;
			  do
			    {
			      iter++;
			      kpi = et->P*pow(z*kconst,nu);
			      double tmp = fmax(lambda*(at + chi*av),lambda2*kpi) + at;
			      diff = fabs(tmp-kconst);
			      kconst = tmp;
			    }
			  while(diff>1.0e-10 && iter<10000);
			  if(iter==10000)
			    {
			      fprintf(logfile,"Failed to converge constraint!\n");
			    }

			  if(fabs(kconst - (lambda*(at+chi*av)+at)) > 1e-10)
			    {
			      if(fabs(kconst - (lambda2*kpi + at)) > 1.0e-10)
				{
				  fprintf(logfile,"Error in constraint calculation!\n");
				  fprintf(logfile,"%0.6f    %0.6f    %0.6f\n",kconst,lambda2*kpi + at, lambda*(at+chi*av)+at);
				}
			    }
			}
		      else
			{
			  double lambda_ = lambda*((double)(iz)/((double)(NZ)-1.0));			      
			  kconst = (1.0+lambda_)*at + lambda_*chi*av;
			}
		      
		      double k=kstar;
		      if(kstar>kconst)
			{
			  k=kconst;
			}
		      et->gk[iz][ii][iv][ia] = k;
		      et->gks[iz][ii][iv][ia] = a;

		      double kpi = et->P*pow(z*k,nu) - (et->r+delta)*k + et->r*at;
		      
		      et->gtk[iz][ii][iv][ia] = tauk*kpi;
		      et->gtk[iz][ii][iv][ia] += tauk2*fmax(0.0,kpi-kbar);
		      et->gyk[iz][ii][iv][ia] = kpi + et->r*av - et->gtk[iz][ii][iv][ia];

		  
		      if(gsl_isnan(et->gyk[iz][ii][iv][ia]))
			{
			  fprintf(logfile,"NAN Capital income!\n");
			  fprintf(logfile,"%0.8f\n",k);
			  fprintf(logfile,"%0.8f\n",at);
			  fprintf(logfile,"%0.8f\n",et->gtk[iz][ii][iv][ia]);
			  fprintf(logfile,"r=%0.8f\n",et->r);
			  fprintf(logfile,"P=%0.8f\n",et->P);
			  fprintf(logfile,"ksyar=%0.8f\n",kstar);
			  fprintf(logfile,"test!\n");
			}			 
		    }
		}
	    }
	}
    }
}

double calc_ktax_evaded(double ktax, int iv, int ivp, double a, double ap)
{
  if(tauk<0)
    return 0.0;
  else
    return fmin(ktax,tauk*fmax(0.0,ap*evasion_grid[ivp] - a*evasion_grid[iv]));
}

void calc_c_l(int j, int ie, double W, double inc_nol_post_saving_evasion, double * c, double * l)
{
    if(j<R)
    {
      // (1) c(1+tauc)+ap = (1-taul)*W*l*zeta_j*e + yk - gta
      // (2) (1+tauc)/((1-taul)*W*zeta_j*e) = (gama/(1-gama))*(1-l)/c
      // (3) (2)--> c*(1+tauc) =  (1-taul)*W*zeta_j*e * (gama/(1-gama)) * (1-l)
      // (4) (1)+(3) --> (1-taul)*W*zeta_j*e*(gama/(1-gama)) * (1-l) + ap = (1-taul)*W*l*zeta_j*e + yk - gta
      // (5) (4) --> (1-taul)*W*zeta_j*e*((gama/(1-gama))*(1-l)-l)* = yk - gta - ap
      // (6) (5) --> (1-taul)*W*zeta_j*e*(gama/(1-gama)) - (1-taul)*W*zeta_j*e*(gama/(1-gama))*l - (1-taul)*W*zeta_j*e*l = yk - gta - ap
      // (7) (6) --> l = ((1-taul)*W*zeta_j*e*(gama/(1-gama)) - yk + gta + ap)/((1+(gama/(1-gama)))*(1-taul)*W*zeta_j*e*(gama/(1-gama)))			  
      double tmp2 = (1.0-taul[ie])*W*zeta[j]*e_grid[ie];      
      *l = (tmp2*gama/(1.0-gama)-inc_nol_post_saving_evasion)/(tmp2*(1.0+gama/(1.0-gama)));
      *c = (1.0 - (*l))*tmp2*(gama/(1.0-gama))/(1.0+tauc);
      
      if( (*l) < 0.0)
	{
	  *l=0.0;
	  *c = (inc_nol_post_saving_evasion)/(1.0+tauc);
	}
      else if( (*l) > 1.0)
	{
	  fprintf(logfile,"Error l>1!\n");
	  free_mem();
	  exit(1);
	  *l=1.0;
	  *c = (inc_nol_post_saving_evasion + tmp2)/(1.0+tauc);
	}
    }
  else
    {
      *c = (inc_nol_post_saving_evasion)/(1.0+tauc);
    }
}

/****************************************************************/
/*   Functions related to updating/solving the Bellman equation */
/****************************************************************/

double bellman_eqn_on_grid(int iap, int ivp, void * p)
{
  bellman_params * p2 = (bellman_params *)p;
  int j = p2->j;
  int ie = p2->ie;
  int iz = p2->iz;
  int ii = p2->ii;
  double W = p2->W;
  double evasion_flow = p2->evasion_flow;
  double inc_nol_post_saving_evasion = p2->inc_nol_post_saving_evasion;
  gsl_spline ******* spline = p2->spline;

  double u =0.0;
  double c = 0.0;
  double l = 0.0;
  calc_c_l(j, ie, W, inc_nol_post_saving_evasion, &c, &l);
  u = util2(c,l);
  p2->c=c;
  p2->l=l;
  
  double EV = 0.0;
  double pd = prob_detect(ii,evasion_flow);
  
  if(j<(R-1))
    {
      for(int iep=0; iep<NE; iep++)
	{
	  for(int iip=0; iip<NI; iip++)
	    {
	      double tmp = spline[j+1][iep][iz][iip][ivp][0]->y[iap];
	      if(pd < 1.0e-10)
		{
		  EV += phi[j]*e_probs[ie][iep]*i_probs[ii][iip]*tmp;
		}
	      else
		{
		  EV += (1.0-pd)*phi[j]*e_probs[ie][iep]*i_probs[ii][iip]*tmp;

		  double tmp2 = spline[j+1][iep][iz][iip][ivp][1]->y[iap];
		  EV += pd*phi[j]*e_probs[ie][iep]*i_probs[ii][iip]*tmp2;
		}
	    }
	}
    }
  else if(j<(J-1))
    {
      
      for(int iip=0; iip<NI; iip++)
	{
	  double tmp = spline[j+1][ie][iz][iip][ivp][0]->y[iap];

	  if(pd < 1.0e-10)
	    {
	      EV += phi[j]*i_probs[ii][iip]*tmp;
	    }
	  else
	    {
	      EV += (1.0-pd)*phi[j]*i_probs[ii][iip]*tmp;
	      
	      double tmp2 = spline[j+1][ie][iz][iip][ivp][1]->y[iap];
	      EV += pd*phi[j]*i_probs[ii][iip]*tmp2;
	    }
	    
	}
    }
     

  double V_new = u + beta*EV;

  if(gsl_isnan(V_new) || gsl_isinf(V_new))
    {
      fprintf(logfile,"\t\t\t%d %0.6f %0.6f %0.6f\n",j,inc_nol_post_saving_evasion, p2->c, p2->l);
      GSL_ERROR("Bellman equation is NaN or Inf!",GSL_ERANGE);
    }

  return -V_new;
}

double bellman_eqn_off_grid(double ap, int ivp, void * p)
{
  bellman_params * p2 = (bellman_params *)p;
  int j = p2->j;
  int ie = p2->ie;
  int iz = p2->iz;
  int ii = p2->ii;
  double W = p2->W;
  double inc_nol_post_saving_evasion = p2->inc_nol_post_saving_evasion;
  double evasion_flow = p2->evasion_flow;
  gsl_spline ******* spline = p2->spline;
  gsl_interp_accel * acc = p2->acc;
  
  double u =0.0;
  double c = 0.0;
  double l = 0.0;
  calc_c_l(j, ie, W, inc_nol_post_saving_evasion, &c, &l);
  u = util2(c,l);
  p2->c=c;
  p2->l=l;

  double pd = prob_detect(ii,evasion_flow);
  double EV = 0.0;

  if(j<(R-1))
    {
      for(int iep=0; iep<NE; iep++)
	{
	  for(int iip=0; iip<NI; iip++)
	    {
	      double tmp = gsl_spline_eval(spline[j+1][iep][iz][iip][ivp][0],ap,acc);

	      if(pd<1.0e-10)
		{
		  EV += phi[j]*e_probs[ie][iep]*i_probs[ii][iip]*tmp;
		}
	      else
		{
		  EV += (1.0-pd)*phi[j]*e_probs[ie][iep]*i_probs[ii][iip]*tmp;
		  
		  double tmp2 = gsl_spline_eval(spline[j+1][iep][iz][iip][ivp][1],ap,acc);
		  EV += pd*phi[j]*e_probs[ie][iep]*i_probs[ii][iip]*tmp2;

		}
	    }
	}
    }
  else if(j<(J-1))
    {
      for(int iip=0; iip<NI; iip++)
	{
	  double tmp = gsl_spline_eval(spline[j+1][ie][iz][iip][ivp][0],ap,acc);

	  if(pd<1.0e-10)
	    {
	      EV += phi[j]*i_probs[ii][iip]*tmp;
	    }
	  else
	    {
	      EV += (1.0-pd)*phi[j]*i_probs[ii][iip]*tmp;
	      
	      double tmp2 = gsl_spline_eval(spline[j+1][ie][iz][iip][ivp][1],ap,acc);
	      EV += pd*phi[j]*i_probs[ii][iip]*tmp2;
	    }
	}
    }

  double V_new = u + beta*EV;

  if(gsl_isnan(V_new) || gsl_isinf(V_new))
    {
      GSL_ERROR("Bellman equation is NaN or Inf!",GSL_ERANGE);
    }

  return -V_new;
}

double bdiffs[NTH];
double fine_grid[NTH][FINE_GRID_SIZE];

int iterate_bellman(eqm_t * et, eqm_t * etp, double * supnorm, int fine_opt, int v)
{
  int abort=0;
  
  *supnorm = 0.0;

  SET_ALL_V(bdiffs,NTH,0.0);
  int j;
  for(j=(J-1); j>=0; j--)
    {
      if(j==(J-1))
	{
	  init_splines(etp,0,v);
	}
	else
	{
	  init_splines(etp,j+1,v);
	}
      
      int ie_iz;
      omp_set_num_threads(NTH);
#pragma omp parallel for private(ie_iz)
      for(ie_iz=0; ie_iz<NE*NZ; ie_iz++)
	{
	  int tn=0;
	  tn = omp_get_thread_num();
	  
	  int ie = ie_iz / ((int)NZ);
	  int iz = ie_iz % ((int)NZ);
	  if(ie>NE || iz > NZ)
	    {
	      fprintf(logfile,"Bad indices!\n");
	      abort=1;
	    }

	  int ii;
	  for(ii=0; ii<NI; ii++)
	    {
	      int iv;
	      int iv_ub = (evasion_type==1 && v==1 ? NV : 1);
	      for(iv=0; iv<iv_ub; iv++)
		{
		  int id;
		  int id_ub = (evasion_type==1 && iv>0) ? ND : 1;
		  for(id=0; id<id_ub; id++)
		    {

		      int ivd = (iv>0 && id==1) ? 0 : iv;
			
		      gsl_interp_accel_reset(acc[tn]);
	      
		      int ia=0;
		      int cnt;

		      double lower_bound = 0.0;
		      double upper_bound = et->a_grid[NA-1];

		      for(cnt=0; cnt<NA; cnt++)
			{
			  if(!abort)
			    {
			      if(cnt % 2==0)
				{
				  ia=cnt/2;
				  if(ia>0)
				    {
				      et->ga[v][j][ie][iz][ii][iv][id][ia] = et->ga[v][j][ie][iz][ii][iv][id][ia-1]+1.0e-6;
				    }
				}
			      else
				{
				  ia=NA-1-cnt/2;
				  if(ia<(NA-1))
				    {
				      et->ga[v][j][ie][iz][ii][iv][id][ia] = et->ga[v][j][ie][iz][ii][iv][id][ia+1]-1.0e-4;
				    }
			      
				}

			      double z=z_grid[iz];
			      if(ii==1)
				{
				  z = 0.0;
				}
			      double at = et->a_grid[ia]*
				(1.0-evasion_grid[iv]);
			      double k = et->gk[iz][ii][iv][ia];
			      double business_inc = et->P*pow(z*k,nu) -
				delta*k - et->r*fmax(k - at, 0.0);
			      double gross_inc_nol = et->gyk[iz][ii][iv][ia] + et->gtk[iz][ii][iv][ia];
			      
			      double inc_nol = et->lump_sum + et->a_grid[ia] + et->gyk[iz][ii][iv][ia] - et->gta[iv][ia];			      

			      double max_inc = et->lump_sum + et->a_grid[ia] + et->gyk[iz][ii][iv][ia] - et->gta[iv][ia];
			      if(j<R)
				{
				  max_inc += (1.0-taul[ie])*et->W*zeta[j]*e_grid[ie];
				}
			      else
				{
				  max_inc += et->ylbar*Phi[ie];
				  inc_nol += et->ylbar*Phi[ie];
				}

			      if(id==1)
				{
				  double ktax_evaded_penalty = calc_ktax_evaded(et->gtk[iz][ii][iv][ia],iv,et->gv[v][j][ie][iz][ii][iv][0][ia],
										et->a_grid[ia],et->ga[v][j][ie][iz][ii][iv][0][ia]);

				  double penalty = penalty_frac_stock * et->a_grid[ia]*evasion_grid[iv];
				  penalty += penalty_frac_tauk * fmin(max_penalty_years,j+1) * ktax_evaded_penalty;

				  if(wealth_tax_type==5)
				    penalty += penalty_frac_taua * et->gta_evaded[iv][ia];
				  else
				    penalty += penalty_frac_taua * fmin(max_penalty_years,j+1) * et->gta_evaded[iv][ia];

				  penalty = fmin(et->a_grid[ia]-et->gta[iv][ia],penalty);
				  
				  inc_nol -= penalty;
				  max_inc -= penalty;

				  if(gsl_isnan(penalty) || gsl_isinf(penalty) || inc_nol - et->lump_sum < -1.0e-10)
				    {
				      
				      fprintf(logfile,"\t\t****%d %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f****\n",
					     j, penalty,et->gta_evaded[iv][ia],inc_nol,et->a_grid[ia],
					     max_inc,et->gyk[iz][ii][iv][ia]);
				      abort=1;
				      break;
				      //GSL_ERROR("Invalid penalty!",GSL_ERANGE);
				    }

				}

			      double upper_bound2 = fmin(upper_bound,max_inc-1.0e-12);
			      bellman_params params = {j,ie,iz,ii,et->W,
						       gross_inc_nol,
						       business_inc,inc_nol,
						       0.0,0.0,0.0,0.0,
						       spline_V[tn],acc[tn]};

			      if(j==(J-1))
				{
				  et->ga[v][j][ie][iz][ii][iv][id][ia]=0.0;
				  et->gv[v][j][ie][iz][ii][iv][id][ia]=0;
				  et->V[v][j][ie][iz][ii][iv][id][ia] = util1(inc_nol/(1.0+tauc));
				}
			      else
				{	      
				  int iap=0;
				  int ivp=0;
				  int iapmax=0;
				  int ivpmax=0;
				  double vmax  = -HUGE_VAL;
				  for(ivp=0; ivp<(v==1 && id==0 ? NV : 1); ivp++)
				    {
				      for(iap=(ivp==0 ? 0 : 1); iap<NA; iap++)
					{
					  double ap = et->a_grid[iap];

					  double ktax_evaded = calc_ktax_evaded(et->gtk[iz][ii][ivd][ia],ivd,ivp,et->a_grid[ia],ap);

					  double evasion_cost = (ivp>0 && ap>1.0e-10)*theta +
					    eta*pow(fmax(0.0,fabs(ap*evasion_grid[ivp]-et->a_grid[ia]*evasion_grid[ivd])),eta2);

					  double upper_bound3 = fmin(upper_bound,upper_bound2+ktax_evaded-evasion_cost);
				      
					  if(ap>upper_bound3)
					    break;

					  params.inc_nol_post_saving_evasion = inc_nol - ap + ktax_evaded - evasion_cost;
					  params.evasion_flow = ap*evasion_grid[ivp] - et->a_grid[ia]*evasion_grid[ivd];
						
					  double tmp = -bellman_eqn_on_grid(iap,ivp,&params);
					  if(tmp>vmax)
					    {
					      vmax=tmp;
					      iapmax=iap;
					      ivpmax=ivp;
					      et->ga[v][j][ie][iz][ii][iv][id][ia] = et->a_grid[iapmax];
					      et->gv[v][j][ie][iz][ii][iv][id][ia] = ivpmax;
					      et->V[v][j][ie][iz][ii][iv][id][ia] = vmax;
					    }
					}
				    }
				
				  if(fine_opt)
				    {
				      if(iapmax==0)
					{
					  linspace(et->a_grid[0],
						   et->a_grid[1],
						   FINE_GRID_SIZE,
						   fine_grid[tn]);

					}
				      else if(iapmax==NA-1)
					{
					  linspace(et->a_grid[NA-2],
						   et->a_grid[NA-1],
						   FINE_GRID_SIZE,
						   fine_grid[tn]);
			  
					}
				      else
					{
					  linspace(et->a_grid[iapmax-1],
						   et->a_grid[iapmax+1],
						   FINE_GRID_SIZE,
						   fine_grid[tn]);
					}

				      if(fine_grid[tn][0]<et->a_grid[0])
					{
					  fine_grid[tn][0] = et->a_grid[0];
					}
				      if(fine_grid[tn][FINE_GRID_SIZE-1]>
					 et->a_grid[NA-1])
					{
					  fine_grid[tn][FINE_GRID_SIZE-1]=
					    et->a_grid[NA-1];
					}
		      
				      int iap2;
				      for(iap2=1; iap2<FINE_GRID_SIZE-1; iap2++)
					{
					  if(fine_grid[tn][iap2]<et->a_grid[0] ||
					     fine_grid[tn][iap2]>et->a_grid[NA-1])
					    {
					      fprintf(logfile,"Fine grid constructed incorrectly!\n");
					      abort=1;
					      break;
					    }

					  double ap = fine_grid[tn][iap2];
					  double ktax_evaded = calc_ktax_evaded(et->gtk[iz][ii][ivd][ia],ivd,ivpmax,et->a_grid[ia],ap);
					  
					  double evasion_cost = (ivpmax>0 && ap>1.0e-10)*theta +
					    eta*pow(fmax(0.0,fabs(ap*evasion_grid[ivpmax]-et->a_grid[ia]*evasion_grid[ivd])),eta2);

					  double upper_bound3 = fmin(upper_bound,upper_bound2+ktax_evaded-evasion_cost);
				  
					  if(ap>upper_bound3)
					    break;

					  params.inc_nol_post_saving_evasion = inc_nol - ap + ktax_evaded - evasion_cost;
					  params.evasion_flow = ap*evasion_grid[ivpmax]-et->a_grid[ia]*evasion_grid[ivd];
				      
					  double tmp = -bellman_eqn_off_grid(ap,ivpmax,&params);
					  if(tmp>vmax)
					    {
					      et->ga[v][j][ie][iz][ii][iv][id][ia] = ap;
					      vmax=tmp;
					      et->V[v][j][ie][iz][ii][iv][id][ia] = vmax;
					    }
					}

				      if(et->ga[v][j][ie][iz][ii][iv][id][ia]-fine_grid[tn][FINE_GRID_SIZE-1]>1.0e-8 ||
					 et->ga[v][j][ie][iz][ii][iv][id][ia]-fine_grid[tn][0]<-1.0e-8)
					{
					  fprintf(logfile,"Error! Continuous optimization went outside of bounds %0.6f %0.6f %0.6f.\n",
					  et->ga[v][j][ie][iz][ii][iv][id][ia],fine_grid[tn][0],fine_grid[tn][FINE_GRID_SIZE-1]);
					  abort=1;
					  break;
					}
				    }
				}

			  
			      int iap0 = gsl_interp_accel_find(acc[tn],et->a_grid,NA,et->ga[v][j][ie][iz][ii][iv][id][ia]);
			      et->giap0[v][j][ie][iz][ii][iv][id][ia] = iap0;

			      
			      if(j==0)
				{
				  double diff = fabs(et->V[v][j][ie][iz][ii][iv][id][ia]-
						     spline_V[tn][j][ie][iz][ii][iv][id]->y[ia]);
				  diff = diff - vf_tol_rel*fabs(et->V[v][j][ie][iz][ii][iv][id][ia]);
				
				  if(diff>bdiffs[tn])
				    {
				      bdiffs[tn] = diff;
				    }
				}

			      if(exploit_monotonicity)
				{
				  if(cnt % 2==0)
				    {
				      lower_bound =
					fmax(lower_bound,et->ga[v][j][ie][iz][ii][iv][id][ia]/bound_mult-2.0);
				    }
				  else
				    {
				      upper_bound =
					fmin(upper_bound,et->ga[v][j][ie][iz][ii][iv][id][ia]*bound_mult);
				    }
				}
			    } 
			}
		    }
		}
	    }
	}
#pragma omp flush(abort)
    }

  int tn;
  for(tn=0; tn<NTH; tn++)
    {
      if(bdiffs[tn]>*supnorm)
	{
	  *supnorm = bdiffs[tn];
	}
    }
           
  if(abort)
    {
      fprintf(logfile,"\nError updating value function!\n");
      return 1;		      
    }
  else
    {
      return 0;
    }
}

int solve_ss_bellman(eqm_t * et)
{
  time_t start, stop;
  time(&start);
  #ifndef __cplusplus
  if(verbose)
    fprintf(logfile,"\tSolving for steady-state value function\n");
  #endif

  double supnorm = +HUGE_VAL;
  int iter = 0;
  int status=0;

  fmin_ftol_rel = 1.0e-11;
  fmin_ftol_abs = 1.0e-11;
  fmin_xtol_rel = 1.0e-6;
  fmin_xtol_abs = 1.0e-6;
  bound_mult = 2.0;
  int fine_opt=1;

  set_capital_policies(et);
  
  //do
    {      
      time_t stop2;
      iter++;

      if(exploit_monotonicity)
	{
	  if(supnorm<1000.0*vf_tol_abs)
	    {
	      bound_mult = 1.8;
	    }
	  if(supnorm<100.0*vf_tol_abs)
	    {
	      bound_mult = 1.6;
	    }
	  if(supnorm<10.0*vf_tol_abs)
	    {
	      bound_mult = 1.4;
	    }
	  if(supnorm<5.0*vf_tol_abs)
	    {
	      bound_mult = 1.2;
	    }
	}

      double dn=0;
      if(iterate_bellman(et,et,&dn,fine_opt,0))
	{
	  status=1;
	}

      double dv=0;
      if(iterate_bellman(et,et,&dv,fine_opt,1))
	{
	  status=1;
	}

      supnorm = fmax(dn,dv);

      if(iter % 1 == 0 && 0)
	{
	  time(&stop2);
	  fprintf(logfile,"\t\tIter %d, opt type = %d, time = %0.2f secs, ||V-V'|| = %0.8g\n",iter,fine_opt,difftime(stop2,start),supnorm);
	}
    }
  //while(supnorm>vf_tol_abs && iter<vf_max_iter);

  time(&stop);

  if(status==1)
    {
      fprintf(logfile,"\tVF failed to converge! Time = %0.2f\n",difftime(stop,start));
      return 1;
    }
  else
    {
      if(verbose)
	fprintf(logfile,"\tVF converged! Time = %0.0f, ||V-V'|| = %0.6e\n",
		difftime(stop,start),supnorm);
      return 0;
    }
}

/****************************************************************/
/*    Functions related to updating the distribution            */
/****************************************************************/

void init_dist(eqm_t * et)
{
  double sum3=0.0;
  
  SET_ALL_V(et->Psi,2*J*NE*NZ*NI*NV*ND*NA,0.0);
  int v, j, ie, iz, ii;
  for(ie=0; ie<NE; ie++)
    {
      for(iz=0; iz<NZ; iz++)
	{
	  et->Psi[1][0][ie][iz][0][0][0][0] = et->frac_v*e_ergodic_dist[ie]*z_ergodic_dist[iz]*pi;
	  et->Psi[1][0][ie][iz][1][0][0][0] = et->frac_v*e_ergodic_dist[ie]*z_ergodic_dist[iz]*(1.0-pi);
	  sum3+= et->Psi[1][0][ie][iz][0][0][0][0] + et->Psi[1][0][ie][iz][1][0][0][0];

	  et->Psi[0][0][ie][iz][0][0][0][0] = (1.0-et->frac_v)*e_ergodic_dist[ie]*z_ergodic_dist[iz]*pi;
	  et->Psi[0][0][ie][iz][1][0][0][0] = (1.0-et->frac_v)*e_ergodic_dist[ie]*z_ergodic_dist[iz]*(1.0-pi);
	  sum3+= et->Psi[0][0][ie][iz][0][0][0][0] + et->Psi[0][0][ie][iz][1][0][0][0];

	}
    }
  
  double sum=1.0;
  double sum2=0.0;

  for(j=1; j<J; j++)
    {
      sum2=sum;
      sum=0.0;
      
      for(v=0; v<2; v++)
	{
	  for(ie=0; ie<NE; ie++)
	    {
	      for(iz=0; iz<NZ; iz++)
		{
		  for(ii=0; ii<NI; ii++)
		    {
		      et->Psi[v][j][ie][iz][ii][0][0][0] = et->Psi[v][j-1][ie][iz][ii][0][0][0]*phi[j-1];
		      sum3 += et->Psi[v][j][ie][iz][ii][0][0][0];
		      sum += et->Psi[v][j][ie][iz][ii][0][0][0];
		      
		    }
		}
	    }
	}
	  
      if(fabs(sum-phi[j-1]*sum2)>1.0e-9)
	{
	  fprintf(logfile,"Error! sum = %0.2f, phi*sum2 = %0.2f\n",sum,phi[j-1]*sum2);
	  double x=9;
	}
    }

}

void update_dist(eqm_t * et, double * supnorm, int check_conv)
{
  time_t start, stop;
  time(&start);  
  
#pragma omp parallel for
  for(int j=0; j<J; j++)
  {
    SET_ALL_V(tmp_Psi[0][j],NE*NZ*NI*NV*ND*NA,0.0);
    SET_ALL_V(tmp_Psi[1][j],NE*NZ*NI*NV*ND*NA,0.0);
  }
  
#pragma omp parallel
  {
    double (*tmp_Psi_private)[J][NE][NZ][NI][NV][ND][NA] = malloc(sizeof(double[2][J][NE][NZ][NI][NV][ND][NA]));
    SET_ALL_V(tmp_Psi_private,2*J*NE*NZ*NI*NV*ND*NA,0.0);
    
#pragma omp for
    for(int j=0; j<J; j++)
      {
	for(int v=0; v<2; v++)
	  {
	    for(int ie_iz=0; ie_iz<NE*NZ; ie_iz++)
	      {	    
		int ie = ie_iz / ((int)NZ);
		int iz = ie_iz % ((int)NZ);
	  
		double phi_ = phi[j];

		for(int ii=0; ii<NI; ii++)
		  {
		    int iv_ub = (evasion_type==1 && v==1 ? NV : 1);
		    for(int iv=0; iv<iv_ub; iv++)
		      {
			int id_ub = (evasion_type==1 && iv>0) ? ND : 1;	    
			for(int id=0; id<id_ub; id++)
			  {
			    for(int ia=0; ia<NA; ia++)
			      {
				double psi = et->Psi[v][j][ie][iz][ii][iv][id][ia];
				int iap0 = et->giap0[v][j][ie][iz][ii][iv][id][ia];
				int iap1 = iap0+1;
				double frac1 = (et->ga[v][j][ie][iz][ii][iv][id][ia]-et->a_grid[iap0])/(et->a_grid[iap1]-et->a_grid[iap0]); 
				double frac0 = 1.0-frac1;
				int ivp = et->gv[v][j][ie][iz][ii][iv][id][ia];
				
				int ivd = (iv>0 && id==1) ? 0 : iv;
				double evasion_flow = et->ga[v][j][ie][iz][ii][iv][id][ia]*evasion_grid[ivp] - et->a_grid[ia]*evasion_grid[ivd];
				double pd = prob_detect(ii,evasion_flow);
				
				if(j<(R-1))
				  {
				    for(int iep = 0; iep<NE; iep++)
				      {
					for(int iip=0; iip<NI; iip++)
					  {
					    tmp_Psi_private[v][j+1][iep][iz][iip][ivp][0][iap0] += psi *
					      frac0*phi_*(1.0-pd)*e_probs[ie][iep]*i_probs[ii][iip];
					
					    tmp_Psi_private[v][j+1][iep][iz][iip][ivp][0][iap1] += psi *
					      frac1*phi_*(1.0-pd)*e_probs[ie][iep]*i_probs[ii][iip];
					
					    tmp_Psi_private[v][j+1][iep][iz][iip][ivp][1][iap0] += psi *
					      frac0*phi_*pd*e_probs[ie][iep]*i_probs[ii][iip];
					
					    tmp_Psi_private[v][j+1][iep][iz][iip][ivp][1][iap1] += psi *
					      frac1*phi_*pd*e_probs[ie][iep]*i_probs[ii][iip];
					
					  }
				      }
				  }
			    
				else if(j<(J-1))
				  {
				    for(int iip=0; iip<NI; iip++)
				      {
					tmp_Psi_private[v][j+1][ie][iz][iip][ivp][0][iap0] += psi*
					  frac0*phi_*(1.0-pd)*i_probs[ii][iip];
				    
					tmp_Psi_private[v][j+1][ie][iz][iip][ivp][0][iap1] += psi*
					  frac1*phi_*(1.0-pd)*i_probs[ii][iip];
				    
					tmp_Psi_private[v][j+1][ie][iz][iip][ivp][1][iap0] += psi*
					  frac0*phi_*pd*i_probs[ii][iip];
				    
					tmp_Psi_private[v][j+1][ie][iz][iip][ivp][1][iap1] += psi*
					  frac1*phi_*pd*i_probs[ii][iip];
				      }
				  }
			    
				for(int iep = 0; iep<NE; iep++)
				  {
				    for(int izp=0; izp<NZ; izp++)
				      {
				
					tmp_Psi_private[0][0][iep][izp][0][0][0][iap0] += psi*
					  frac0*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*pi*(1.0-et->frac_v);
				    
					tmp_Psi_private[0][0][iep][izp][0][0][0][iap1] += psi*
					  frac1*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*pi*(1.0-et->frac_v);
				    
					tmp_Psi_private[0][0][iep][izp][1][0][0][iap0] += psi*
					  frac0*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*(1.0-pi)*(1.0-et->frac_v);
				    
					tmp_Psi_private[0][0][iep][izp][1][0][0][iap1] += psi*
					  frac1*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*(1.0-pi)*(1.0-et->frac_v);

					tmp_Psi_private[1][0][iep][izp][0][0][0][iap0] += psi*
					  frac0*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*pi*(et->frac_v);
				    
					tmp_Psi_private[1][0][iep][izp][0][0][0][iap1] += psi*
					  frac1*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*pi*(et->frac_v);
				    
					tmp_Psi_private[1][0][iep][izp][1][0][0][iap0] += psi*
					  frac0*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*(1.0-pi)*(et->frac_v);
				    
					tmp_Psi_private[1][0][iep][izp][1][0][0][iap1] += psi*
					  frac1*(1.0-phi_)*e_birth_probs[ie][iep]*
					  z_probs[iz][izp]*(1.0-pi)*(et->frac_v);

				      }
				  }
			      }
			  }
		      }
		  }
	      }
	  }
      }
    
#pragma omp critical
    {
      for(int v=0; v<2; v++)
	{
	  for(int j=0; j<J; j++)
	    {
	      int ie;
	      for(ie=0; ie<NE; ie++)
		{
		  int iz;
		  for(iz=0; iz<NZ; iz++)
		    {
		      int ii;
		      for(ii=0; ii<NI; ii++)
			{
			  int iv;
			  int iv_ub = (evasion_type==1 && v==1 ? NV : 1);
			  for(iv=0; iv<iv_ub; iv++)
			    {
			      int id;
			      int id_ub = (evasion_type==1 && iv>0) ? ND : 1;
			      for(id=0; id<id_ub; id++)
				{
				  int ia;
				  for(ia=0; ia<NA; ia++)
				    {
				      tmp_Psi[v][j][ie][iz][ii][iv][id][ia] += tmp_Psi_private[v][j][ie][iz][ii][iv][id][ia];
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }

    free(tmp_Psi_private);
  }

  if(check_conv)
    {
      *supnorm = 0.0;
      double test0=0.0;
      double test1=0.0;

      for(int v=0; v<2; v++)
	{
	  for(int j=0; j<J; j++)
	    {
	      int ie;
	      for(ie=0; ie<NE; ie++)
		{
		  int iz;
		  for(iz=0; iz<NZ; iz++)
		    {
		      int ii;
		      for(ii=0; ii<NI; ii++)
			{
			  int iv;
			  int iv_ub = (evasion_type==1 && v==1 ? NV : 1);
			  for(iv=0; iv<iv_ub; iv++)
			    {
			      int id;
			      int id_ub=(evasion_type==1 && iv>0) ? ND : 1;
			      for(id=0; id<id_ub; id++)
				{
				  int ia;
				  for(ia=0; ia<NA; ia++)
				    {
				  
				      double tmp = fabs(tmp_Psi[v][j][ie][iz][ii][iv][id][ia] -
							et->Psi[v][j][ie][iz][ii][iv][id][ia]);
				  
				      test0 += (j==0)*et->Psi[v][j][ie][iz][ii][iv][id][ia];
				      test1 += (j==0)*tmp_Psi[v][j][ie][iz][ii][iv][id][ia];
					
				      if(tmp>*supnorm)
					{
					  *supnorm = tmp;
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}

	  
      if(fabs(test0-test1)>1.0e-5)
	{
	  fprintf(logfile,"Unstable household mass! Phi0[0] = %0.6f, Phi1[0] = %0.6f\n",test0,test1);
	}
    }

  time(&stop);
  //fprintf(logfile,"\t\t\t%0.1f\n",difftime(stop,start));
}
	 

int solve_ss_dist(eqm_t * et)
{
  time_t start, stop;
  time(&start);
  #ifndef __cplusplus
  if(verbose)
    fprintf(logfile,"\tSolving for steady-state distribution\n");
  #endif

  double supnorm = +HUGE_VAL;
  int iter = 0;

  if(dist_max_iter>0)
    {
      do
	{      
	  time_t stop2;
	  iter++;

	  if(iter>=200 && iter % 10==0)
	    update_dist(et,&supnorm,1);
	  else
	    update_dist(et,&supnorm,0);
	    
	  memcpy( (double *)et->Psi,
		  (double *)tmp_Psi,
		  //(double *)tmp_Psi2[0],
		  sizeof(double)*2*J*NE*NZ*NI*NV*NA*ND );

	  if(iter>=200 && iter % 10 == 0 && 1)
	    {
	      time(&stop2);
	      fprintf(logfile,"\t\tIter %d, time = %0.2f secs, ||D-D'|| = %0.2e\n",
		      iter,difftime(stop2,start),supnorm);
	    }
	}
      while(supnorm>dist_tol && iter<dist_max_iter);
      
      time(&stop);

      if(iter==dist_max_iter)
	{
	  fprintf(logfile,"\tDistribution failed to converge!\n Time = %0.2f\n",
		  difftime(stop,start));
	  return 1;
	}
      else
	{
	  if(verbose)
	    fprintf(logfile,"\tDist converged! Time = %0.2f, iter = %d, ||D-D'|| = %0.2e\n",
		    difftime(stop,start),iter,supnorm);
      
	  return 0;
	}
    }
  else
    {
      return 0;
    }
}

/****************************************************************/
/*    Functions related to aggregation                          */
/****************************************************************/

double excess_k_demand(double r, void * params)
{
  eqm_t * et = (eqm_t *)params;
  
  et->r = r;
  set_capital_policies(et);
  
  double Q = 0.0;
  double Kd = 0.0;
  double Ks = 0.0;
  int v, j, ie, iz, ii, ia, iv, id;
  for(v=0; v<2; v++)
    {
      for(j=0; j<J; j++)
	{
	  for(ie=0; ie<NE; ie++)
	    {
	      for(iz=0; iz<NZ; iz++)
		{
		  for(ii=0; ii<NI; ii++)
		    {
		      for(iv=0; iv< (evasion_type==1 && v==1 ? NV : 1); iv++)
			{
			  for(id=0; id< (evasion_type==1 && iv>0 ? 2 : 1); id++)
			    {
			      for(ia=0; ia<NA; ia++)
				{
				  Kd += et->Psi[v][j][ie][iz][ii][iv][id][ia]*et->gk[iz][ii][iv][ia];
				  Ks += et->Psi[v][j][ie][iz][ii][iv][id][ia]*et->gks[iz][ii][iv][ia];

				  double z = z_grid[iz];
				  if(ii==1)
				    {
				      z=0.0;
				    }
				  Q += et->Psi[v][j][ie][iz][ii][iv][id][ia]*pow(z*et->gk[iz][ii][iv][ia],nu);
				}
			    }
			}
		    }
		}
	    }
	}
    }

  Q = pow(Q,1.0/nu);
  double tmp = ((et->r+delta)/alpha2)*pow(Q,-alpha)*pow(et->L,alpha+alpha2-1.0);
  tmp = pow(tmp,1.0/(alpha2-1.0));
  Kd += tmp;
  double diff = Kd-Ks;
  
  return diff;
}

int clear_k_market(eqm_t * et, double * r_new)
{
  double r0 = et->r;
  double r = r0;
  double r_lo, r_hi;
  
  double diff = excess_k_demand(r0,et);
  if(diff>0.0)
    {
      r_lo=r0;
      r_hi=r0+0.005;

      if(excess_k_demand(r_hi,et)>0.0)
	{
	  r_hi += 0.005;
	}
      if(excess_k_demand(r_hi,et)>0.0)
	{
	  r_hi += 0.01;
	}
      if(excess_k_demand(r_hi,et)>0.0)
	{
	  r_hi += 0.015;
	}
      if(excess_k_demand(r_hi,et)>0.0)
	{
	  r_hi += 0.02;
	}
      if(excess_k_demand(r_hi,et)>0.0)
	{
	  r_hi += 0.03;
	}
      if(excess_k_demand(r_hi,et)>0.0)
	{
	  r_hi += 0.05;
	}

    }
  else
    {
      r_hi=r0;
      r_lo=r0-0.005;

      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.005;
	}
      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.005;
	}
      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.0025;
	}
      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.0025;
	}
      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.0015;
	}
      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.0015;
	}
      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.00075;
	}
      if(excess_k_demand(r_lo,et)<0.0)
	{
	  r_lo -= 0.00075;
	}

    }
  
  gsl_function F;
  F.function = &excess_k_demand;
  F.params = et;

  gsl_root_fsolver * s = gsl_root_fsolver_alloc (gsl_root_fsolver_brent);
  if(gsl_root_fsolver_set (s, &F, r_lo, r_hi))
    {
      return 1;
    }

  int iter=0;
  int status=0;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate (s);
      r = gsl_root_fsolver_root (s);
      r_lo = gsl_root_fsolver_x_lower (s);
      r_hi = gsl_root_fsolver_x_upper (s);
      status = gsl_root_test_interval (r_lo, r_hi, root_tol_abs, root_tol_rel);
    }
  while (status == GSL_CONTINUE && iter < max_root_iter);

  if(status != GSL_SUCCESS)
    {
      fprintf(logfile,"Failed to find capital market-clearing interest rate!\n");
    }
  
  gsl_root_fsolver_free(s);
  
  et->r = r0;
  set_capital_policies(et);
  *r_new = r;
  
  return status;
}


int compare_fun(const void * a, const void * b)
{
  double * da = (double *)a;
  double * db = (double *)b;
  return da[0]-db[0];
}

int aggregate(eqm_t * et, double update_speed, int pe, int verbose, int clear_k_mkt_flag, double * diff)
{

  time_t start, stop;
  time(&start);
	 
  double sumA = 0.0;
  double sumA_rep = 0.0;
  double sumA_hid = 0.0;
  double L = 0.0;
  double Kd = 0.0;
  double Ks = 0.0;
  double Q = 0.0;
  double ylbar = 0.0;
  double Y = 0.0;
  double mass_not_retired = 0.0;
  double total_mass = 0.0;
  double mass_newborn = 0.0;
  double ltaxR=0.0;
  double ktaxR = 0.0;
  double wtaxR = 0.0;
  double ktaxR_lost = 0.0;
  double wtaxR_lost = 0.0;
  double ctaxR = 0.0;
  double SS = 0.0;
  double wealth_tax_mass = 0.0;
  double C = 0.0;
  double EI = 0.0;
  double welfare=0.0;
  double welfare_newborn=0.0;
  double approval=0.0;
  double approval_newborn=0.0;
  double mass_w_tax_shelter=0.0;
  
  double detection_penalty = 0.0;
  double detection_mass = 0.0;
    
  int v, j, ie, iz, ii, ia, iv, id;
  for(v=0; v<2; v++)
    {
      for(j=0; j<J; j++)
	{
	  for(ie=0; ie<NE; ie++)
	    {
	      for(iz=0; iz<NZ; iz++)
		{
		  for(ii=0; ii<NI; ii++)
		    {
		      for(iv=0; iv<(evasion_type==1 && v==1 ? NV : 1); iv++)
			{
			  for(id=0; id<(evasion_type ==1 && iv>0 ? 2 : 1); id++)
			    {
			      int ivd = (iv>0 && id==1) ? 0 : iv;
			  
			      for(ia=0; ia<NA; ia++)
				{
				  double psi = et->Psi[v][j][ie][iz][ii][iv][id][ia];

				  int ivp = et->gv[v][j][ie][iz][ii][iv][id][ia];
				  double ap = et->ga[v][j][ie][iz][ii][iv][id][ia];
			      
				  double z=z_grid[iz];
				  if(ii==1)
				    {
				      z = 0.0;					      
				    }
			      
				  double ktax_evaded = calc_ktax_evaded(et->gtk[iz][ii][ivd][ia],ivd,et->gv[v][j][ie][iz][ii][iv][id][ia],et->a_grid[ia],et->ga[v][j][ie][iz][ii][iv][id][ia]);
				  //double wtax_evaded = et->gta_evaded[ivd][ia];

				  double inc_nol = et->lump_sum + et->a_grid[ia] + et->gyk[iz][ii][iv][ia] - et->gta[iv][ia];			      
				  if(j>=R)
				    {
				      inc_nol += et->ylbar*Phi[ie];
				    }

				  if(id==1)
				    {
				      detection_mass += psi;
				      double ktax_evaded_penalty = calc_ktax_evaded(et->gtk[iz][ii][iv][ia],iv,et->gv[v][j][ie][iz][ii][iv][0][ia],
										    et->a_grid[ia],et->ga[v][j][ie][iz][ii][iv][0][ia]);

				      double penalty = penalty_frac_stock * et->a_grid[ia]*evasion_grid[iv];
				      penalty += penalty_frac_tauk * fmin(max_penalty_years,j+1) * ktax_evaded_penalty;

				      if(wealth_tax_type==5)
					penalty += penalty_frac_taua * et->gta_evaded[iv][ia];
				      else
					penalty += penalty_frac_taua * fmin(max_penalty_years,j+1) * et->gta_evaded[iv][ia];

				      penalty = fmin(penalty,et->a_grid[ia]-et->gta[iv][ia]);

				      inc_nol -= penalty;
				  
				      if(!detection_revenue_flag)
					detection_penalty += penalty*psi;
				    }

				  double evasion_cost = (ivp>0 && ap>1.0e-10)*theta +
				    eta*pow(fmax(0.0,fabs(ap*evasion_grid[ivp]-et->a_grid[ia]*evasion_grid[ivd])),eta2);

				  double inc_nol_post_saving_evasion = inc_nol - ap + ktax_evaded - evasion_cost;
				  double c=0.0;
				  double l=0.0;

				  if(j==(J-1))
				    {
				      c=inc_nol/(1.0+tauc);
				    }
				  else
				    {
				      calc_c_l(j,ie,et->W,inc_nol_post_saving_evasion,&c,&l);
				    }
				  double gtc = tauc*c;
				  double gtl = taul[ie] * et->W * zeta[j] * e_grid[ie] * l;
			      
				  sumA += et->a_grid[ia] * psi;
				  sumA_rep += et->a_grid[ia] * (1.0-evasion_grid[ivd]) * psi;
				  sumA_hid += et->a_grid[ia] * evasion_grid[ivd] * psi;
				  L += l * zeta[j] * e_grid[ie] * psi;
				  C += c * psi;
				  Kd += et->gk[iz][ii][iv][ia] * psi;
				  Ks += et->gks[iz][ii][iv][ia] * psi;
				  ylbar += et->W*l*zeta[j]*e_grid[ie] * psi;
			  
				  Q += pow(z*et->gk[iz][ii][ivd][ia],nu) * psi;
				  EI += et->P*pow(z*et->gk[iz][ii][ivd][ia],nu) * psi;
		      
				  total_mass += psi;

				  ltaxR += psi * gtl;
				  ktaxR += psi * et->gtk[iz][ii][ivd][ia];
				  ctaxR += psi * gtc;
				  wtaxR += psi * et->gta[ivd][ia];

				  ktaxR -= ktax_evaded * psi;
				  ktaxR_lost += ktax_evaded * psi;
				  wtaxR_lost += et->gta_evaded[ivd][ia]*psi;
		      
				  // compute wealth tax revenue
				  if(et->gta[ivd][ia]>1.0e-9)
				    {
				      wealth_tax_mass += psi;
				    }

				  if(ivd>0)
				    {
				      mass_w_tax_shelter += psi;
				    }
			  
				  if(j>=R)
				    {
				      SS += et->ylbar * Phi[ie] * psi;
				    }
				  else
				    {
				      mass_not_retired += psi;
				    }
				  
				  if(j==0)
				    mass_newborn += psi;

				  welfare += et->V[v][j][ie][iz][ii][iv][id][ia]*psi;
				  if(j==0)
				    {
				      welfare_newborn += et->V[v][j][ie][iz][ii][iv][id][ia]*psi;
				    }

				  if(!benchmark_eqm)
				    {
				      if(et->V[v][j][ie][iz][ii][iv][id][ia]>ss0->V[v][j][ie][iz][ii][iv][id][ia])
					{
					  approval += ss0->Psi[v][j][ie][iz][ii][iv][id][ia];
					  if(j==0)
					    {
					      approval_newborn += ss0->Psi[v][j][ie][iz][ii][iv][id][ia];
					    }
					}
				    }
				}
			    }
			}
		    }
	      
		}
	    }
	}
    }

  if(gsl_isnan(L) || gsl_isinf(L) || fabs(L)<1.0e-5)
    {
      L=et->L;
    }
  
  Q = pow(Q,1.0/nu);
  double K = ((et->r+delta)/alpha2)*pow(Q,-alpha)*pow(L,alpha+alpha2-1.0);
  K = pow(K,1.0/(alpha2-1.0));
  Y = pow(Q,alpha)*pow(K,alpha2)*pow(L,1.0-alpha-alpha2);
  ylbar = ylbar/mass_not_retired;
  
  welfare = welfare/total_mass;
  welfare_newborn = welfare_newborn/mass_newborn;
  approval = approval/total_mass;
  approval_newborn = approval_newborn/mass_newborn;
  
  double G = ltaxR+ktaxR+ctaxR+wtaxR+detection_penalty-SS;
  double lump_sum = (G-et->G)/total_mass;

  calc_wealth_tax_evasion_elast(0,et);
  
  double Wp = (1.0-alpha-alpha2)*pow(Q,alpha)*pow(K,alpha2)*pow(L,-alpha-alpha2);
  double rp=et->r;
  double tauap=taua;

  taul_bar=ltaxR/(et->W*et->L);
  double taulp=taul_bar;

  if(taua_clear_gbc)
    {
      tauap = (SS+et->G-ltaxR-ktaxR-ctaxR-detection_penalty)/sumA_rep;
    }
  else if(taul_clear_gbc)
    {
      taulp = (SS + et->G - ktaxR - ctaxR - wtaxR - detection_penalty)/(et->W*et->L);
    }

  if(fabs(Kd+K-Ks)>root_tol_abs+root_tol_rel*fmin(Kd+K,Ks) && clear_k_mkt_flag && pe==0)
    {
      if(clear_k_market(et,&rp))
	{
	  fprintf(logfile,"Failed to clear capital market!\n");
	  return 1;
	}
    }
  
  double diff1 = Wp-et->W;
  double diff2 = rp-et->r;
  double diff3 = 0.0;

  if(benchmark_eqm)
    {
      diff3=0.0;
    }
  else if(taua_clear_gbc)
    {
      diff3=tauap-taua;
    }
  else if(taul_clear_gbc)
    {
      diff3=taulp-taul_bar;
    }
  else if(lump_sum)
    {
      diff3 = lump_sum-et->lump_sum;
    }

  et->detection_revenue = 100*detection_penalty/Y;
  et->detection_rate = 100*detection_mass/mass_w_tax_shelter;
    
  if(verbose==2)
    {
      fprintf(logfile,"\tQ  = %0.6g, Q' = %0.6g\n",et->Q,Q);
      fprintf(logfile,"\tL  = %0.6g, L' = %0.6g\n",et->L,L);
      fprintf(logfile,"\tK  = %0.6g, K' = %0.6g\n",et->K,K);
      fprintf(logfile,"\tW = %0.6g, W' = %0.6g\n",et->W,Wp);
      fprintf(logfile,"\tr = %0.6g, r' = %0.6g\n",et->r,rp);
      fprintf(logfile,"\tP = %0.6g, P' = %0.6g\n",et->P,alpha * pow(Q,alpha-nu) * pow(K,alpha2) * pow(L,1.0-alpha-alpha2));
      if(taua_clear_gbc)
	{
	  fprintf(logfile,"\ttaua = %0.6g, taua' = %0.6g\n",taua,tauap);
	}
      
      fprintf(logfile,"\tA = %0.6g\n",sumA);
      fprintf(logfile,"\tA (reported) = %0.6g\n",sumA_rep);
      fprintf(logfile,"\tA (hidden) = %0.6g\n",sumA_hid);
      fprintf(logfile,"\tY = %0.6g\n",Y);
      fprintf(logfile,"\tA_rep/Y = %0.6g\n",sumA_rep/Y);
      fprintf(logfile,"\tW/Y = %0.6g\n",et->W*L/Y);
      fprintf(logfile,"\tEI/Y = %0.6g\n",EI/Y);
      fprintf(logfile,"\tavg. labor income = %0.6g\n\n",ylbar);
      
      fprintf(logfile,"\tTotal tax revenue (pct Y) = %0.6g\n",100.*(ltaxR+ktaxR+ctaxR+wtaxR+detection_penalty)/Y);
      fprintf(logfile,"\tCapital income tax revenue (pct Y) = %0.6g\n",100*ktaxR/Y);
      fprintf(logfile,"\tWealth tax revenue (pct Y) = %0.6g\n",100*wtaxR/Y);
      fprintf(logfile,"\tLost capital income tax revenue (pct Y) = %0.6g\n",100*ktaxR_lost/Y);
      fprintf(logfile,"\tLost Wealth tax revenue (pct Y) = %0.6g\n",100*wtaxR_lost/Y);
      fprintf(logfile,"\tSS outlays = %0.6g\n",SS);
      fprintf(logfile,"\tPublic goods = %0.6g\n",G);
      if(lump_sum)
	fprintf(logfile,"\tLump sum = %0.6g\n",lump_sum);
      fprintf(logfile,"\tDetection penalty (pct Y) = %0.6g\n",100*detection_penalty/Y);     
      
      time(&stop);
      fprintf(logfile,"\n\tAggregation complete. Time = %0.0f\n",difftime(stop,start));

    }
  else if(verbose==1 && update_speed>1e-10)
    {
      if(pe==2)
	{
	  fprintf(logfile,"\tImplied prices: W = %0.4f, T = %0.6e\n",Wp,lump_sum);
	}
      else if(taua_clear_gbc)
	{
	  fprintf(logfile,"\tImplied prices: W = %0.4f, r = %0.2e, taua = %0.2e\n",Wp,rp,tauap);
	}
      else if(taul_clear_gbc)
	{
	  fprintf(logfile,"\tImplied prices: W = %0.4f, r = %0.2e, taul = %0.2e\n",Wp,rp,taulp);
	}
      else if(lump_sum)
	{
	  fprintf(logfile,"\tImplied prices: W = %0.4f, r = %0.2e, T = %0.6e\n",Wp,rp,lump_sum);
	}
      else
	{
	  fprintf(logfile,"\tImplied prices: W = %0.4f, r = %0.2e\n",Wp,rp);
	}
    }

  // update the aggregates
  if(pe!=1)
    {
      double W = update_speed*Wp + (1.0-update_speed)*et->W;
      double r = (update_speed)*rp + (1.0-update_speed)*et->r;
      double taua_ = (update_speed)*tauap + (1.0-update_speed)*taua;
      double taul_ = (update_speed)*taulp + (1.0-update_speed)*taul_bar;

      et->W=W;

      if(pe==0)
	{
	  et->r = r;
	}

      if(taua_clear_gbc)
	{
	  taua=taua_;
	}
      else if(taul_clear_gbc)
	{
	  for(int ie=0; ie<NE; ie++)
	    {
	      taul[ie] = (taul_/taul_bar) * taul[ie];
	    }
	}
      
      et->C = C;
      et->welfare = welfare;
      et->welfare_newborn = welfare_newborn;
      et->approval = approval;
      et->approval_newborn = approval_newborn;

      Q = pow( (W/(1.0-alpha-alpha2)) / pow(L,-alpha-alpha2) / pow(K,alpha2) , 1.0/alpha);
      K = ((et->r+delta)/alpha2)*pow(Q,-alpha)*pow(L,alpha2+alpha-1.0);
      K = pow(K,1.0/(alpha2-1.0));
      
      et->Q = Q;
      et->L = L;
      et->K = K;
      et->P = alpha * pow(Q,alpha-nu) * pow(L,1.0-alpha-alpha2) * pow(K,alpha2);
      
      if(benchmark_eqm)
	{
	  ylbar = ylbar*W/et->W;
	  et->ylbar = ylbar;
	  et->G = G;
	  et->lump_sum=0.0;
	}
      else
	{
	  ylbar = et->ylbar;
	  if(taua_clear_gbc || taul_clear_gbc)
	    {
	      G = et->G;
	      lump_sum = 0.0;
	    }
	  else if(public_goods)
	    {
	      et->G=G;
	      lump_sum=0.0;
	    }
	  else if(lump_sum)
	    {
	      et->lump_sum = update_speed*lump_sum + (1.0-update_speed)*et->lump_sum;
	      G = et->G;
	    }
	}
    }
  else
    {
      et->Q=Q;
      et->L=L;
      et->K=K;
    }
  
  et->Y = Y;
  et->wealth_tax_pct=100.0*wealth_tax_mass/total_mass;
  et->tax_shelter_pct=100.0*mass_w_tax_shelter/total_mass;
  et->wealth_tax_revenue=wtaxR;
  et->k_tax_revenue=ktaxR;
  et->w_tax_revenue_lost=wtaxR_lost;
  et->k_tax_revenue_lost=ktaxR_lost;
  et->sumA = sumA;
  et->sumA_rep = sumA_rep;
  et->sumA_hid = sumA_hid;
  et->Kd = K+Kd;
  et->Ks = Ks;

  if(benchmark_eqm)
    {
      *diff = sqrt((diff1*diff1 + 300*diff2*diff2)/2.0);
    }
  else if(taua_clear_gbc || taul_clear_gbc)
    {
      *diff = sqrt((diff1*diff1 + 300*diff2*diff2 + diff3*diff3)/3.0);
    }
  else if (pe==0 && lump_sum)
    {
      *diff = sqrt((diff1*diff1 + 300*diff2*diff2 + 500*diff3*diff3)/3.0);
    }
  else if (pe==2 && lump_sum)
    {
      *diff = sqrt((diff1*diff1 + 500*diff3*diff3)/2.0);
    }
  else
    {
      *diff = sqrt((diff1*diff1 + 300*diff2*diff2)/2.0);
    }
  return 0;
}

/****************************************************************/
/*    Functions related to solving for stationary equilibria    */
/****************************************************************/

int solve_ss(eqm_t * et, int pe)
{
  #ifndef __cplusplus
  #endif
  
  time_t start, stop;
  time(&start);
  int iter=0;
  int status=0;
  double diff = HUGE_VAL;
  double last_diff=diff;

  int flag=4;

  double diff_tol = wage_tol;
  if(pe==1)
    {
      diff_tol = 1.0e-7;
    }
  
  double update_speed = 0.25;
  
  do
    {
      iter++;
      
      //linebreak2();
      if(pe==0 && verbose)
	{
	  if(taua_clear_gbc)
	    {
	      fprintf(logfile,"\n\tIter %d. Price guess: W = %0.4f, r = %0.2e, taua = %0.2e\n",iter,et->W,et->r,taua);
	    }
	  else if(taul_clear_gbc)
	    {
	      fprintf(logfile,"\n\tIter %d. Price guess: W = %0.4f, r = %0.2f, taul = %0.2e\n",iter,et->W,et->r,taul_bar);
	    }
	  else if(lump_sum)
	    {
	      fprintf(logfile,"\n\tIter %d. Price guess: W = %0.4f, r = %0.2e, transfer = %0.2e\n",iter,et->W,et->r,et->lump_sum);
	    }
    	  else
	    {
	      fprintf(logfile,"\n\tIter %d. Price guess: W = %0.4f, r = %0.2e\n",iter,et->W,et->r);
	    }
	}
      else if(pe==2 && verbose)
	{
	  if(lump_sum)
	    {
	      fprintf(logfile,"\n\tIter %d. Price guess: W = %0.4f, transfer = %0.2e\n",iter,et->W,et->lump_sum);
	    }
	  else
	    {
	      fprintf(logfile,"Error!\n");
	      return 1;
	    }
	}
      
      if(benchmark_eqm)
	{
	  initialize_asset_grid(et->W,et->a_grid);
	}
      
      if(solve_ss_bellman(et))
	{
	  status=1;
	  break;
	}

      if(solve_ss_dist(et) && dist_max_iter>1)
	{
	  status=1;
	  break;
	}
      
      if(aggregate(et,update_speed,pe,verbose,1,&diff))
	{
	  status=1;
	  break;
	}

      write_eqm_t_csv(et,csv);
      
      if(pe!=1 && verbose)
	{
	  fprintf(logfile,"\tPrice vector RMSE = %+0.8f\n",diff);
	}

      if(diff>last_diff && (diff<3*diff_tol || (diff<5*diff_tol && last_diff<1.5*diff_tol)) && iter>18)
	{
	  break;
	}

      last_diff=diff;
    }
  while(iter<wage_max_iter && !(fabs(diff)<=diff_tol && flag==4) && status==0);

  time(&stop);
  
  if((0 && iter==wage_max_iter) || status==1)
    {
      fprintf(logfile,"\nFailed to solve for stationary equilibrium! Time = %0.2f.\n",difftime(stop,start));
      return 1;
    }
  else
    {
      fprintf(logfile,"Stationary equilibrium found! Time = %0.2f.\n",difftime(stop,start));

      //fprintf(logfile,"\tStoring equilbrium data\n");

      return 0;
    }
}

void store_bin_eqm_t(eqm_t * et, char const * fname)
{
  //fprintf(logfile,"Storing equilibrium in binary\n");
  
  char * fname2 = concat(fname,".bin");
  FILE * file = fopen(fname2,"wb");

  if(file)
    {  
      fwrite(&(et->L),sizeof(double),1,file);
      fwrite(&(et->Q),sizeof(double),1,file);
      fwrite(&(et->K),sizeof(double),1,file);
      fwrite(&(et->W),sizeof(double),1,file);
      fwrite(&(et->r),sizeof(double),1,file);
      fwrite(&(et->P),sizeof(double),1,file);
      fwrite(&(et->ylbar),sizeof(double),1,file);
      fwrite(&(et->lump_sum),sizeof(double),1,file);
      fwrite(&(et->G),sizeof(double),1,file);
      fwrite(&taua,sizeof(double),1,file);
      fwrite(&tauk,sizeof(double),1,file);
      fwrite(&taul_bar,sizeof(double),1,file);
      fwrite(&abar,sizeof(double),1,file);

      fwrite( (double *)et->a_grid, NA*sizeof(double), 1, file);
      fwrite( (double *)et->ga, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      fwrite( (int *)et->gv, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(int), 1, file);
      fwrite( (double *)et->V, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      fwrite( (double *)et->Psi, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
  
      fclose(file);
    }
  else
    {
      fprintf(logfile,"Couldn't open file to store eqm_t struct!\n");
    }
}

int load_bin_eqm_t(eqm_t * et, char const * fname)
{
  fprintf(logfile,"Loading equilibrium from binary file %s\n",fname);
    
  FILE * file = fopen(fname,"rb");

  if(file)
    {
      int tmp=0;
      tmp += fread(&(et->L),sizeof(double),1,file);
      tmp += fread(&(et->Q),sizeof(double),1,file);
      tmp += fread(&(et->K),sizeof(double),1,file);
      tmp += fread(&(et->W),sizeof(double),1,file);
      tmp += fread(&(et->r),sizeof(double),1,file);
      tmp += fread(&(et->P),sizeof(double),1,file);
      tmp += fread(&(et->ylbar),sizeof(double),1,file);
      tmp += fread(&(et->lump_sum),sizeof(double),1,file);
      tmp += fread(&(et->G),sizeof(double),1,file);
      tmp += fread(&taua,sizeof(double),1,file);
      tmp += fread(&tauk,sizeof(double),1,file);
      tmp += fread(&taul_bar,sizeof(double),1,file);
      tmp += fread(&abar,sizeof(double),1,file);
      
      tmp += fread( (double *)et->a_grid, NA*sizeof(double), 1, file);
      
      tmp += fread( (double *)et->ga, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      tmp += fread( (int *)et->gv, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(int), 1, file);
      tmp += fread( (double *)et->V, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      tmp += fread( (double *)et->Psi, 2*J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      
      fclose(file);
      
      if(tmp!=18)
	{
	  fprintf(logfile,"Failed to read eqm binary data!\n");
	  return 1;
	}
      else
	{
	  return 0;
	}
    }
    else
    {
      fprintf(logfile,"Couldn't open file to load eqm_t struct!\n");
      return 1;
    }
}

void write_eqm_t_csv(const eqm_t * et, char const * fname)
{
  FILE * file = fopen(fname,"w");
  if(file != NULL)
    {
      fprintf(file,"Y,A,A_rep,A_hid,L,Q,K,W,P,r,taua,abar,tauk,taul,WtaxRev,KtaxRev,WtaxRev_lost,KtaxRev_lost,G,lump_sum,ylbar,welfare,welfare0,approval,approval0,wtax_evasion_elast\n");
      fprintf(file,"%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f,%0.8f\n",
	      et->Y,et->sumA/et->Y,et->sumA_rep/et->Y,et->sumA_hid/et->Y,
	      et->L,et->Q,et->K,et->W,et->P,et->r,
	      taua,abar,tauk,taul_bar,
	      100.0*et->wealth_tax_revenue/et->Y,
	      100.0*et->k_tax_revenue/et->Y,
	      100.0*et->w_tax_revenue_lost/et->Y,
	      100.0*et->k_tax_revenue_lost/et->Y,
	      et->G,et->lump_sum,et->ylbar,
	      et->welfare,
	      et->welfare_newborn,
	      et->approval,et->approval_newborn,et->wealth_tax_evasion_elast);
      fclose(file);
    }
  else
    {
      fprintf(logfile,"Failed to open CSV file %s!\n",fname);
      fprintf(logfile,"Fopen error code: %s\n",strerror(errno));
    }
}

int load_benchmark_eqm_t(eqm_t * et)
{
  fprintf(logfile,"Loading data from calibrated benchmark equilibrium baseline models\n");

  FILE * file = fopen("output/ss0_evasion.bin","rb");

  int tmp=0;
  if(file)
    {
      tmp += fread(&(et->L),sizeof(double),1,file);
      tmp += fread(&(et->Q),sizeof(double),1,file);
      tmp += fread(&(et->K),sizeof(double),1,file);
      tmp += fread(&(et->W),sizeof(double),1,file);
      tmp += fread(&(et->r),sizeof(double),1,file);
      tmp += fread(&(et->P),sizeof(double),1,file);
      tmp += fread(&(et->ylbar),sizeof(double),1,file);
      tmp += fread(&(et->lump_sum),sizeof(double),1,file);
      tmp += fread(&(et->G),sizeof(double),1,file);
      tmp += fread(&taua,sizeof(double),1,file);
      tmp += fread(&tauk,sizeof(double),1,file);
      tmp += fread(&taul_bar,sizeof(double),1,file);
      tmp += fread(&abar,sizeof(double),1,file);
      
      tmp += fread( (double *)et->a_grid, NA*sizeof(double), 1, file);
      
      tmp += fread( (double *)(&(et->ga[1])), J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      tmp += fread( (int *)(&(et->gv[1])), J*NE*NZ*NI*NA*NV*ND*sizeof(int), 1, file);
      tmp += fread( (double *)(&(et->V[1])), J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      tmp += fread( (double *)(&(et->Psi[1])), J*NE*NZ*NI*NA*NV*ND*sizeof(double), 1, file);
      
      fclose(file);
    }
  if(file && tmp!=18)
    {
      fprintf(logfile,"Failed to read eqm binary data!\n");
      return 1;
    }  
  else if(!file)
    {
      fprintf(logfile,"Couldn't open file with evasion eqm data!\n");
      return 1;
    }

  for(int j=0; j<J; j++)
    {
      for(int ie=0; ie<NE; ie++)
	{
	  for(int iz=0; iz<NZ; iz++)
	    {
	      for(int ii=0; ii<NI; ii++)
		{
		  for(int ia=0; ia<NA; ia++)
		    {
		      et->ga[0][j][ie][iz][ii][0][0][ia] = et->ga[1][j][ie][iz][ii][0][0][ia];
		      et->gv[0][j][ie][iz][ii][0][0][ia] = 0;
		      et->V[0][j][ie][iz][ii][0][0][ia] = et->V[1][j][ie][iz][ii][0][0][ia];
		      et->Psi[0][j][ie][iz][ii][0][0][ia] = 0.0;
		      for(int iv=0; iv<NV; iv++)
			{
			  for(int id=0; id<ND; id++)
			    {
			      et->Psi[0][j][ie][iz][ii][0][0][ia] += (1.0-et->frac_v)*et->Psi[1][j][ie][iz][ii][iv][id][ia];
			      et->Psi[1][j][ie][iz][ii][iv][id][ia]  = et->Psi[1][j][ie][iz][ii][iv][id][ia] * et->frac_v;
			    }
			}
		      
		    }
		}
	    }
	}
    }

  return 0;
}



/****************************************************************/
/*    Functions related to testing                              */
/****************************************************************/

void calc_wealth_tax_evasion_elast(int verbose_, eqm_t * et)
{
  double taxable_wealth_total_bench=0.0;
  double taxable_wealth_bench[8] = {0.0};
  double taxable_wealth_total_new=0.0;
  double taua_[8] = {taua,taua+0.01,taua+0.02,taua+0.03,taua+0.04,taua+0.05,taua+0.06,taua+0.07};
  double abar_[8] = {abar*32./50.,abar,abar*5.,abar*10.,abar*20.,abar*50.,abar*100.,abar*200.};

  for(int vv=0; vv<2; vv++)
    {
      for(int j=0; j<J; j++)
	{
	  for(int ie=0; ie<NE; ie++)
	    {
	      for(int iz=0; iz<NZ; iz++)
		{
		  for(int ii=0; ii<NI; ii++)
		    {
		      for(int iv=0; iv<(evasion_type==1 && vv==1 ? NV : 1); iv++)
			{
			  for(int id=0; id<(evasion_type==1 && iv>0 ? 2 : 1); id++)
			    {
			      for(int ia=0; ia<NA; ia++)
				{
				  double a = ss0->a_grid[ia];
				  double v = 1.0-evasion_grid[iv];
			  
				  if(wealth_tax_type==0 || wealth_tax_type==5 || wealth_tax_type==6)
				    {
				      taxable_wealth_total_bench += fmax(0.0,a*v-abar) * ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];
				      taxable_wealth_total_new += fmax(0.0,a*v-abar) * et->Psi[vv][j][ie][iz][ii][iv][id][ia];
				    }
				  else if(wealth_tax_type==1)
				    {
				      taxable_wealth_total_bench += fmax(0.0,a*v-abar) * ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];
				      taxable_wealth_total_new += fmax(0.0,a*v-abar) * et->Psi[vv][j][ie][iz][ii][iv][id][ia];
				      taxable_wealth_bench[1] += fmax(0.0,a*v-abar*20) * ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];
				      taxable_wealth_bench[0] += (fmax(0.0,a*v-abar)-fmax(0.0,a*v-abar*20)) * ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];

				    }
				  else if(wealth_tax_type==2)
				    {
				      taxable_wealth_total_bench += fmax(0.0,a*v-abar) * ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];
				      taxable_wealth_total_new += fmax(0.0,a*v-abar) * et->Psi[vv][j][ie][iz][ii][iv][id][ia];

				      for(int ibar=0; ibar<8; ibar++)
					{
					  if(ibar<8-1)
					    {
					      taxable_wealth_bench[ibar] += fmax(0.0, fmin(a*v,abar_[ibar+1])-abar_[ibar]) *
						ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];
					    }
					  else
					    {
					      taxable_wealth_bench[ibar] += fmax(0.0,a*v-abar_[ibar]) * ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];
					    }
					}

				    }
				  else if(wealth_tax_type==3)
				    {
				      taxable_wealth_total_bench += fmax(0.0,fmin(abar2,a*v)-abar) * ss0->Psi[vv][j][ie][iz][ii][iv][id][ia];
				      taxable_wealth_total_new += fmax(0.0,fmin(abar2,a*v)-abar) * et->Psi[vv][j][ie][iz][ii][iv][id][ia];
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
  
  double chg_taxable_wealth = log(taxable_wealth_total_new)-log(taxable_wealth_total_bench);
  double avg_tax=0.0;
  
  if(wealth_tax_type==0)
    {
      avg_tax = taua;
    }
  else if(wealth_tax_type==1)
    {
      avg_tax = (0.02*taxable_wealth_bench[0] + 0.03*taxable_wealth_bench[1])/taxable_wealth_total_bench;
    }
  else if(wealth_tax_type==2)
    {
      for(int ibar=0; ibar<8; ibar++)
	{
	  avg_tax += taua_[ibar]*taxable_wealth_bench[ibar];
	  //fprintf(logfile,"\t%0.8f %0.8f\n",taua_[ibar],taxable_wealth_bench[ibar]);
	}
      avg_tax = avg_tax/taxable_wealth_total_bench;
    }
  else if(wealth_tax_type==3)
    {
      avg_tax = 0.02;
    }
  else if(wealth_tax_type==6)
    {
      avg_tax = 0.015;
    }

  avg_tax = log(1.0-avg_tax)-log(1.0);

  if(verbose_)
    {
      fprintf(logfile,"\nAvg tax change: %0.8f\n",avg_tax);
      fprintf(logfile,"Pct change in taxable wealth: %0.8f\n",chg_taxable_wealth);
      fprintf(logfile,"Total elasticity: %0.8f\n\n",chg_taxable_wealth/avg_tax);
      fprintf(logfile,"Evasion elasticity: %0.8f\n\n",chg_taxable_wealth/avg_tax -  1.19);
    }

  et->wealth_tax_evasion_elast = chg_taxable_wealth/avg_tax;
}

  


#endif

