#ifndef __CALIBRATION_C__
#define __CALIBRATION_C__

#include "calibration.h"

/****************************************************************/
/*    Some useful utilities                                     */
/****************************************************************/

void discretize_z()
{
  int n = NZ;
  int i,j;
  double mup= 3.0;
  double mdown = 3.0;

  if(constraint_flag)
    mup=5;
  
  double ucond_std = sqrt(sig_z*sig_z/(1.0-rho_z*rho_z));
  double lo = -mdown*ucond_std;
  double hi = mup*ucond_std;
  double d = (hi-lo)/(n-1.0);
  linspace(lo,hi,n,z_grid);
  
  for(i=0; i<n; i++)
    {
      double x = z_grid[i];
      
      for(j=0; j<n; j++)
	{
	  double y = z_grid[j];
	  
	  if(j==0)
	    {
	      z_probs[i][j] = gsl_cdf_ugaussian_P( (y + d/2.0 - rho_z*x) / sig_z);
	    }
	  else if(j==(n-1))
	    {
	      z_probs[i][j] = 1.0 - gsl_cdf_ugaussian_P( (y - d/2.0 - rho_z*x) / sig_z);
	    }
	  else
	    {
	      z_probs[i][j] = (gsl_cdf_ugaussian_P( (y + d/2.0 - rho_z*x) / sig_z) -
			  gsl_cdf_ugaussian_P( (y - d/2.0 - rho_z*x) / sig_z));
	    }
	}
    }  
  
  for(i=0; i<n; i++)
    {
      z_grid[i] = exp(z_grid[i]);
      if(constraint_flag && i==(n-1))
	{
	  if(evasion_type==0)
	    {
	      z_grid[i] = 22.218455752773007;
	    }
	  else
	    {
	      z_grid[i] = 23.347457528048778;
	    }
	     
	}
      
      double sum = 0.0;
      
      for(j=0; j<n; j++)
	{
	  sum += z_probs[i][j];
	}
      if(fabs(sum-1.0)>1.0e-10)
	{
	  fprintf(logfile,"warning: transition probabilities do not sum to 1 for state %d!\n",i);
	}
    }

  double diff=HUGE_VAL;
  SET_ALL_V(z_ergodic_dist,NZ,1.0/((int)(NZ)));
  do
    {
      double tmp[NZ]={0.0};
      diff=0.0;
      for(i=0; i<NZ; i++)
	{
	  for(j=0; j<NZ; j++)
	    {
	      tmp[i] += z_probs[j][i]*z_ergodic_dist[j];
	    }
	  if(fabs(tmp[i]-z_ergodic_dist[i])>diff)
	    {
	      diff=fabs(tmp[i]-z_ergodic_dist[i]);
	    }
	}      
      memcpy(z_ergodic_dist,tmp,sizeof(double)*NZ);
    }
  while(diff>1.0e-12);

}

void discretize_e()
{
  int n = NE;
  int i,j;
  double mup = 3.0;
  double mdown = 3.0;
  double ucond_std = sig1_e/sqrt(1.0-rho1_e*rho1_e);
  double lo = -mdown*ucond_std;
  double hi = mup*ucond_std;
  double d = (hi-lo)/(n-1.0);
  linspace(lo,hi,n,e_grid);
  
  for(i=0; i<n; i++)
    {
      double x = e_grid[i];
      
      for(j=0; j<n; j++)
	{
	  double y = e_grid[j];
	  
	  if(j==0)
	    {
	      e_probs[i][j] = gsl_cdf_ugaussian_P( (y + d/2.0 - rho1_e*x) / sig1_e);
	    }
	  else if(j==(n-1))
	    {
	      e_probs[i][j] = 1.0 - gsl_cdf_ugaussian_P( (y - d/2.0 - rho1_e*x) / sig1_e);
	    }
	  else
	    {
	      e_probs[i][j] = (gsl_cdf_ugaussian_P( (y + d/2.0 - rho1_e*x) / sig1_e) -
			  gsl_cdf_ugaussian_P( (y - d/2.0 - rho1_e*x) / sig1_e));
	    }

	  if(j==0)
	    {
	      e_birth_probs[i][j] = gsl_cdf_ugaussian_P( (y + d/2.0 - rho2_e*x) / sig2_e);
	    }
	  else if(j==(n-1))
	    {
	      e_birth_probs[i][j] = 1.0 - gsl_cdf_ugaussian_P( (y - d/2.0 - rho2_e*x) / sig2_e);
	    }
	  else
	    {
	      e_birth_probs[i][j] = (gsl_cdf_ugaussian_P( (y + d/2.0 - rho2_e*x) / sig2_e) -
			  gsl_cdf_ugaussian_P( (y - d/2.0 - rho2_e*x) / sig2_e));
	    }

	}
    }  
  
  for(i=0; i<n; i++)
    {
      e_grid[i] = exp(e_grid[i]);
	    
      double sum = 0.0;
      
      for(j=0; j<n; j++)
	{
	  sum += e_probs[i][j];
	}
      if(fabs(sum-1.0)>1.0e-10)
	{
	  fprintf(logfile,"warning: transition probabilities do not sum to 1 for state %d!\n",i);
	}
    }

    for(i=0; i<n; i++)
    {
      double sum = 0.0;
      
      for(j=0; j<n; j++)
	{
	  sum += e_birth_probs[i][j];
	}
      if(fabs(sum-1.0)>1.0e-10)
	{
	  fprintf(logfile,"warning: transition probabilities do not sum to 1 for state %d!\n",i);
	}
    }


    double diff=HUGE_VAL;
    SET_ALL_V(e_ergodic_dist,NE,0.0);
    double tmp1[J][NE]={{0.0}};
    SET_ALL_V(tmp1[0],NE,1.0/((int)(NE)));
    do
      {
	double tmp2[J][NE]={{0.0}};
	diff=0.0;
	for(int age=0; age<J; age++)
	  {
	    for(i=0; i<NE; i++)
	      {
		for(j=0; j<NE; j++)
		  {
		    if(age<(J-1))
		      {
			if(age<(R-1))
			  {
			    tmp2[age+1][i] += e_probs[j][i]*tmp1[age][j]*phi[age];
			  }
			else if(i==j)
			  {
			    tmp2[age+1][i] += tmp1[age][i]*phi[age];
			  }
		      }
		    tmp2[0][i] += e_birth_probs[j][i]*tmp1[age][j]*(1.0-phi[age]);
		  }
	      }
	  }
	for(int age=0; age<J; age++)
	  {
	    for(i=0; i<NE; i++)
	      {
		if(fabs(tmp1[age][i]-tmp2[age][i])>diff)
		  {
		    diff=fabs(tmp1[age][i]-tmp2[age][i]);
		  }
	      }
	  }
	memcpy(tmp1,tmp2,sizeof(double)*J*NE);
      }
    while(diff>1.0e-12);

    double sum=0.0;
    for(int age=0; age<J; age++)
      {
	for(i=0; i<NE; i++)
	  {
	    e_ergodic_dist[i]+=tmp1[age][i];
	    sum+=tmp1[age][i];
	  }
      }
    for(i=0; i<NE; i++)
      {
	e_ergodic_dist[i] = e_ergodic_dist[i]/sum;
      }
}

/****************************************************************/
/*    Parameter values (calibrated elsewhere)                   */
/****************************************************************/

void set_params(int evasion_type_)
{
  fprintf(logfile,"Setting and checking parameter values...\n");

  evasion_type=evasion_type_;
  
  // preferences
  sigma = 4.0;

  if(evasion_type==0)
    {
      if(constraint_flag)
	{
	  beta = 0.95792734672921165;
	  gama = 0.44268781887422776;
	}
      else	
	{
	  // flat taul
	  //beta = 0.96931893488594656;
	  //gama = 0.43504254396315795;

	  // prog taul, labor share = 0.6
	  //beta = 0.97367360653240609;
	  //gama = 0.42913860042472707;

	  // prog taul, labor share = 0.58
	  beta = 0.96960738778840472;
	  gama = 0.43230748778721129;
	}
    }
  else if(evasion_type==1)
    {
      if(constraint_flag)
	{
	  beta = 0.95665963228378059;
	  gama = 0.44406038159396188;
	}
      else
	{
	  // flat taul
	  //beta = 0.97514911477187793;
	  //gama = 0.43187713127031535;

	  // prog taul, labor share = 0.6
	  //beta = 0.97962200628012586;
	  //gama = 0.42659410138603765;

	  // prog taul, labor share = 0.58
	  beta=0.97549773509010362;
	  gama = 0.42908235664034916;
	}
    }
  
  // OLG stuff 
  //double az = 0.5/30.;

  double phi_[] = {0.999064143747091000, 0.999057437351439000, 0.999047139252070000,
		   0.999028513848316000, 0.999002128606662000, 0.998970763757824000,
		   0.998937236843630000, 0.998900839709677000, 0.998862511711195000,
		   0.998819889500737000, 0.998765300260856000, 0.998697769246064000,
		   0.998622532351874000, 0.998539241380058000, 0.998442850308492000,
		   0.998336586169898000, 0.998207374825142000, 0.998037528945133000,
		   0.997822653735056000, 0.997576926602050000, 0.997324218973517000,
		   0.997068772325292000, 0.996794620528817000, 0.996494883671402000,
		   0.996169835096225000, 0.995823254343122000, 0.995465271640568000,
		   0.995096642058342000, 0.994716361630707000, 0.994315713178366000,
		   0.993882585782557000, 0.993410578928887000, 0.992905301041901000,
		   0.992374176625162000, 0.991819550283253000, 0.991233213804662000,
		   0.990602982230484000, 0.989915066398680000, 0.989137337543070000,
		   0.988242183811962000, 0.987190308049321000, 0.985988594591617000,
		   0.984710191376507000, 0.983398867771029000, 0.981995487585663000,
		   0.980451723560690000, 0.978705510497093000, 0.976724572479724000,
		   0.974471978843212000, 0.971939291805028000, 0.969180349260568000,
		   0.966224979609251000, 0.962747562676668000, 0.958864040672779000,
		   0.954589251428842000, 0.949853785336017000, 0.944555081427097000,
		   0.938727524131536000, 0.932235985994338000, 0.924182154238224000,
		   0.000000000000000000};

  SET_ALL_V(zeta,J,0.0);
  int j;
  for(j=0; j<J; j++)
    {
      phi[j] = phi_[j];
      if(j<R)
	{
	  zeta[j] = -(pow((double)j,2.0)/1800.0 -(double)j/30.0);
	  zeta[j] = exp(zeta[j]);
	}
    }

  // employment states
  rho1_e = 0.93712872266769409;
  sig1_e = 0.20193409919738770;
  rho2_e = 0.18484359979629517;
  sig2_e = 0.56884044408798218;
  discretize_e();

  e_ergodic_dist_cum[0] = e_ergodic_dist[0];
  int ie;
  for(ie=0; ie<NE; ie++)
    {
      if(ie>0)
	{
	  e_ergodic_dist_cum[ie] = e_ergodic_dist_cum[ie-1]+e_ergodic_dist[ie];
	}

      e_probs_cum[ie][0] = e_probs[ie][0];
      e_birth_probs_cum[ie][0] = e_birth_probs[ie][0];
      int iep;
      for(iep=1; iep<NE; iep++)
	{
	  e_probs_cum[ie][iep] = e_probs_cum[ie][iep-1]+e_probs[ie][iep];
	  e_birth_probs_cum[ie][iep] = e_birth_probs_cum[ie][iep-1]+e_probs[ie][iep];
	}

    }

  // entrepreneurial ability states
  rho_z = 0.1;

  if(evasion_type==0)
    {
      if(constraint_flag)
	{
	  sig_z = 0.38238900899887085;
	}
      else
	{
	  //sig_z = 0.40469038996252732; // flat taul
	  //sig_z = 0.40417448376203285; // prog taul, labor share = 0.6
	  sig_z = 0.38238901391392599;
	}
    }
  else if(evasion_type==1)
    {
      if(constraint_flag)
	{
	  sig_z = 0.38238900899887085;
	}
      else
	{
	  //sig_z = 0.44210022798784226; // flat taul
	  //sig_z = 0.44254610748972162; // prog taul, labor share = 0.6
	  sig_z = 0.41706789919573989; // prog taul, labor share = 0.58	  
	}
    }
  
  discretize_z();
  pi=8.7399996817111969E-002;
  double pi_entry=2.2590890526771545E-002;
  double pi_exit = 8.1000000238418579E-002;
  double i_probs_[NI][NI] = {{1.0-pi_exit, pi_exit},
			     {pi_entry,1.0-pi_entry}};

  int ii, iip;
  for(ii=0; ii<NI; ii++)
    {
      for(iip=0; iip<NI; iip++)
	{
	  i_probs[ii][iip] = i_probs_[ii][iip];
	  if(iip==0)
	    {
	      i_probs_cum[ii][iip]=i_probs[ii][iip];
	    }
	  else
	    {
	      i_probs_cum[ii][iip]=i_probs[ii][iip]+i_probs_cum[ii][iip-1];
	    }
	}
    }
  
  z_ergodic_dist_cum[0] = z_ergodic_dist[0];
  int iz;
  for(iz=0; iz<NZ; iz++)
    {
      if(iz>0)
	{
	  z_ergodic_dist_cum[iz] = z_ergodic_dist_cum[iz-1]+z_ergodic_dist[iz];
	}

      z_probs_cum[iz][0] = z_probs[iz][0];
      int izp;
      for(izp=1; izp<NZ; izp++)
	{
	  z_probs_cum[iz][izp] = z_probs_cum[iz][izp-1]+z_probs[iz][izp];
	}

    }

  // financial markets
  delta = 5.0000000745058060e-2;
 
  // production
  alpha2 = 7.1000002324581146E-002;
  //alpha = 0.4-alpha2;
  alpha = 0.42 - alpha2;
  nu = 0.89999997615814209;
  
  // taxes
  tauc = 7.5000002980232239E-002;
  taua = 0.0;
  abar = 0.0;
  kbar = 0.0;
  tauk2 = 0.0;

  lambda=0;
  lambda2=0;
  if(evasion_type==0)
    {
      tauk = 0.25;
      if(constraint_flag)
	{
	  lambda = +HUGE_VAL;
	}
      else
	{
	  //lambda = 6*0.35191770664228800; // flat taul
	  //lambda = 0.35128012791982371*6; // prog taul, labor share = 0.6
	  lambda = 0.34042467824755857*6;
	}
    }
  else if(evasion_type==1)
    {
      tauk=0.25;

      if(constraint_flag)
	{
	  lambda=+HUGE_VAL;
	}
      else
	{
	  //lambda =  2.1157284445420297; // flat taul
	  //lambda = 0.34624538481515987*6; // prog taul, labor share = 0.6
	  lambda = 0.33416905664610919*6; // prog taul, labor share = 0.58
	}
    }
  
  if(taul_flag==0)
    {
      taul[0] = 0.0;
      taul[1] = 0.0;
      taul[2] = 0.125;
      taul[3] = 0.244;
      taul[4] = 0.272;
    }
  else
    {
      for(int ie=0; ie<NE; ie++)
	{
	  taul[ie] = 0.224;
	}
    }
  
  abar=0.0;
  abar2=0.0;

  // social security... compute average efficiency units contingent on ending in a particular state
  double ybar[NE];
  for(ie=0; ie<NE; ie++)
    {
      double tmp[NE] = {0.0};
      
      tmp[ie]=1.0;
      ybar[ie]=0.0;

      for(j=R-1; j>=0; j--)
	{
	  int ie2;
	  double tmp2[NE] = {0.0};
	  for(ie2=0; ie2<NE; ie2++)
	    {
	      ybar[ie] += tmp[ie2] * zeta[j] * e_grid[ie2];

	      double sum = 0.0;
	      int ie4;
	      for(ie4=0; ie4<NE; ie4++)
		{
		  sum += e_probs[ie4][ie2];
		}
	      
	      int ie3;
	      for(ie3=0; ie3<NE; ie3++)
		{
		  tmp2[ie3] += tmp[ie2]*e_probs[ie3][ie2]/sum;
		}
	    }
	  for(ie2=0; ie2<NE; ie2++)
	    {
	      tmp[ie2]=tmp2[ie2];
	    }
	}
      ybar[ie]=ybar[ie]/R;
    }

  double mean_ybar=0.0;
  for(ie=0; ie<NE; ie++)
    {
      mean_ybar+=ybar[ie]*e_ergodic_dist[ie];
    }

  for(ie=0; ie<NE; ie++)
    {
      double tmp = ybar[ie]/mean_ybar;
      if(tmp <= 0.3)
	{
	  Phi[ie] = 0.9*tmp;
	}
      else if(tmp<=2.0)
	{
	  Phi[ie] = 0.27+0.32*(tmp-0.3);
	}
      else if(tmp<=4.1)
	{
	  Phi[ie] = 0.81+0.15*(tmp-2.0);
	}
      else
	{
	  Phi[ie]=1.13;
	}
     
    }
  
  theta=0.0;
  eta=0.0;
  eta2=0.0;
  chi=1.0;
  p1n=0;
  p1e=0;
  p2n=0;
  p2e=0;
  penalty_frac_stock=0.0;
  penalty_frac_tauk=0;
  penalty_frac_taua=0;
  max_penalty_years=0;
  
  if(evasion_type==1)
    {
      if(constraint_flag)
	{
	  theta = 1.0327821855903083;
	  eta = 0.22619538477411605;
	  chi=1.0;
	  eta2=1.0;
	  p2n=2.0000000949949026E-003/1000;
	  p2e=0.34717764805068441/1000;
	}
      else
	{
	  // flat taul
	  //theta=1.1156962066771572;
	  //eta=0.10017364664422068;
	  //chi=1.0;
	  //eta2=1.0;
	  //p2n=2.0000000949949026E-003/1000;
	  //p2e=0.34717764805068441/1000;

	  // prog taul, labor share = 0.6
	  //theta = 1.1084151318436195;
	  //eta = 0.10349749162752282;
	  //p2n = 2.0000000949949026E-003/1000;
	  //p2e = 0.33704013546096434/1000;
	  //chi=1.0;
	  //eta2=1.0;
	  
	  // prog taul, labor share = 0.58
	  theta = 1.1486886780691947;
	  eta = 0.10830660893029871;
	  p2n = 0.008/1000;
	  p2e = 0.24441908216101330/1000;
	  chi=1.0;
	  eta2=1.0;

	}
      penalty_frac_stock = 0.5;
      penalty_frac_tauk = 1.75;
      penalty_frac_taua = 1.75;
      max_penalty_years=3;
    }

  if(evasion_type>0 && NV==12)
    {
      double evasion_grid_[NV] = {0.0000000000000000,
				  1.9999999552965164E-002,
				  3.9999999105930328E-002,
				  5.9999998658895493E-002,
				  7.9999998211860657E-002,
				  0.10000000149011612,
				  0.25000000000000000,
				  0.40000000596046448,
				  0.55000001192092896,
				  0.69999998807907104,
				  0.85000002384185791,
				  1.0000000000000000};
      //double evasion_grid_[NV] = {0.0};
      for(int iv=0; iv<NV; iv++)
	{
	  evasion_grid[iv]=evasion_grid_[iv];
	}
    }
  else
    {
      evasion_grid[0]=0.0;
    }
  
  asset_grid_exp = 5.0;
  asset_grid_ub_mult = 60000;
  root_tol_abs = 1.0e-8;
  root_tol_rel = 1.0e-8;
  max_root_iter=100;
  fmin_ftol_abs = 1.0e-11;
  fmin_ftol_rel = 1.0e-11;
  fmin_xtol_abs = 1.0e-6;
  fmin_xtol_rel = 1.0e-6;
  fmin_max_iter = 1000;
  vf_tol_abs = 5.0e-8;
  vf_tol_rel = 5.0e-8;
  vf_max_iter = 25;
  bound_mult = 2.0;
  dist_max_iter = 2000;
  dist_tol = 1.0e-11;
  wage_max_iter=50;
  wage_tol = 1.0e-4;
  max_trans_iter=75;
  verbose = 0;
  exploit_monotonicity=0;
  
  double quantiles_[NQUANTILES] = {0.9,0.95,0.99,0.995,0.999,0.9999};
  int iq;
  for(iq=0; iq<NQUANTILES; iq++)
    {
      quantiles[iq] = quantiles_[iq];
    }
  
  return;

}

double interp(gsl_interp_accel * acc, const double *xa, const double *ya, int n, double x)
{
  double x0=0.0;
  double x1=0.0;
  double xd=0.0;
  double q0=0.0;
  double q1=0.0;
  double retval=0.0;

  int ix = gsl_interp_accel_find(acc, xa, n, x);

  if(ix==0)
    {
      x0 = xa[0];
      x1 = xa[1];
      xd = x1-x0;
      q0 = ya[0];
      q1 = ya[1];
    }
  else if(ix==n-1)
    {
      x0 = xa[n-2];
      x1 = xa[n-1];
      xd = x1-x0;
      q0 = ya[n-2];
      q1 = ya[n-1];
    }
  else
    {
      x0 = xa[ix];
      x1 = xa[ix+1];
      xd = x1-x0;
      q0 = ya[ix];
      q1 = ya[ix+1];
    }

  retval = ( q0*(x1-x) + q1*(x-x0) ) / xd;
  return retval;
}

#endif
