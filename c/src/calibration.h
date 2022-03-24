#ifndef __CALIBRATION_H__
#define __CALIBRATION_H__

#include <unistd.h>
#include <getopt.h>
#include <ctype.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <limits.h>
#include <memory.h>
#include <errno.h>
#include <omp.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_roots.h>


/****************************************************************/
/*    Declaration of constants and parameters                   */
/****************************************************************/

int sens_type;
int benchmark_eqm; // flag indicating we are solving benchmark, and hence need to store some additional stuff
int taua_clear_gbc;
int taul_clear_gbc;
int wealth_tax_type;
int public_goods;
int lump_sum;
int constraint_flag;
int detection_revenue_flag;
int taul_flag;

// preferences
double sigma; // risk aversion/IES
double beta; // discount factor
double gama; // consumption share in utility
double sigma2; // frisch

// OLG stuff 
#define J 61 // lifespan
#define R 41 // retirement age
double phi[J]; // survival probabilities
double zeta[J]; // deterministic income growth

// employment states
#define NE 5 // number of employment shocks
double rho1_e; // LF persistence
double sig1_e; // LF dispersion
double rho2_e; // LF persistence
double sig2_e; // LF dispersion
double e_grid[NE]; // state grid
double e_ergodic_dist[NE]; // ergodic distribution
double e_ergodic_dist_cum[NE]; // ergodic distribution
double e_birth_probs[NE][NE]; // transition matrix
double e_birth_probs_cum[NE][NE]; // transition matrix
double e_probs[NE][NE]; // transition matrix
double e_probs_cum[NE][NE]; // transition matrix

// entrepreneurial ability states
#define NZ 7 // number of persistent ability states
#define NI 2 // number of entrepreneurial shocks
double rho_z; // IG persistence
double sig_z; // IG dispersion
double z_grid[NZ]; // persistent state grid
double z_ergodic_dist[NZ]; // ergodic distribution
double z_ergodic_dist_cum[NZ]; // ergodic distribution
double z_probs[NZ][NZ]; // inheritance transition matrix
double z_probs_cum[NZ][NZ]; // inheritance transition matrix
double pi;
double i_probs[NI][NI]; // shock transition probabilities
double i_probs_cum[NI][NI]; // shock transition probabilities

// financial markets
#define NA 200
double lambda;
double lambda2;
double delta;

// production
double alpha; // capital share
double alpha2; // corporate capital share
double nu; // elasticity of substitution

// taxes
double tauc; // consumption tax
double tauk; // capital income tax
double tauk2; // capital income tax
double g; // gov't spending
double taul_bar;
double taul[NE]; // labor income tax
double Phi[NE]; // SS benefits
double M50_threshold;
double taua;
double abar;
double abar2;
double kbar;

// evasion
int evasion_type;
double theta;
double eta;
double eta2;
double chi;
double p1n;
double p1e;
double p2n;
double p2e;
double penalty_frac_stock;
double penalty_frac_tauk;
double penalty_frac_taua;
int max_penalty_years;
#define NV 1
#define ND 1
double evasion_grid[NV];

// computational stuff
#ifdef _OPENMP
#define NTH 20
#else
#define NTH 1
#endif
double asset_grid_ub_mult;
double asset_grid_exp;
double root_tol_abs;
double root_tol_rel;
int max_root_iter;
double fmin_xtol_abs;
double fmin_xtol_rel;
double fmin_ftol_abs;
double fmin_ftol_rel;
int fmin_max_iter;
double vf_tol_abs;
double vf_tol_rel;
int vf_max_iter;
double bound_mult;
double dist_tol;
int dist_max_iter;
int wage_max_iter;
double wage_tol;
int max_trans_iter;
#define FINE_GRID_SIZE 51
#define NQUANTILES 6
double quantiles[NQUANTILES];
int verbose;
int exploit_monotonicity;
char fname0[128];
char pref1[128];
char csv[128];
int write_binary;
FILE * logfile;
#define NT 100

/****************************************************************/
/*    Declaration of calibration-related functions              */
/****************************************************************/
void set_params(int evasion_type_);

/****************************************************************/
/*    Simple helper functions                                   */
/****************************************************************/

static inline void open_log(int lognum)
{
  char tmp[12];
  sprintf(tmp,"log%u.txt",lognum);
  logfile = fopen(tmp,"wb");
}

static inline void close_log()
{
  if(logfile != NULL && logfile != stdout)
    {
      fclose(logfile);
    }
}

static inline char* concat(const char *s1, const char *s2)
{
  char *result = (char *)malloc(strlen(s1)+strlen(s2)+1);
  strcpy(result, s1);
  strcat(result, s2);
  return result;
}

static inline void linebreak()
{
  printf("\n////////////////////////////////////////////////////////////////////////////\n\n");
}

static inline void linebreak2()
{
  printf("\n----------------------------------------------------------------------------\n");
}

static inline void i_1d_to_3d(int nx, int ny, int nz, int i1d, int * x, int * y, int * z)
{
  *z = i1d/(nx*ny);
  i1d -= ((*z)*nx*ny);
  *y = i1d/nx;
  *x = i1d%nx;
}

static inline double util1(double c)
{
  return pow(c,gama*(1.0-sigma))/(1.0-sigma);
}

static inline double util2(double c, double l)
{
  return pow(c,gama*(1.0-sigma))*pow((1.0-l),(1.0-gama)*(1.0-sigma))/(1.0-sigma);
}

static inline double util3(double c, double l)
{
  return pow(c,1.0-sigma)/(1.0-sigma) - gama*pow(l,1.0+sigma2)/(1.0+sigma2);
}

static inline double prob_detect(int type, double x)
{
  if(evasion_type==0)
    {
      return 0.0;
    }
  //else if(type == 0 && x>1.0e-11)
  else if(type == 0)
    {
      //return 1.0/(1.0+p1e*exp(p2e*fmax(x,0.0)));
      return tanh(p2e*fmax(x,0.0));
    }
  //else if(type == 1 && x>1.0e-11)
  else if(type == 1)
    {
      //return 1.0/(1.0+p1n/x);
      //return 1.0/(1.0+p1n*exp(p2n*fmax(x,0.0)));
      return tanh(p2n*fmax(x,0.0));
    }  
  else
    {
      return 0.0;
    }
}

static inline void linspace(double lo, double hi, int n, double * v)
{
  double d=(hi-lo)/(n-1.0);
  v[0]=lo;
  int i=0;
  for(i=1;i<n;i++)
    {
      v[i] = v[i-1]+d;
    }
}

static inline void expspace(double lo, double hi, int n, double ex, double * v)
{
  linspace(0.0,pow(hi-lo,1.0/ex),n,v);
  int i;
  for(i=0;i<n;i++)
    {
      v[i] = pow(v[i],ex)+lo;
    }
  return;
}

static inline void reverse(double * v, int n)
{
  double * tmp = (double *)malloc(n*sizeof(double));
  memcpy(tmp,v,n*sizeof(double));
  int i;
  for(i=0; i<n; i++)
    {
      v[i]=tmp[n-1-i];
    }
  free(tmp);
}

static inline void set_all_v(double * v, int n, double x)
{
  int i;
  for(i=0; i<n; i++)
    {
      v[i] = x;
    }
}
#define SET_ALL_V(v,n,x) set_all_v( (double *)(v), (n), (x) )

static inline double sum(const double * v, int n)
{
  int i;
  double mysum = 0.0;
  for(i=0; i<n; i++)
    {
      mysum = mysum + v[i];
    }
  return mysum;
}
#define SUM(v,n) (sum ( (double *)(v), (n) ))

static inline double maxv(const double * v, int n)
{
  int i;
  double max = -HUGE_VAL;
  for(i=0; i<n; i++)
    {
      max = fmax(max,v[i]);
    }
  return max;
}
#define MAXV(v,n) (maxv ( (double *)(v), (n) ))

static inline double dot_prod(const double * v1, const double * v2, int n)
{
  int i;
  double sum = 0.0;
  for(i=0; i<n; i++)
    {
      sum = sum + v1[i]*v2[i];
    }
  return sum;
}
#define DOT_PROD(v1,v2,n) (dot_prod( (double *)(v1), (double *)(v2), (n) ))

// linear interpolation
double interp(gsl_interp_accel * acc, const double *xa, const double *ya, int n, double x);


#endif
