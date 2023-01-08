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

#define J 61 // lifespan
#define R 41 // retirement age
#define NE 5 // number of employment shocks
#define NZ 7 // number of persistent ability states
#define NI 2 // number of entrepreneurial shocks
#define NA 200
#define NV 12
#define ND 2

#ifdef _OPENMP
#define NTH 61
#else
#define NTH 1
#endif

#define FINE_GRID_SIZE 51
#define NQUANTILES 6
#define NT 100

/****************************************************************/
/*    Declaration of constants and parameters                   */
/****************************************************************/

#ifndef EXTERN
#define EXTERN extern
#endif

EXTERN int sens_type;
EXTERN int benchmark_eqm; // flag indicating we are solving benchmark, and hence need to store some additional stuff
EXTERN int taua_clear_gbc;
EXTERN int taul_clear_gbc;
EXTERN int wealth_tax_type;
EXTERN int public_goods;
EXTERN int lump_sum;
EXTERN int constraint_flag;
EXTERN int detection_revenue_flag;
EXTERN int taul_flag;

// preferences
EXTERN double sigma; // risk aversion/IES
EXTERN double beta; // discount factor
EXTERN double gama; // consumption share in utility
EXTERN double sigma2; // frisch

// OLG stuff 
EXTERN double phi[J]; // survival probabilities
EXTERN double zeta[J]; // deterministic income growth

// employment states
EXTERN double rho1_e; // LF persistence
EXTERN double sig1_e; // LF dispersion
EXTERN double rho2_e; // LF persistence
EXTERN double sig2_e; // LF dispersion
EXTERN double e_grid[NE]; // state grid
EXTERN double e_ergodic_dist[NE]; // ergodic distribution
EXTERN double e_ergodic_dist_cum[NE]; // ergodic distribution
EXTERN double e_birth_probs[NE][NE]; // transition matrix
EXTERN double e_birth_probs_cum[NE][NE]; // transition matrix
EXTERN double e_probs[NE][NE]; // transition matrix
EXTERN double e_probs_cum[NE][NE]; // transition matrix

// entrepreneurial ability states
EXTERN double rho_z; // IG persistence
EXTERN double sig_z; // IG dispersion
EXTERN double z_grid[NZ]; // persistent state grid
EXTERN double z_ergodic_dist[NZ]; // ergodic distribution
EXTERN double z_ergodic_dist_cum[NZ]; // ergodic distribution
EXTERN double z_probs[NZ][NZ]; // inheritance transition matrix
EXTERN double z_probs_cum[NZ][NZ]; // inheritance transition matrix
EXTERN double pi;
EXTERN double i_probs[NI][NI]; // shock transition probabilities
EXTERN double i_probs_cum[NI][NI]; // shock transition probabilities

// financial markets
EXTERN double lambda;
EXTERN double lambda2;
EXTERN double delta;

// production
EXTERN double alpha; // capital share
EXTERN double alpha2; // corporate capital share
EXTERN double nu; // elasticity of substitution

// taxes
EXTERN double tauc; // consumption tax
EXTERN double tauk; // capital income tax
EXTERN double tauk2; // capital income tax
EXTERN double g; // gov't spending
EXTERN double taul_bar;
EXTERN double taul[NE]; // labor income tax
EXTERN double Phi[NE]; // SS benefits
EXTERN double M50_threshold;
EXTERN double taua;
EXTERN double abar;
EXTERN double abar2;
EXTERN double kbar;

// evasion
EXTERN int evasion_type;
EXTERN double theta;
EXTERN double eta;
EXTERN double eta2;
EXTERN double chi;
EXTERN double p1n;
EXTERN double p1e;
EXTERN double p2n;
EXTERN double p2e;
EXTERN double penalty_frac_stock;
EXTERN double penalty_frac_tauk;
EXTERN double penalty_frac_taua;
EXTERN int max_penalty_years;
EXTERN double evasion_grid[NV];

// computational stuff
EXTERN double asset_grid_ub_mult;
EXTERN double asset_grid_exp;
EXTERN double root_tol_abs;
EXTERN double root_tol_rel;
EXTERN int max_root_iter;
EXTERN double fmin_xtol_abs;
EXTERN double fmin_xtol_rel;
EXTERN double fmin_ftol_abs;
EXTERN double fmin_ftol_rel;
EXTERN int fmin_max_iter;
EXTERN double vf_tol_abs;
EXTERN double vf_tol_rel;
EXTERN int vf_max_iter;
EXTERN double bound_mult;
EXTERN double dist_tol;
EXTERN int dist_max_iter;
EXTERN int wage_max_iter;
EXTERN double wage_tol;
EXTERN int max_trans_iter;
EXTERN double quantiles[NQUANTILES];
EXTERN int verbose;
EXTERN int exploit_monotonicity;
EXTERN char fname0[128];
EXTERN char pref1[128];
EXTERN char csv[128];
EXTERN int write_binary;
EXTERN FILE * logfile;

/****************************************************************/
/*    Simple helper functions                                   */
/****************************************************************/

void set_params(int evasion_type_);

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

#endif
