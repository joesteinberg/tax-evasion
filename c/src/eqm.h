#ifndef __EQM_H__
#define __EQM_H__

#include "calibration.h"

/****************************************************************/
/*    Object declarations                                       */
/****************************************************************/

typedef struct
{  
  // aggregates
  double L; // labor supply
  double Q; // intermediate aggregate
  double K; // corp capital
  double Y; // GDP;
  double W; // wage
  double r; // interest rate
  double P; // = alpha * Q^{alpha-nu} * L^{1-alpha}
  double ylbar; // average labor income
  double lump_sum;
  double G; // government consumption = total revenues - SS payments
  double C; // total C
  double wealth_tax_pct;
  double tax_shelter_pct;
  double wealth_tax_revenue;
  double w_tax_revenue_lost;
  double k_tax_revenue;
  double k_tax_revenue_lost;
  double gini;
  double p90_share;
  double p99_share;
  double p999_share;
  double p9999_share;
  double gini_total;
  double p90_share_total;
  double p99_share_total;
  double p999_share_total;
  double p9999_share_total;
  double sumA;
  double sumA_rep;
  double sumA_hid;
  double Kd;
  double Ks;
  double welfare;
  double welfare_newborn;
  double approval;
  double approval_newborn;
  double detection_rate;
  double detection_revenue;
  double wealth_tax_evasion_elast;
  
  double a_grid[NA];
  
  double gyk[NZ][NI][NV][NA]; // capital income policy
  double gk[NZ][NI][NV][NA]; // capital demand policy
  double gks[NZ][NI][NV][NA]; // capital supply policy
  double gtk[NZ][NI][NV][NA]; // capital income tax
  double gta[NV][NA]; // wealth tax
  double gta_evaded[NV][NA]; // wealth tax evaded
  
  double ga[J][NE][NZ][NI][NV][ND][NA]; // final saving policy
  int gv[J][NE][NZ][NI][NV][ND][NA]; // evasion policy
  //double gtk_evaded[J][NE][NZ][NI][NV][ND][NA]; // detection probability
  int giap0[J][NE][NZ][NI][NV][ND][NA]; // saving policy
  double V[J][NE][NZ][NI][NV][ND][NA]; // value function
  
  // distribution stuff
  double Psi[J][NE][NZ][NI][NV][ND][NA];
  double reported_wealth_dist[NQUANTILES];
  double total_wealth_dist[NQUANTILES];
  double shelter_pct_by_wealth_dist[NQUANTILES+1];
  double conceal_pct_by_wealth_dist[NQUANTILES+1];
  double shelter_pct_by_income_dist[NQUANTILES+1];
  double evasion_pct_by_income_dist[NQUANTILES+1];
  double conceal_pct_by_income_dist[NQUANTILES+1];

  double shelter_pct_e;
  double shelter_pct_n;
  double conceal_share_e;
  double conceal_share_n;
  double evasion_share_e;
  double evasion_share_n;
  double conceal_frac_e;
  double conceal_frac_n;
  double evasion_frac_e;
  double evasion_frac_n;  
  
}eqm_t;

typedef struct
{
  int j;
  int ie;
  int iz;
  int ii;
  double W;
  double gross_inc_nol;
  double business_inc;
  double evasion_flow;
  double inc_nol_post_saving_evasion;
  double pd;
  double c;
  double l;
  gsl_spline ******* spline;
  gsl_interp_accel * acc;

}bellman_params;

#ifndef EXTERN2
#define EXTERN2 extern
#endif

//double EW[J][NE][NZ][NI][NV][ND][NA];
EXTERN2 gsl_spline ******* spline_V[NTH];
EXTERN2 gsl_interp_accel * acc[NTH];

EXTERN2 double tmp_Psi[J][NE][NZ][NI][NV][ND][NA];
EXTERN2 double income_cdf_w_shelter_pct[J*NE*NZ*NI*NV*ND*NA][6];

EXTERN2 eqm_t * ss0;
EXTERN2 eqm_t * ss1;
EXTERN2 eqm_t * eqm_trans[NT];

//double penalty_frac;

#define ivd_flag 1

/****************************************************************/
/*    Function declarations                                     */
/****************************************************************/

void alloc_mem();
void free_mem();
void alloc_trans_mem();
void free_trans_mem();

void store_bin_eqm_t(eqm_t * et, char const * fname);
int load_bin_eqm_t(eqm_t * et, char const * fname);
void store_txt_eqm_t(eqm_t * et, char const * fname);
int load_txt_eqm_t(eqm_t * et, char const * fname);
void copy_eqm_t(eqm_t * dest, const eqm_t * src);

void initialize_asset_grid(double W, double * a_grid);
void init_eqm_t(eqm_t * et);

void set_wealth_tax_params(const eqm_t * et, int wealth_tax_type);
void set_capital_policies(eqm_t * et);
int iterate_bellman(eqm_t * et, eqm_t * etp, double * supnorm, int fine_opt);
int solve_ss_bellman(eqm_t * et);

void init_dist(eqm_t * et);
void update_dist(eqm_t * et, double * supnorm, int check_conv);
int solve_ss_dist(eqm_t * et);

double excess_k_demand(double r, void * params);
int clear_k_market(eqm_t * et, double * r_new);
int aggregate(eqm_t * et, double update_speed, int pe, int verbose, int clear_k_mkt_flag, double * diff);
int solve_ss(eqm_t * et, int pe);

int solve_trans();

void lf_dyn(eqm_t * et);
void write_asset_grid(eqm_t * et);
void write_vf_dist(eqm_t * et, char const * prefix);
void write_eqm_t_csv(const eqm_t * et, char const * fname);

void disp_trans();
void write_trans_csv(char const * fname);

void print_wealth_grid(eqm_t * et);
void test_saving_policy(eqm_t * et);
void test_capital_policy(eqm_t * et);
void calc_wealth_tax_evasion_elast(int verbose_, eqm_t * et);
int write_reg_data(eqm_t * et);

#endif
