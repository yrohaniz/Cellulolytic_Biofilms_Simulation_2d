//
// Created by yrohanizadegan on 20/09/17.
//

#ifndef CELLULOLYTIC_BIOFILMS_SIMULATION_CEBISI_H
#define CELLULOLYTIC_BIOFILMS_SIMULATION_CEBISI_H

void params_r();
extern unsigned long NTHRDS;
extern double INIT_T, END_T;
extern unsigned long long SEED_init;
extern unsigned long long SEED_fin;
extern double INIT_M, INIT_C;
extern double INIT_PHI;
extern unsigned long BMASS_INIT_LOC;
extern double DOMAIN_L, DOMAIN_W;
extern double LOWER_BMASS;
extern double UPPER_BMASS;
extern double T_PRTN;
extern double BMASS_ATTACH_RATE;
extern unsigned long N_coarse, M_coarse;
extern unsigned long N_fine, M_fine;
extern unsigned long N, M;
extern const int LT, UP, RT, DN;
extern const unsigned long long MODUL;
extern const unsigned long long MULTIPR2;
extern const unsigned long long INC2;
extern const unsigned long long MULTIPR3;
extern const unsigned long long INC3;
extern long t_num_invs;
extern double step;
extern double dlt_x, dlt_x_2;
extern double inv_dlt_x;
extern double sqrtstep;
unsigned long long lcg1(unsigned long long sd);
unsigned long long lcg2(unsigned long long sd);
unsigned long long lcg3(unsigned long long sd);
void parseed(unsigned long long rndnum1, unsigned long long rndnum2, unsigned long long **rnd1,
             unsigned long long **rnd2, unsigned long long *mult2_n, unsigned long long *mult3_n,
             unsigned long long *inc2_n, unsigned long long *inc3_n);
double D_M(double m);
double G_C(double c);
double F_C(double c);
double Euler_phi(double y_k_o, double m, double c, double g, double xi);
double Euler_biomass(double y_k_o, double flx, double c, double fi);
double Euler_carb(double y_k_o, double m);
double flux(double m1, double m2);
double *PRNG_uniform(unsigned long long sd);
void MP_PRNG(double *Z1, double *Z2, unsigned long long a, unsigned long long b);

#endif //CELLULOLYTIC_BIOFILMS_SIMULATION_CEBISI_H
