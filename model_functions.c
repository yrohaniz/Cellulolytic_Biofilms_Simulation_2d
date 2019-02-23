//
// Created by yrohanizadegan on 20/09/17.
//

#include <math.h>
#include "CeBiSi.h"

//Biomass-diffusion parameters
const int alfa = 4;
//const int beta = 4;
const double delta = 1.0e-06;

//Carbon uptake parameters
const double gamma_c = 0.4;
const double kappa_c = 0.01;

//Biomass net growth parameters
const double mu = 1.0;
const double kappa_m = 0.01;
const double lambda = 0.42;

//Biomass-diffusion coefficient
double D_M(double m) {
    return delta * (pow(m/(1.0-m), alfa));
}

//Uptake rate of carbon substrate
double G_C(double c) {
    return gamma_c * c / (kappa_c + c);
}

//Net biomass growth rate
double F_C(double c) {
    return (mu * c / (kappa_m + c)) - lambda;
}
