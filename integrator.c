//
// Created by yrohanizadegan on 22/09/17.
//
#include <math.h>
#include "CeBiSi.h"
#include <stdio.h>

double sigma = 0.00001;

//The impulse function
double impulse(double x) {
    return 1/(exp((x-1.00)/sigma)+1.0) - 1/(exp((x-0.99)/sigma)+1.0);
}

//The recipe for Euler-Maruyama method
double Euler_phi(double y_k_o, double m, double c, double g, double xi) {
    double dw_t = g * sqrtstep; // Wiener increment, attachment function
    return y_k_o + xi*(1.0 - m)*c*dw_t;
}

double Euler_biomass(double y_k_o, double flx, double c, double fi) {
    //Euler-Maruyama approximation for biomass
    return y_k_o + ((flx*inv_dlt_x + F_C(c)*y_k_o) + (BMASS_ATTACH_RATE)*((1.0 - y_k_o)*c)*impulse(fi))*step;
}

double Euler_carb(double y_k_o, double m) {
    return y_k_o - G_C(y_k_o)*m*step; //Euler approximation for carbon
}
