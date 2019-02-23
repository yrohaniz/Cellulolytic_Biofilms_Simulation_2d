//
// Created by yrohanizadegan on 22/09/17.
//

#include "CeBiSi.h"

//Flux across the cell wall (:=arithmetic average of diffusion values
//at the centers of two adjacent cells times the central difference).
//Here d1 and d2 represent the diffusion function values at the adjacent
//centers and m1 and m2 represent the biomass at these centers.
double flux(double m1, double m2) {
    return (D_M(m1) + D_M(m2)) * 0.5 * (m2 - m1) * inv_dlt_x;
}
