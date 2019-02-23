//
// Created by yrohanizadegan on 28/05/18.
//

/*This function intializes the random number arrays used by the Marsaglia-polar generator.
 * Here lcg2 and lcg3 are used to fill up the arrays and the multipliers and the increments
 * used by the linear congruential generators are prepared for leapfrogging in the parallel
 * region.*/

#include "CeBiSi.h"

void parseed(unsigned long long rndnum1, unsigned long long rndnum2, unsigned long long **rnd1,
             unsigned long long **rnd2, unsigned long long *mult2_n, unsigned long long *mult3_n,
             unsigned long long *inc2_n, unsigned long long *inc3_n) {
    unsigned long j;
    *mult2_n = 1; //Initialize to one for iterative multiplication
    *mult3_n = 1;
    *inc2_n = 0; //Initialize to zero for iterative summation
    *inc3_n = 0;
    for (j=0; j<((N*M)/2); j++) { //For theses simulations loop through half the grid size
        rndnum1 = lcg2(rndnum1); //Use lcg2 to generate an rn and update the value of rndnum1
        rndnum2 = lcg3(rndnum2); //Use lcg3 to generate an rn and update the value of rndnum2
        rnd1[j][0] = rndnum1; //Store the generated rn in the arrays
        rnd2[j][0] = rndnum2;
        *inc2_n += INC2*(*mult2_n); //Update the value of increments
        *inc3_n += INC3*(*mult3_n);
        *mult2_n = *mult2_n * MULTIPR2; //Update the value of multipliers
        *mult3_n = *mult3_n * MULTIPR3;
    }
    *inc2_n = *inc2_n % MODUL; //Calculate the final value of increments and multipliers using the leapfrogging method formulas
    *mult2_n = *mult2_n % MODUL;
    *inc3_n = *inc3_n % MODUL;
    *mult3_n = *mult3_n % MODUL;
}
