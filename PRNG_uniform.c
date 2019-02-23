//
// Created by yrohanizadegan on 27/11/17.
//

#include <stdlib.h>
#include <stdio.h>
#include "CeBiSi.h"

/* This function uses the linear congruential generator and the Lagged-Fibonacci
 * generator (LFG) algorithms to produce random numbers (U(0,1)) in the unit
 * interval [0,1] where the period of the prng is determined by values of
 * 's' and 'r' defined in the LFG algorithm. The reason to use LFG is for achieving
 * greater period of random number generation. The seed must be such that at least
 * one of the initialization values (indices 0 to r-1 in random_num array) is odd.*/

static int r = 127; // Declaring and initializing r and s values used in LFG
static int s = 97;

double *PRNG_uniform(unsigned long long sd) {
    double *uni_rand; // Declaration of an array which holds long random numbers
    int i;
    unsigned long long *random_num, tmp; // Declaration of an array which holds uniform PRNGs and a temporary variable

    if (s > r) {
        fprintf(stderr, "s must be less than %d (0<s<r)\n", r);
        exit(1);
    } // Error raised if an 's' value greater than 'r' is passed to the function
    if ((uni_rand = (double *) malloc((r-s)*sizeof(double))) == NULL) {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    } // Initialize the returning array as a double
    if ((random_num = (unsigned long long *) malloc(r*sizeof(unsigned long long))) == NULL) {
        fprintf(stderr,"malloc failed\n");
        exit(1);
    } // Initialize random_num to generate integer (>0) random numbers

    tmp = sd; // Initializing tmp to the seed value input through either BM_PRNG or MP_PRNG
    for (i=0; i<r; i++) {
        random_num[i] = tmp; // Assigning value in tmp to the array component 'i'
        tmp = lcg1(random_num[i]);// Use lcg1 to generate the next PRN
    }

    //Use LFG algorithm to generate random numbers and convert them to U(0,1)
    for (i=0; i<(r-s); i++) {
        uni_rand[i] = ((random_num[i] + random_num[i+(r-s)]) % MODUL) / (double)MODUL;
    }
    free(random_num);
    return uni_rand;
}
