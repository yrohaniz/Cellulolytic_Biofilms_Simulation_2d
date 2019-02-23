//
// Created by yrohanizadegan on 02/03/18.
//

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "CeBiSi.h"

/* This file reads the simulation parameters from the parameter file
 * 'params' */

#define NUM_CHARS 100 //Initialize the size of string variables used below.

//These functions go through the input string and pass the numeric value of the numerical characters
static unsigned long find_ul(char *p, unsigned long *c) {
    ++*c; //update the value of counter which keeps track of input parameters
    while (*p) { // While there are more characters to process...
        if (isdigit(*p)) { // Upon finding a digit, ...
            unsigned long val = strtoul(p, &p, 10); // Read a number, ...
            return val;
        } else { // Otherwise, move on to the next character.
            p++;
        }
    }
}

static unsigned long long find_ull(char *p, unsigned long *c) {
    ++*c;
    while (*p) {
        if (isdigit(*p)) {
            unsigned long long val = strtoull(p, &p, 10);
            return val;
        } else {
            p++;
        }
    }
}

static double find_fl(char *p, unsigned long *c) {
    ++*c;
    while (*p) {
        if (isdigit(*p)) {
            double val = strtod(p, &p);
            return val;
        } else {
            p++;
        }
    }
}

void params_r() {
    FILE *f_in;
    unsigned long n=0;
    unsigned long count_param=0;
    char str[NUM_CHARS]; // Declaration of an string to arbitrary length

    f_in = fopen("params", "r");

    if (fgets(str, NUM_CHARS, f_in) != NULL) { // Check whether params is opened and have contents
        fgets(str, NUM_CHARS, f_in); // Ignore the first line in params
        fgets(str, NUM_CHARS, f_in); // Read the next line which is number of parameters
        n = find_ul(str, &count_param);
        printf("Total number of parameters: %lu\n", n); // Print to the screen the number of parameters
        fgets(str, NUM_CHARS, f_in);
        NTHRDS = find_ul(str, &count_param); // Pick the value of the number of threads
        printf("# of threads: %lu\n", NTHRDS);
        fgets(str, NUM_CHARS, f_in);
        INIT_T = find_fl(str, &count_param); // Pick the value of the start time
        printf("Initial time: %3.2f\n", INIT_T);
        fgets(str, NUM_CHARS, f_in);
        END_T = find_fl(str, &count_param); // Pick the value of the stop time
        printf("Final time: %3.2f\n", END_T);
        fgets(str, NUM_CHARS, f_in);
        SEED_init = find_ull(str, &count_param); // Pick the value of the initial seed for initializing the PRNG algorithm
        printf("Initial seed is: %llu\n", SEED_init);
        fgets(str, NUM_CHARS, f_in);
        SEED_fin = find_ull(str, &count_param); // Pick the value of the final seed for initializing the PRNG algorithm
        printf("Final seed is: %llu\n", SEED_fin);
        fgets(str, NUM_CHARS, f_in);
        INIT_M = find_fl(str, &count_param); // Pick the value of the initial biomass density
        printf("Initial bmass density: %3.2f\n", INIT_M);
        fgets(str, NUM_CHARS, f_in);
        INIT_C = find_fl(str, &count_param); // Pick the value of the initial carbon concentration
        printf("Initial carb concentration: %3.2f\n", INIT_C);
        fgets(str, NUM_CHARS, f_in);
        INIT_PHI = find_fl(str, &count_param); // Pick the value of the initial attachment function value
        printf("Initial value of attachment function: %f\n", INIT_PHI);
        fgets(str, NUM_CHARS, f_in);
        BMASS_INIT_LOC = find_ul(str, &count_param); // Pick the value of the initial biomass deposit location
        printf("Location of the initial bmass deposit is: %ld\n", BMASS_INIT_LOC);
        fgets(str, NUM_CHARS, f_in);
        DOMAIN_L = find_fl(str, &count_param); // Pick the value of the domain length
        printf("Domain length: %2.1f\n", DOMAIN_L);
        fgets(str, NUM_CHARS, f_in);
        DOMAIN_W = find_fl(str, &count_param); // Pick the value of the domain width
        printf("Domain width: %2.1f\n", DOMAIN_W);
        fgets(str, NUM_CHARS, f_in);
        LOWER_BMASS = 1.0 - 1.0/find_fl(str, &count_param); // Pick the value of the biomass concentration lower limit
        printf("Lower limit of biomass density: %10.9f\n", LOWER_BMASS);
        fgets(str, NUM_CHARS, f_in);
        UPPER_BMASS = 1.0 - 1.0/find_fl(str, &count_param); // Pick the value of the biomass concentration upper limit
        printf("Upper limit of biomass density: %10.9f\n", UPPER_BMASS);
        fgets(str, NUM_CHARS, f_in);
        T_PRTN = find_fl(str, &count_param); // Pick the value of the data recording time increments
        printf("Recording time increments: %3.2f\n", T_PRTN);
        fgets(str, NUM_CHARS, f_in);
        BMASS_ATTACH_RATE = find_fl(str, &count_param); // Pick the value of the biomass density attachment rate
        printf("Biomass density attachment rate: %f\n", BMASS_ATTACH_RATE);
        fgets(str, NUM_CHARS, f_in);
        N_coarse = find_ul(str, &count_param); // Pick the value of the coarsest grid size in x direction
        printf("Coarsest grid size in x direction: %lu\n", N_coarse);
        fgets(str, NUM_CHARS, f_in);
        M_coarse = find_ul(str, &count_param); // Pick the value of the coarsest grid size in y direction
        printf("Coarsest grid size in y direction: %lu\n", M_coarse);
        fgets(str, NUM_CHARS, f_in);
        N_fine = find_ul(str, &count_param); // Pick the value of the finest grid size in x direction
        printf("Finest grid size in x direction: %lu\n", N_fine);
        fgets(str, NUM_CHARS, f_in);
        M_fine = find_ul(str, &count_param); // Pick the value of the finest grid size in y direction
        printf("Finest grid size in y direction: %lu\n", M_fine);
        fscanf(f_in, "%s", str);
        if (feof(f_in) && n == (count_param-1)) { // Check whether the end of file is true and if not exit and print warning to screen
            printf("Parameter reading complete.\n\n");
        } else {
            printf("ERROR in file!CHECK params!\n");
            exit(1);
        }
    }

    fclose(f_in); // Close the params file
}
