//
// Created by yousefrohan on 12/28/17.
//

/* This file contains three linear congruential algorithms for pseudo random
 * number generation. Each algorithm has its unique multiplier and increment
 * values whereas the modulus is the same for all and also throughout the
 * entire program.*/

static const unsigned long long MULTIPR1 = 9219741426499971445; //Multiplier in lcg1 (L'Ecuyer table 4)
static const unsigned long long INC1 = 493; //Increment in lcg1 (L'Ecuyer table 4) has to be odd
const unsigned long long MODUL = 9223372036854775808; //Modulus in all LCGs (global variable)(2^63)
const unsigned long long MULTIPR2 = 2806196910506780709; //Multiplier in lcg2 (L'Ecuyer table 4)
const unsigned long long INC2 = 3407; //Increment in lcg2 (L'Ecuyer table 4) has to be odd
const unsigned long long MULTIPR3 = 3249286849523012805; //Multiplier in lcg3 (L'Ecuyer table 4)
const unsigned long long INC3 = 1; //Increment in lcg3 (L'Ecuyer table 4) has to be odd

unsigned long long lcg1(unsigned long long sd) {
    return (MULTIPR1 * sd + INC1) % MODUL;
}

unsigned long long lcg2(unsigned long long sd) {
    return (MULTIPR2 * sd + INC2) % MODUL;
}

unsigned long long lcg3(unsigned long long sd) {
    return (MULTIPR3 * sd + INC3) % MODUL;
}
