/*******************************************************************************
* 
* RndGen.h: A header for random generator function => r in [0,1)
* - implementing additive lagged Fibonacci generator
* - Monte Carlo methods for statistical physics: Newman & Barkema
*
* This program is made freely available with the understanding that every 
* copy of this file must include this header and that it comes without any 
* WITHOUT ANY WARRANTY.
*
*******************************************************************************/

#ifndef RNDGEN_H
#define RNDGEN_H

/** for 64-bit integers: #define rng_conv 5.421010862427522e-20 **/
/** for 32-bit integers: #define rng_conv 2.3283064365387e-10   **/

#define BIT64

// #define drandom (rng_conv*(rng_ia[rng_p=rng_mod[rng_p]] += rng_ia[rng_pp=rng_mod[rng_pp]]))

#define rng_a 2416
#define rng_c 374441
#define rng_m 1771875
#define rng_conv1 2423.96743336861

#ifdef BIT64
#define rng_conv2 5.421010862427522e-20
#else
#define rng_conv2 2.3283064365387e-10
#endif

unsigned long rng_ia[55];
int rng_p,rng_pp;
int rng_mod[55];

void seed(unsigned long i) {

    /** Use that to seed the additive generator.  Also setup the mod array **/

    for (int n=0; n<55; n++) {
        rng_ia[n] = (unsigned long)(rng_conv1*(double)(i=(rng_a*i+rng_c)%rng_m));
        rng_mod[n] = n-1;
    }
    rng_mod[0] = 54;

    rng_p = 0;
    rng_pp = 24;

    /** Run off ten thousand random numbers, just to get things started **/

    for (int n=0; n<10000; n++) rng_ia[rng_p=rng_mod[rng_p]] += rng_ia[rng_pp=rng_mod[rng_pp]];
}

double drandom() {
    return rng_conv2*(double)(rng_ia[rng_p=rng_mod[rng_p]] += rng_ia[rng_pp=rng_mod[rng_pp]]);
}

#endif
