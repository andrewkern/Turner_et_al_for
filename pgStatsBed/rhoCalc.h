#ifndef __RHO_CALC__
#define __RHO_CALC__

#define CHR_2L 1
#define CHR_2R 2
#define CHR_3L 3
#define CHR_3R 4
#define CHR_X  5

double calc_r(int chromosome, int position, int males_dont_recombine, 
              int windowsize);

double calc_rho(int chromosome, int position, int males_dont_recombine, 
                int windowsize, int populationsize);

#endif

