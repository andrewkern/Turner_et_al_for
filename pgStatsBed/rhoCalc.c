#include <stdlib.h>
#include <stdio.h>
#include "rhoCalc.h"

#define CHR_2L 1
#define CHR_2R 2
#define CHR_3L 3
#define CHR_3R 4
#define CHR_X  5

double calc_r(int chromosome, int position, int males_dont_recombine, 
              int windowsize)
{
	double location = (double)(position) / 1000000.0;

	double a = 0.0;
	double b = 0.0;
	double c = 0.0;

	double lo_cutoff = 0.0;
	double hi_cutoff = 0.0;

	double r;

	switch(chromosome)
	{
		case CHR_2L:
			a = -0.04197;
			b = 0.65386;
			c = 1.69634;
			lo_cutoff = 300000;
			hi_cutoff = 17840000;
			break;
		case CHR_2R:
			a = -0.022218;
			b = 0.66925;
			c = -1.278363;
			lo_cutoff = 2050000;
			//lo_cutoff = 4231522;
			hi_cutoff = 20660000;
			break;
		case CHR_3L:
			a = -0.019554;
			b = 0.244726;
			c = 2.638370;
			lo_cutoff = 860000;
			hi_cutoff = 19450000;
			break;
		case CHR_3R:
			a = -0.012588;
			b = 0.502048;
			c = -1.781510;
			lo_cutoff = 4530000;
			//lo_cutoff = 8256340;
			hi_cutoff = 25740000;
			break;
		case CHR_X:
			a = -0.029259;
			b = 0.602636;
			c = 1.132521;
			lo_cutoff = 1530000;
			hi_cutoff = 22320000;
			break;
		default:
			fprintf(stderr, "ERROR: unknown chromosome %d\n", chromosome);
			exit(1);	
	}

	if(position < lo_cutoff || position > hi_cutoff)
		return 0;
	
	r = ((a * location * location) + (b * location) + c) * (windowsize / 1000000.0);

	if(males_dont_recombine)
	{
		if(chromosome == CHR_X)
			r *= 2.0/3.0;
		else
			r *= 0.5;	
	}

	if(r < 0.0)
		return 0.0;

	return r;
}

double calc_rho(int chromosome, int position, int males_dont_recombine, 
                int windowsize, int populationsize)
{
	double r = calc_r(chromosome, position, males_dont_recombine, windowsize);
	
	if(chromosome == CHR_X)
		return 3.0 * populationsize * r;
	
	return 4.0 * populationsize * r;
}

