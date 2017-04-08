#include <stdio.h>
/* The header file tells us what fitness function we are calling
and what the parameter structure for this function is.
*/
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>

#include "ptapso.h"
#include "ptapso_maxphase.h"

#include "inspiral_pso_fitness.h"
#include "inspiral_source.h"
#include "inspiral_chirp.h"
#include "random.h"
#include "sky.h"

#include "settings.h"

/*
fitFuncVal = ptapsotestfunc_gsl(xVec,P)
A benchmark test function for PTAPSO
that computes the Rastrigin fitness function for
each row of xVec.  The fitness value is returned in fitFuncVal.
xVec is standardized, that is 0<=xVec(i,j)<=1.
The values used to convert xVec(i,j)
internally before computing fitness are given in P.rmin and
P.rangeVec:
xVec(j) -> xVec(j)*rangevec(j)+rmin(j).
fitFuncVal = infty if the point xVec falls
outside the hypercube defined by 0<=xVec(j)<=1.
The real coordinates are returned in P.realCoord.
Soumya D. Mohanty, Jan 2016
- Derived from ptapsotestfunc.c. Converts to gsl_vector inputs and
  uses the interface needed by GSL multi-dimensional local minimization routines.
*/

double ptapso_func(gsl_vector *xVec, void  *inParamsPointer){

	unsigned int validPt;
    unsigned int lpc;
	//! [Cast fit func params]
	struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;
	//! [Cast fit func params]
	unsigned int ncols = inParams->nDim;
	double rangeVec, rmin, x;
	double fitFuncVal;
	//! [Cast special params]
	/* Shows how to retrieve special parameters (dummy ones are supplied for this particular
	fitness function).
	*/
	struct ptapso_func_params *splParams = (struct ptapso_func_params *)inParams->splParams;

	/* This fitness function knows what fields are given in the special parameters struct */

	s2rvector(xVec,inParams->rmin,inParams->rangeVec,inParams->realCoord);

	validPt = chkstdsrchrng(xVec);

	if (validPt){
		inParams->fitEvalFlag = 1;
		fitFuncVal = 0;

		double ra = gsl_vector_get(inParams->realCoord, 0);
		double dec = gsl_vector_get(inParams->realCoord, 1);
		double chirp_time_0 = gsl_vector_get(inParams->realCoord, 2);
		double chirp_time_1_5 = gsl_vector_get(inParams->realCoord, 3);

		chirp_time_t chirp_time;
		CF_CT_compute(splParams->f_low, chirp_time_0, chirp_time_1_5, &chirp_time);

		/* The network statistic requires the time of arrival to be zero
		   in order for the matched filtering to work correctly. */


		sky_t sky;
		sky.ra = ra;
		sky.dec = dec;

		/* DANGER */
		double polarization_angle = 0.0;

		coherent_network_statistic(
				splParams->network,
				splParams->net_asd,
				splParams->f_low,
				splParams->f_high,
				&chirp_time,
				&sky,
				polarization_angle,
				splParams->signals,
				splParams->workspace,
				&fitFuncVal);
		/* The statistic is larger for better matches, but PSO is finding
		   minimums, so multiply by -1.0. */
		fitFuncVal *= -1.0;
    }
	else{
		fitFuncVal=GSL_POSINF;
		inParams->fitEvalFlag = 0;
	}
   return fitFuncVal;
}

int ptapso_estimate_ra_dec_t0_t1_5(char *pso_settings_filename, ptapso_fun_params_t *splParams, gslseed_t seed, pso_result_t* result) {
	/* Estimate right-ascension and declination */
	unsigned int nDim = 4, lpc;
	/* [0] = RA
	   [1] = Declination
	   [2] = Chirp time 0
	   [3] = Chirp time 1.5
	 */
	/* Shihan said his PSO went from 0 to PI for dec */
	double rmin[4] = {-M_PI, 	-0.5*M_PI, 	0.0, 		0.0};
	double rmax[4] = {M_PI, 	0.5*M_PI, 	43.4673, 	1.0840};
	double rangeVec[4];

	/* Error handling off */
	gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();

	/* Initialize random number generator */
	gsl_rng *rngGen = gsl_rng_alloc(gsl_rng_taus);
	/* Soumya version gsl_rng_set(rngGen,2571971); */
	gsl_rng_set(rngGen, seed);

    /* Allocate fitness function parameter struct.
	 */
	struct fitFuncParams *inParams = ffparam_alloc(nDim);

	/* Load fitness function parameter struct */
	for (lpc = 0; lpc < nDim; lpc++){
		rangeVec[lpc]=rmax[lpc]-rmin[lpc];
		gsl_vector_set(inParams->rmin,lpc,rmin[lpc]);
		gsl_vector_set(inParams->rangeVec,lpc,rangeVec[lpc]);
	}
	/* Set up pointer to fitness function. Use the prototype
	declaration given in the header file for the fitness function. */
	double (*fitfunc)(gsl_vector *, void *) = ptapso_func;
	/* Set up special parameters, if any, needed by the fitness function used.
	   These should be provided in a structure that should be defined in
	   the fitness function's header file.
	 */
	/* ============================================ */
	/* Pass on the special parameters through the generic fitness function parameter
	struct */
	inParams->splParams = splParams;

	/* Set up storage for output from ptapso. */
	struct returnData *psoResults = returnData_alloc(nDim);


	/* nelder-meade method .. look up */
	/* Set up the pso parameter structure.*/

	/* Load the pso settings */
	settings_file_t *settings_file = settings_file_open(pso_settings_filename);
	if (settings_file == NULL) {
		printf("Error opening the PSO settings file (%s). Aborting.\n", pso_settings_filename);
		abort();
	}

	struct psoParamStruct psoParams;
	psoParams.popsize = atoi(settings_file_get_value(settings_file, "popsize"));;
	psoParams.maxSteps= atoi(settings_file_get_value(settings_file, "maxSteps"));;
	psoParams.c1 = atof(settings_file_get_value(settings_file, "c1"));;
	psoParams.c2 = atof(settings_file_get_value(settings_file, "c2"));;
	psoParams.max_velocity = atof(settings_file_get_value(settings_file, "max_velocity"));;
	psoParams.dcLaw_a = atof(settings_file_get_value(settings_file, "dcLaw_a"));;
	psoParams.dcLaw_b = atof(settings_file_get_value(settings_file, "dcLaw_b"));;
	psoParams.dcLaw_c = psoParams.maxSteps;
	psoParams.dcLaw_d = atof(settings_file_get_value(settings_file, "dcLaw_d"));;
	psoParams.locMinIter = atof(settings_file_get_value(settings_file, "locMinIter"));
	psoParams.locMinStpSz = atof(settings_file_get_value(settings_file, "locMinStpSz"));
	psoParams.rngGen = rngGen;
	psoParams.debugDumpFile = NULL; /*fopen("ptapso_dump.txt","w"); */

	settings_file_close(settings_file);

	ptapso(nDim, fitfunc, inParams, &psoParams, psoResults);

	/*fclose(psoParams.debugDumpFile);*/

	/* Information returned by PSO */
	/*
	printf("Total number of iterations %zu\n", psoResults->totalIterations);
	printf("Total number of function evaluations %zu\n", psoResults->totalFuncEvals);
	printf("Best Location found: \n");
	for (lpc = 0; lpc < nDim; lpc++){
		printf("%f, ",gsl_vector_get(psoResults->bestLocation,lpc));
	}
	printf("\n");
	printf("Best Fitness Value: %f\n", psoResults->bestFitVal);
	*/

	/* convert values to function ranges, instead of pso ranges */
	s2rvector(psoResults->bestLocation,inParams->rmin,inParams->rangeVec,inParams->realCoord);
	result->ra = gsl_vector_get(inParams->realCoord, 0);
	result->dec = gsl_vector_get(inParams->realCoord, 1);
	result->chirp_t0 = gsl_vector_get(inParams->realCoord, 2);
	result->chirp_t1_5 = gsl_vector_get(inParams->realCoord, 3);
	result->snr = -1.0 * psoResults->bestFitVal;

	/* Free allocated memory */
	ffparam_free(inParams);
	returnData_free(psoResults);
	gsl_rng_free(rngGen);

	return 0;
}
