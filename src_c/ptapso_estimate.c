#include <stdio.h>
#include "ptapso.h"
/* The header file tells us what fitness function we are calling
and what the parameter structure for this function is.
*/
#include "ptapso_func.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include "ptapso_maxphase.h"
#include "random.h"
#include "ptapso_estimate.h"

int ptapso_estimate(ptapso_fun_params_t *splParams, gslseed_t seed, size_t max_steps, pso_result_t* result) {
	/* Estimate right-ascension and declination */
	unsigned int nDim = 2, lpc;
	/* [0] = RA
	   [1] = Declination */
	double rmin[2] = {-M_PI, -0.5*M_PI};
	double rmax[2] = {M_PI, 0.5*M_PI};
	double rangeVec[2];

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
	struct psoParamStruct psoParams;
	psoParams.popsize=40;
	psoParams.maxSteps= max_steps;
	psoParams.c1=2;
	psoParams.c2=2;
	psoParams.max_velocity = 0.2;
	psoParams.dcLaw_a = 0.9;
	psoParams.dcLaw_b = 0.4;
	psoParams.dcLaw_c = psoParams.maxSteps;
	psoParams.dcLaw_d = 0.2;
	psoParams.locMinIter = 0; /* original value was 10 */
	psoParams.locMinStpSz = 0.01;
	psoParams.rngGen = rngGen;
	psoParams.debugDumpFile = NULL; /*fopen("ptapso_dump.txt","w"); */

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
	result->snr = -1.0 * psoResults->bestFitVal;

	/* Free allocated memory */
	ffparam_free(inParams);
	returnData_free(psoResults);
	gsl_rng_free(rngGen);

	return 0;
}
