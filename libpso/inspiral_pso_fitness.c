#include <assert.h>
#include <stdio.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>

#include "pso.h"
#include "ptapso_maxphase.h"

#include "inspiral_pso_fitness.h"
#include "inspiral_chirp.h"
#include "inspiral_chirp_time.h"
#include "random.h"
#include "sky.h"

#include "settings_file.h"

#include "parallel.h"

pso_fitness_function_parameters_t* pso_fitness_function_parameters_alloc(
		double f_low, double f_high, detector_network_t* network, network_strain_half_fft_t *network_strain)
{
	assert(network != NULL);
	assert(network_strain != NULL);

	size_t i;
	pso_fitness_function_parameters_t *params = (pso_fitness_function_parameters_t*) malloc( sizeof(pso_fitness_function_parameters_t) );
	if (params == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for the pso_fitness_function_parameters_t. Exiting.\n");
		exit(-1);
	}

	params->workspace = (coherent_network_workspace_t**) malloc( parallel_get_max_threads() * sizeof(coherent_network_workspace_t*) );
	if (params->workspace == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for params->workspace. Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < parallel_get_max_threads(); i++) {
		params->workspace[i] = CN_workspace_alloc(
				network_strain->num_time_samples, network, network->detector[0]->asd->len,
				f_low, f_high);
	}

	/* Setup the parameter structure for the pso fitness function */
	params->f_low = f_low;
	params->f_high = f_high;
	params->network = network;
	params->network_strain = network_strain;

	fprintf(stderr, "Number of threads: %lu\n", parallel_get_max_threads());

	return params;
}

void pso_fitness_function_parameters_free(pso_fitness_function_parameters_t *params) {
	assert(params != NULL);

	size_t i;
	assert(params != NULL);
	for (i = 0; i < parallel_get_max_threads(); i++) {
		CN_workspace_free(params->workspace[i]);
	}

	assert(params->workspace);
	free(params->workspace);
	params->workspace = NULL;

	free(params);
}

/* this routine was written for the PSO code. */
void CN_template_chirp_time(double f_low, double chirp_time0, double chirp_time1_5, inspiral_chirp_time_t *ct) {
	assert(ct != NULL);

	/* f_low and these are used to compute the required chirp times. */
	ct->chirp_time0 = chirp_time0;
	ct->chirp_time1_5 = chirp_time1_5;

	double calculated_reduced_mass = Chirp_Calc_CalculatedReducedMass(f_low, chirp_time0, chirp_time1_5);
	double calculated_total_mass = Chirp_Calc_CalculatedTotalMass(f_low, chirp_time0, chirp_time1_5);
	double multi_fac_cal = Chirp_Calc_MultiFacCal(f_low, calculated_total_mass);
	double s_mass_ratio_cal = Chirp_Calc_SMassRatioCal(calculated_reduced_mass, calculated_total_mass);

	ct->chirp_time1 =  Chirp_Calc_Time1(f_low, multi_fac_cal, s_mass_ratio_cal);
	ct->chirp_time2 =  Chirp_Calc_Time2(f_low, multi_fac_cal, s_mass_ratio_cal);

	double calc_tchirp = Chirp_Calc_TChirp(ct->chirp_time0, ct->chirp_time1, ct->chirp_time1_5, ct->chirp_time2);

	/* careful that this doesn't use time of arrival because the network statistic wants it as 0 */
	ct->tc = calc_tchirp;
}

double pso_fitness_function(gsl_vector *xVec, void  *inParamsPointer){
	assert(xVec != NULL);
	assert(inParamsPointer != NULL);

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
	struct pso_fitness_function_parameters_s *splParams = (struct pso_fitness_function_parameters_s *)inParams->splParams;

	/* This fitness function knows what fields are given in the special parameters struct */

	gsl_vector *realCoord = inParams->realCoord[parallel_get_thread_num()];

	s2rvector(xVec,inParams->rmin,inParams->rangeVec,realCoord);

	validPt = chkstdsrchrng(xVec);

	if (validPt){
		inParams->fitEvalFlag[parallel_get_thread_num()] = 1;
		fitFuncVal = 0;

		double ra = gsl_vector_get(realCoord, 0);
		double dec = gsl_vector_get(realCoord, 1);
		double chirp_time_0 = gsl_vector_get(realCoord, 2);
		double chirp_time_1_5 = gsl_vector_get(realCoord, 3);

		inspiral_chirp_time_t chirp_time;
		CN_template_chirp_time(splParams->f_low, chirp_time_0, chirp_time_1_5, &chirp_time);

		/* The network statistic requires the time of arrival to be zero
		   in order for the matched filtering to work correctly. */


		sky_t sky;
		sky.ra = ra;
		sky.dec = dec;

		coherent_network_statistic(
				splParams->network,
				splParams->f_low,
				splParams->f_high,
				&chirp_time,
				&sky,
				splParams->network_strain,
				splParams->workspace[parallel_get_thread_num()],
				&fitFuncVal);
		/* The statistic is larger for better matches, but PSO is finding
		   minimums, so multiply by -1.0. */
		fitFuncVal *= -1.0;
    }
	else{
		fitFuncVal=GSL_POSINF;
		inParams->fitEvalFlag[parallel_get_thread_num()] = 0;
	}

   return fitFuncVal;
}

int pso_estimate_parameters(char *pso_settings_filename, pso_fitness_function_parameters_t *splParams, gslseed_t seed, pso_result_t* result) {
	assert(pso_settings_filename != NULL);
	assert(splParams != NULL);
	assert(result != NULL);

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
	double (*fitfunc)(gsl_vector *, void *) = pso_fitness_function;
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

	gbestpso(nDim, fitfunc, inParams, &psoParams, psoResults);

	/* convert values to function ranges, instead of pso ranges */
	// use the 0 index to convert the value
	s2rvector(psoResults->bestLocation,inParams->rmin,inParams->rangeVec,inParams->realCoord[0]);
	result->ra = gsl_vector_get(inParams->realCoord[0], 0);
	result->dec = gsl_vector_get(inParams->realCoord[0], 1);
	result->chirp_t0 = gsl_vector_get(inParams->realCoord[0], 2);
	result->chirp_t1_5 = gsl_vector_get(inParams->realCoord[0], 3);

	/* PSO finds minimums but we want the largest network statistic */
	result->snr = -1.0 * psoResults->bestFitVal;

	/* Free allocated memory */
	ffparam_free(inParams);
	returnData_free(psoResults);
	gsl_rng_free(rngGen);

	return 0;
}
