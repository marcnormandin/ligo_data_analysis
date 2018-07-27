#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

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

double pso_fitness_function(gsl_vector *xVec, void  *inParamsPointer){
	assert(xVec != NULL);
	assert(inParamsPointer != NULL);

	unsigned int validPt;
	//! [Cast fit func params]
	struct fitFuncParams *inParams = (struct fitFuncParams *)inParamsPointer;
	//! [Cast fit func params]
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

		// unused, but should be
		int network_statistic_index;

		coherent_network_statistic(
				splParams->network,
				splParams->f_low,
				splParams->f_high,
				&chirp_time,
				&sky,
				splParams->network_strain,
				splParams->workspace[parallel_get_thread_num()],
				&fitFuncVal,
				&network_statistic_index,
				NULL);
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

pso_ranges_t* pso_ranges_alloc(const char *pso_settings_filename) {
	pso_ranges_t *r = (pso_ranges_t*) malloc( sizeof(pso_ranges_t) );
	if (r == 0) {
		printf("Unable to allocate memory for pso_ranges_t structure. Aborting.\n");
		abort();
	}

	/* Read the settings since we need them from file. */
	settings_file_t *settings_file = settings_file_open(pso_settings_filename);
	if (settings_file == NULL) {
		printf("Error opening the PSO settings file (%s). Aborting.\n", pso_settings_filename);
		abort();
	}

	/* First we need the number of dimensions */
	r->nDim = atoi(settings_file_get_value(settings_file, "search_num_dim"));
	if (r->nDim <= 0) {
		printf("search_num_dim (%d) in file (%s) must be >= 1. Aborting.\n", r->nDim, pso_settings_filename );
		abort();
	}

	r->min = (double*) malloc( r->nDim * sizeof(double) );
	if (r->min == 0) {
		printf("Unable to allocate memory for pso_ranges_t structure. Aborting.\n");
		abort();
	}

	r->max = (double*) malloc( r->nDim * sizeof(double) );
	if (r->max == 0) {
		printf("Unable to allocate memory for pso_ranges_t structure. Aborting.\n");
		abort();
	}

	r->rangeVec = (double*) malloc( r->nDim * sizeof(double) );
	if (r->rangeVec == 0) {
		printf("Unable to allocate memory for pso_ranges_t structure. Aborting.\n");
		abort();
	}

	settings_file_close(settings_file);

	return r;
}

void pso_ranges_free( pso_ranges_t* r) {
	assert( r != NULL );
	assert( r->min != NULL);
	assert( r->max != NULL);

	free( r->min );
	free( r->max );
	free( r->rangeVec );
	free( r );
}

void pso_ranges_init(const char *pso_settings_filename, pso_ranges_t* r) {
	assert( pso_settings_filename );
	assert( r );

	/* Read the settings since we need them from file. */
	settings_file_t *settings_file = settings_file_open(pso_settings_filename);
	if (settings_file == NULL) {
		printf("Error opening the PSO settings file (%s). Aborting.\n", pso_settings_filename);
		abort();
	}

	r->min[0] = atof(settings_file_get_value(settings_file, "search_ra_min"));
	r->max[0] = atof(settings_file_get_value(settings_file, "search_ra_max"));

	r->min[1] = atof(settings_file_get_value(settings_file, "search_dec_min"));
	r->max[1] = atof(settings_file_get_value(settings_file, "search_dec_max"));

	r->min[2] = atof(settings_file_get_value(settings_file, "search_chirp_t_0_min"));
	r->max[2] = atof(settings_file_get_value(settings_file, "search_chirp_t_0_max"));

	r->min[3] = atof(settings_file_get_value(settings_file, "search_chirp_t_1_5_min"));
	r->max[3] = atof(settings_file_get_value(settings_file, "search_chirp_t_1_5_max"));

	int i;
	for (i = 0; i < r->nDim; i++) {
		r->rangeVec[i] = r->max[i] - r->min[i];
	}

	settings_file_close(settings_file);
}

double convert_domain_pso_to_ff( pso_ranges_t *ranges, returnData_t *res, int index ) {
	double x = gsl_vector_get( res->bestLocation, index );
	return ( x*ranges->rangeVec[index] + ranges->min[index] );
}

/* Convert from PSO doman to the FF domain. */
void return_data_to_pso_results( pso_ranges_t *ranges, returnData_t *from, pso_result_t *to) {
	to->ra = convert_domain_pso_to_ff( ranges, from, 0 );
	to->dec = convert_domain_pso_to_ff( ranges, from, 1 );
	to->chirp_t0 = convert_domain_pso_to_ff( ranges, from, 2 );
	to->chirp_t1_5 = convert_domain_pso_to_ff( ranges, from, 3 );

	/* PSO finds minimums but we want the largest network statistic */
	to->snr = -1.0 * from->bestFitVal;

	to->total_iterations = from->totalIterations;
	to->total_func_evals = from->totalFuncEvals;
	to->computation_time_secs = from->computationTimeSecs;
}

int pso_estimate_parameters(const char *pso_settings_filename, pso_fitness_function_parameters_t *splParams, current_result_callback_params_t *callback_params, gslseed_t seed, pso_result_t* result) {
	assert(pso_settings_filename != NULL);
	assert(splParams != NULL);
	assert(result != NULL);

	/* Estimate right-ascension, declination, and chirp times. */
	// !Fixme THIS SHOULD BE READ FROM THE PSO SETTINGS FILE
	unsigned int nDim = 4, lpc;
	/* [0] = RA
	   [1] = Declination
	   [2] = Chirp time 0
	   [3] = Chirp time 1.5
	 */
	/* Shihan said his PSO went from 0 to PI for dec */
	//double rmin[4] = {-M_PI, 	-0.5*M_PI, 	0.0, 		0.0};
	//double rmax[4] = {M_PI, 	0.5*M_PI, 	43.4673, 	1.0840};
	pso_ranges_t *pso_ranges = pso_ranges_alloc( pso_settings_filename );
	pso_ranges_init( pso_settings_filename, pso_ranges );

	/* Error handling off */
	//gsl_error_handler_t *old_handler = gsl_set_error_handler_off ();

	/* Initialize random number generator */
	gsl_rng *rngGen = gsl_rng_alloc(gsl_rng_taus);
	/* Soumya version gsl_rng_set(rngGen,2571971); */
	gsl_rng_set(rngGen, seed);

    /* Allocate fitness function parameter struct.
	 */
	struct fitFuncParams *inParams = ffparam_alloc(nDim);

	/* Load fitness function parameter struct */
	for (lpc = 0; lpc < nDim; lpc++){
		gsl_vector_set(inParams->rmin,lpc,pso_ranges->min[lpc]);
		gsl_vector_set(inParams->rangeVec,lpc,pso_ranges->rangeVec[lpc]);
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
	} else {
		printf("Opened the PSO settings file (%s) for reading...\n", pso_settings_filename);
	}

	struct psoParamStruct psoParams;
	psoParams.popsize = atoi(settings_file_get_value(settings_file, "popsize"));
	psoParams.maxSteps= atoi(settings_file_get_value(settings_file, "maxSteps"));
	psoParams.c1 = atof(settings_file_get_value(settings_file, "c1"));
	psoParams.c2 = atof(settings_file_get_value(settings_file, "c2"));
	psoParams.max_velocity = atof(settings_file_get_value(settings_file, "max_velocity"));
	psoParams.dcLaw_a = atof(settings_file_get_value(settings_file, "dcLaw_a"));
	psoParams.dcLaw_b = atof(settings_file_get_value(settings_file, "dcLaw_b"));
	psoParams.dcLaw_c = psoParams.maxSteps;
	psoParams.dcLaw_d = atof(settings_file_get_value(settings_file, "dcLaw_d"));;
	psoParams.locMinIter = atof(settings_file_get_value(settings_file, "locMinIter"));
	psoParams.locMinStpSz = atof(settings_file_get_value(settings_file, "locMinStpSz"));
	psoParams.rngGen = rngGen;
	psoParams.debugDumpFile = NULL; /*fopen("ptapso_dump.txt","w"); */

	const char *pso_version_p = settings_file_get_value(settings_file, "pso_version");
	char *pso_version;
	pso_version = malloc( sizeof(char) * (strlen(pso_version_p)+1) );
	strcpy(pso_version, pso_version_p);

	// Print the settings as a diagnostic
	printf("Using the following PSO settings:\n");
	settings_file_print(settings_file);

	settings_file_close(settings_file);
	printf("Closed the PSO settings file.\n");

	/* Now call the desired PSO implementation */
	if (strcmp(pso_version, "lbest")==0) {
		lbestpso(nDim, fitfunc, inParams, callback_params, &psoParams, psoResults);
	} else if (strcmp(pso_version, "gbest")==0) {
		gbestpso(nDim, fitfunc, inParams, callback_params, &psoParams, psoResults);
	} else if (strcmp(pso_version, "spso")==0) {
		spso(nDim, fitfunc, inParams, callback_params, &psoParams, psoResults);
	} else {
		fprintf(stderr, "Error. pso_version in the pso settings file must be 'lbest', 'gbest', or 'spso'. Exiting.\n");
		exit(-1);
	}

	free(pso_version);

	return_data_to_pso_results( pso_ranges, psoResults, result );

	/* Free allocated memory */
	pso_ranges_free( pso_ranges );
	ffparam_free(inParams);
	returnData_free(psoResults);
	gsl_rng_free(rngGen);

	return 0;
}
