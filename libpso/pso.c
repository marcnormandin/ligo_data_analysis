/*
 * pso.c
 *
 *  Created on: May 2, 2017
 *      Author: marcnormandin
 */

#include <stdio.h>
#include <stddef.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>

#include "ptapso_maxphase.h"
#include "parallel.h"

#include "pso.h"

/*Dummy function to wrap supplied function so that the
  call to the GSL local optimizer is compatible */
double dummyfitfunc(const gsl_vector *xVec, void *dffParams){
	struct dummyFitFuncParam *dfp = (struct dummyFitFuncParam *)dffParams;
	double (*trueFitFunc)(gsl_vector *, void *) = dfp->trufuncPr;
	gsl_vector *xVec2 = (gsl_vector *)xVec;
	double funcVal = trueFitFunc(xVec2,dfp->trufuncParam);

	/* added by Marc */
	return funcVal;
}

/*! Initializer of particle position, velocity, and other properties. */
void initPsoParticles(struct particleInfo *p, size_t nDim, gsl_rng *rngGen){

	double rngNum;
	size_t lpCoord;

	particleinfo_alloc(p,nDim);

	for (lpCoord = 0; lpCoord < nDim; lpCoord++){
		rngNum = gsl_rng_uniform(rngGen);
		gsl_vector_set(p->partCoord,lpCoord,rngNum);
	}

	for (lpCoord = 0; lpCoord < nDim; lpCoord++){
		rngNum = gsl_rng_uniform(rngGen);
		rngNum = - gsl_vector_get(p->partCoord,lpCoord)
			     + rngNum;
		gsl_vector_set(p->partVel,lpCoord,rngNum);
	}

	if(gsl_vector_memcpy(p->partPbest, p->partCoord))
		printf("Error in copying vectors\n");

	p->partSnrPbest = GSL_POSINF;
	p->partSnrCurr = 0;
	p->partSnrLbest = GSL_POSINF;
	p->partInertia = 0;
	p->partFitEvals = 0;
}


/*! Allocate storage for members of particleInfo structure */
void particleinfo_alloc(struct particleInfo *p, size_t nDim){
	p->partCoord = gsl_vector_alloc(nDim); /* Current coordinates */
	p->partVel = gsl_vector_alloc(nDim);  /* Current velocity */
	p->partPbest = gsl_vector_alloc(nDim); /* Coordinates of pbest */
	p->partLocalBest = gsl_vector_alloc(nDim); /* Coordinates of neighborhood best */
}

/*! Free the storage assigned to members of particleInfo structure */
void particleinfo_free(struct particleInfo *p){
	gsl_vector_free(p->partCoord);
	gsl_vector_free(p->partVel);
	gsl_vector_free(p->partPbest);
	gsl_vector_free(p->partLocalBest);
}

/*! Allocate storage for returnData struct members */
struct returnData * returnData_alloc(size_t nDim){
	struct returnData *psoResults = (struct returnData *)malloc(sizeof(struct returnData));
	gsl_vector *bestLocation = gsl_vector_alloc(nDim);
	psoResults->bestLocation = bestLocation;
	return psoResults;
}

/*! Free storage assigned to returnData struct members */
void returnData_free(struct returnData *psoResults){
	gsl_vector_free(psoResults->bestLocation);
	free(psoResults);
}


/*! Dump information stored in particleInfo struct */
void particleinfo_fwrite(FILE *outF, struct particleInfo *p){

	size_t nDim = p->partCoord->size;
	size_t lpc;

	fprintf(outF,"Particle locations in standardized coordinates\n");
	for (lpc = 0; lpc < nDim; lpc++){
		fprintf(outF,"%f ",gsl_vector_get(p->partCoord,lpc));
	}
	fprintf(outF,"\n");
	fprintf(outF,"Particle velocities in standardized coordinates\n");
	for (lpc = 0; lpc < nDim; lpc++){
		fprintf(outF,"%f ",gsl_vector_get(p->partVel,lpc));
	}
	fprintf(outF,"\n -------- \n");
}

/*! Dump particleInfo struct array information as a matrix
with all information pertaining to one particle in a row.
*/
void particleInfoDump(FILE *outF, struct particleInfo *p, size_t popsize){
	if (outF == NULL) {
		return;
	}

	size_t nDim = p[0].partCoord->size;
	size_t lpParticles, lpCoord;

	for (lpParticles = 0; lpParticles < popsize; lpParticles++){
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partCoord,lpCoord));
		}
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partVel,lpCoord));
		}
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partPbest,lpCoord));
		}
		fprintf(outF,"%lf ",p[lpParticles].partSnrPbest);
		fprintf(outF,"%lf ",p[lpParticles].partSnrCurr);
		fprintf(outF,"%lf ",p[lpParticles].partSnrLbest);
		fprintf(outF,"%lf ",p[lpParticles].partInertia);
		for(lpCoord = 0; lpCoord < nDim; lpCoord++){
			fprintf(outF,"%lf ",gsl_vector_get(p[lpParticles].partLocalBest,lpCoord));
		}
		fprintf(outF,"X ");
		fprintf(outF,"%zu ",p[lpParticles].partFitEvals);
		fprintf(outF,"\n");
	}
}

