//
//  spso.cpp
//  spso
//
//  Created by Ethan Fahnestock on 7/18/18.
//  Some code has been taken from Marc Normandy
//
// NOTES
// - This algorithm assigns the particle inertias, but does not force the c1 value (used as C in SPSO)
//     just checks that the value matches the expected value
// - Current boundary condition is reflecting walls as per the paper. But I think we should add free flying BCs
//
// UNSURE
//
// RECCOMENDED VALUES (FROM PAPER)
// Number of particles: 40
// Number of neighbors (excluding self): 3
// W = 1/(2ln(2))
// C (uses c1 in psoParams) = 0.5 + ln(2)

#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stddef.h>

#include "pso.h" //for initializing variables and
#include "ptapso_maxphase.h" //for fitfunc struct
#include "parallel.h" //for parallel additions

#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_matrix.h> //used to store neighborhoods

void spso(size_t nDim, /*!< Number of search dimensions */
          fitness_function_ptr fitfunc, /*!< Pointer to Fitness function */
          void *ffParams, /*!< Fitness function parameter structure */
          current_result_callback_params_t *callback_params, /* Pointer to callback function parameter structure */
          struct psoParamStruct *psoParams, /*!< PSO parameter structure */
          struct returnData *psoResults /*!< Output structure */) {
    
    clock_t time_start = clock();
    
    
    /*
     Random numbers can be generated on the fly or read from a file. If a
     valid file is specified, generation of random numbers is overriden.
     */
    gsl_rng *rngGen = psoParams->rngGen;
    
    /* Number of particles */
    const size_t popsize = psoParams->popsize;
    /* Number of iterations */
    const size_t maxSteps = psoParams->maxSteps;
    
    size_t lpParticles, lpPsoIter, lpDimentions; // to loop thru particles, pso
    
    //check particle c value
    //if (gsl_fcmp(psoParams->c1, 0.5 + gsl_sf_log(2.0), 0.0001) == 0) {
    //if (psoParams->c1 != 0.5 + gsl_sf_log(2.0)) {
    //    fprintf(stderr, "ERROR: C1 (or c in this algorithm) is not set to the specified value in SPSO pg 7.\n");
    //    abort();
    //}
    
    
    //particles
    struct particleInfo pop[popsize];
    //neighborhoods
    const int k = 3; //number of neighbors (excluding self) that each particle has
    gsl_matrix_int *neighborhoods = gsl_matrix_int_alloc(k,psoParams->popsize); //neighbors. Collumn n contains the indexes for particle n's neighbors
    
    /* Variables needed to find and track best */
    gsl_vector *bestCoord = gsl_vector_alloc(nDim);
    double currBestFitVal = GSL_POSINF;
    double gbestFitVal = GSL_POSINF;
    gsl_vector *gbestCoord = gsl_vector_alloc(nDim);
    size_t bestfitParticle;
    
    gsl_vector *partSnrCurrCol = gsl_vector_alloc(popsize); //gets updated to the most recent set of particle fitness evaluations.
    
    int rand_neighbor; //For selecting a random neighbor
    size_t lpNbrs; /* Loop counter over nearest neighbors */
    double nbrFitVal; /*Fitness of a neighbor */
    size_t nbrIndex; // Index of a neighbor
    size_t lbestPart; /* local best particle */
    double lbestFit; /* Fitness of local best particle */
    
    //vectors for velocity update (center  of gravity)
    gsl_vector *G = gsl_vector_alloc(nDim); //center of gravity for velocity update
    double G_x_mag; //magnitude of ||G-x||
    double xprim_mag; //magnitude of x'
    gsl_vector *x_prime = gsl_vector_alloc(nDim); //xprime vector
    double rng_normal_distro; //to store random numbers from normal distribution
    double radius; //stores the chosen radius of the hypersphere
    
    double scaled_vel; //scaled velocity for confinement
    
    /* Variables needed to keep track of number of fitness function evaluations */
    unsigned char computeOK;
    size_t funcCount;
    
    //initialize all particle values (vel,position...)
    for (lpParticles=0; lpParticles < popsize; lpParticles++) {
        initPsoParticles(&pop[lpParticles], nDim, rngGen);
        pop[lpParticles].partInertia = 1.0/(2.0*gsl_sf_log(2.0)); //set inertias to value specified on pg 7
    }
    
    //start PSO loop
    for (lpPsoIter=1; lpPsoIter<maxSteps; lpPsoIter++) {
        
#ifdef HAVE_OMP
#pragma omp parallel for
#endif
        for (lpParticles = 0; lpParticles < popsize; lpParticles++){
            //Calculate G
            
            /* Evaluate fitness */
            pop[lpParticles].partSnrCurr = fitfunc(pop[lpParticles].partCoord,ffParams);
            //fprintf(stderr, "Done evaluating the fitness function...\n");
            /* Separately store all fitness values -- needed to find best particle */
            gsl_vector_set(partSnrCurrCol,lpParticles,pop[lpParticles].partSnrCurr);
            /* Check if fitness function was actually evaluated or not */
            computeOK = ((struct fitFuncParams *)ffParams)->fitEvalFlag[parallel_get_thread_num()];
            funcCount = 0;
            if (computeOK){
                /* Increment fitness function evaluation count */
                funcCount = 1;
            }
            pop[lpParticles].partFitEvals += funcCount;
            /* Update pbest fitness and coordinates if needed */
            if (pop[lpParticles].partSnrPbest > pop[lpParticles].partSnrCurr){ //minimize the function
                pop[lpParticles].partSnrPbest = pop[lpParticles].partSnrCurr;
                
                /* This is causing a segfault using openmp */
                gsl_vector_memcpy(pop[lpParticles].partPbest,pop[lpParticles].partCoord);
            }
            
        } //end paralell fitness loop
        
        /* Find the best particle in the current iteration */
        bestfitParticle = gsl_vector_min_index(partSnrCurrCol);
        currBestFitVal = pop[bestfitParticle].partSnrCurr;
        if (currBestFitVal >= gbestFitVal || lpPsoIter==1){ // if this is the first iteration (need to create neighborhoods) or fitness does not improve
            //assigns particles their new random neighbors
            for (lpParticles=0;lpParticles<(neighborhoods->size2); lpParticles++) {
                //for every particle
                size_t lpNeighbors;
                for (lpNeighbors=0; lpNeighbors < (neighborhoods->size1); lpNeighbors++) {
                    //for each of the particles neighbors (excluding self)
                    //assign random integer
                    rand_neighbor = (int)(gsl_rng_uniform_int(rngGen, neighborhoods->size2));
                    gsl_matrix_int_set(neighborhoods,lpNeighbors,lpParticles, rand_neighbor);
                }
            } //end particle for loop
            
        } else { // if PSO did improve during an iteration...
            /* Update gbest */
            gbestFitVal = pop[bestfitParticle].partSnrCurr;
            gsl_vector_memcpy(gbestCoord,pop[bestfitParticle].partCoord);
        }
        
        //update lbest values
        for (lpParticles = 0; lpParticles < popsize; lpParticles++) {
            lbestPart = lpParticles; //start with the current particle's new fitness, as all particles are the neighbors of themselves.
            lbestFit = gsl_vector_get(partSnrCurrCol, lbestPart); // get particles fitness
            for (lpNbrs = 0; lpNbrs < k; lpNbrs++) { //loop over all of the neighbors
                nbrIndex = (size_t)(gsl_matrix_int_get(neighborhoods, lpNbrs, lpParticles)); //get index of neighbor particle from matrix. Convert to size_t to match lbestPart
                nbrFitVal = gsl_vector_get(partSnrCurrCol, nbrIndex);
                if (nbrFitVal < lbestFit) { // if this neighbor has the best so far
                    lbestPart = nbrIndex;
                    lbestFit = nbrFitVal; //assign to next best
                } //end  if
            } // end for loop over neighbors
            if (lbestFit < pop[lpParticles].partSnrLbest) { //if we have discovered a position that is better in the neighborhood
                pop[lpParticles].partSnrLbest =  lbestFit;
                gsl_vector_memcpy(pop[lpParticles].partLocalBest, pop[lbestPart].partCoord);
            }
        } //end particle neighbor sharing loop
        
        for (lpParticles = 0; lpParticles < popsize; lpParticles++) {
            
            if (gsl_vector_equal(pop[lpParticles].partPbest , pop[lpParticles].partLocalBest)) { //if the local best is the pbest (G calc changes)
                gsl_vector_memcpy(G, pop[lpParticles].partPbest); //copy pbest to G to begin calculation
                gsl_vector_sub(G, pop[lpParticles].partCoord); //subtract position from pbest
                gsl_vector_scale(G, (psoParams->c1)/2.0); //scale vector by c/2
                gsl_vector_add(G, pop[lpParticles].partCoord); // add to the current position
            } else {
                //if lbest != pbest use appropriate G calculation
                gsl_vector_memcpy(G, pop[lpParticles].partCoord); //copy current position to begin calculation
                gsl_vector_scale(G, -2.0); //negitive to accout for subtraction
                gsl_vector_add(G, pop[lpParticles].partLocalBest); // add local best
                gsl_vector_add(G, pop[lpParticles].partPbest); //add personal best
                gsl_vector_scale(G, (psoParams->c1)/3.0); //scale by c/3
                gsl_vector_add(G, pop[lpParticles].partCoord); //add position, finishing g calculation
            } //end else
            
            G_x_mag = 0;
            for (lpDimentions=0; lpDimentions<nDim; lpDimentions++) { //calculate the ||G-pos|| for radius
                G_x_mag += gsl_pow_2(gsl_vector_get(G, lpDimentions)-gsl_vector_get(pop[lpParticles].partCoord, lpDimentions));
            }
            G_x_mag = sqrt(G_x_mag); //finish magnitude calculation
            
            //http://mathworld.wolfram.com/HyperspherePointPicking.html - Explains hypersphere point picking
            //Does not specify std of distro, but assuming 1
            xprim_mag = 0;
            for (lpDimentions = 0; lpDimentions < nDim; lpDimentions++) {
                //assign each index to random number for gauss distro
                rng_normal_distro = gsl_ran_gaussian(rngGen, 1.0);
                gsl_vector_set(x_prime, lpDimentions, rng_normal_distro);
                xprim_mag += gsl_pow_2(rng_normal_distro); // sum to find magnitude of vector
            }
            xprim_mag = sqrt(xprim_mag); //find the magnitude (sqrt of squared components)
            
            radius = gsl_rng_uniform(rngGen) * G_x_mag; //radius  is <= ||G-x|| uniform
            
            gsl_vector_scale(x_prime, radius/xprim_mag); // scale to point on hypersphere w chosen radius
            gsl_vector_add(x_prime, G); //move center of hypersphere to G. This finailizes x_prime
            
            //update velocity
            gsl_vector_scale(pop[lpParticles].partVel, pop[lpParticles].partInertia); //multiply  by  inertia
            gsl_vector_add(pop[lpParticles].partVel, x_prime); //add x_prime
            gsl_vector_sub(pop[lpParticles].partVel, pop[lpParticles].partCoord); //subtract position. Finishes vel update
            
            //update position
            gsl_vector_add(pop[lpParticles].partCoord, pop[lpParticles].partVel);
            
            
            //apply confinement. Absorbing
            for (lpDimentions = 0; lpDimentions<nDim; lpDimentions++) {
                double cur_xi = gsl_vector_get(pop[lpParticles].partCoord, lpDimentions);
                if (cur_xi > 1.0) { //if particle is above upper bound
                    gsl_vector_set(pop[lpParticles].partCoord, lpDimentions, 1.0);
                    scaled_vel = gsl_vector_get(pop[lpParticles].partVel, lpDimentions) * -0.5;
                    gsl_vector_set(pop[lpParticles].partVel, lpDimentions, scaled_vel);
                } else if  (cur_xi < 0.0) { //if particle is below lower bound
                    gsl_vector_set(pop[lpParticles].partCoord, lpDimentions, 0.0);
                    scaled_vel = gsl_vector_get(pop[lpParticles].partVel, lpDimentions) * -0.5;
                    gsl_vector_set(pop[lpParticles].partVel, lpDimentions, scaled_vel);
                }  //end else if for bc
            } // end constraint for loop
            
        }// end particle loop
        
    } //end pso loop
    
    //prep all of those outputs
    psoResults->bestFitVal = gbestFitVal;
    gsl_vector_memcpy(psoResults->bestLocation, gbestCoord);
    
    psoResults->computationTimeSecs = ((double) (clock() - time_start)) / CLOCKS_PER_SEC;
    psoResults->totalIterations = lpPsoIter-1;
    psoResults->totalFuncEvals = 0;
    
    for (lpParticles = 0; lpParticles < popsize; lpParticles ++){
        psoResults->totalFuncEvals += pop[lpParticles].partFitEvals;
    }
    
    //free all the vectors!
    gsl_vector_free(G);
    gsl_vector_free(x_prime);
    gsl_vector_free(bestCoord);
    gsl_vector_free(partSnrCurrCol);
    gsl_matrix_int_free(neighborhoods);
    /* Deallocate members of pop */
    for(lpParticles = 0; lpParticles < popsize; lpParticles++){
        particleinfo_free(&pop[lpParticles]);
    }
} //end function def


