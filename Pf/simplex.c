#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "pmatrices.h"
#include "simplex.h"
#include "util.h"

//#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
#define SIMPLEX_VERBOSE 0
#define JIGGLE_FACTOR 1.0e-9
//#define START_FACTOR 0.85

simplex *newSimplex(
					int      dim,
					double   *prams,
					double   *lowerBounds,
					double   *upperBounds,
					double   yDiffRequested,
					int      nMaxEvaluationsAllowed,
					double   (*funk)(double []),
					int      *compStarts,
					double   startFactor)
{
	simplex  *aSimp;
	int i,j;
	//double temp;

	if(startFactor < 0.0001 || startFactor > 0.4) {
		printf("startFactor should be between 0.0001 and 0.4\n");
		exit(1);
	}
	aSimp = malloc(sizeof(simplex));
	if(!aSimp) {
		printf("failed to malloc simplex. exiting...\n");
		exit(1);
	}
	aSimp->dim = dim;
	aSimp->prams = prams;
	//printf("inside simplex, prams is at address %li, prams[0] (%li) is %f\n", 
	//	   (long int)prams, (long int) &(prams[0]), prams[0]);
	aSimp->lowerBounds = lowerBounds;
	aSimp->upperBounds = upperBounds;
	aSimp->yDiffRequested = yDiffRequested;
	aSimp->nMaxEvaluationsAllowed = nMaxEvaluationsAllowed;
	aSimp->funk = funk;
	aSimp->compStarts = compStarts;

	aSimp->p = pdmatrix(dim + 1, dim);
	aSimp->y = pdvector(dim + 1);
	aSimp->funkEvaluationCount = 0;
	aSimp->pSum = pdvector(dim);
	aSimp->pTry = pdvector(dim);

	for(j = 0; j < dim; j++) {  // for each parameter
		if(upperBounds[j] <= lowerBounds[j]) {
			printf("simplex.c: problem with bounds for (zero-based) parameter %i\n", j);
			exit(1);
		}
		for(i = 0; i < dim + 1; i++) { // for each point
			if(startFactor > 0.05) {
				aSimp->p[i][j] = prams[j] * ((0.95 - startFactor) + (ranDoubleUpToOne() * startFactor));
			} else if(startFactor > 0.005) {
				aSimp->p[i][j] = prams[j] * ((0.99 - startFactor) + (ranDoubleUpToOne() * startFactor));
			} else {
				aSimp->p[i][j] = prams[j] * ((0.999 - startFactor) + (ranDoubleUpToOne() * startFactor));
			}
		}
	}
	for(i = 0; i < dim; i++) {
		aSimp->p[i + 1][i] = prams[i];
	}

	for(i = 0; i < dim + 1; i++) {
		aSimp->y[i] = funk(aSimp->p[i]);
		//printf("assigning y[%i] to %f using params %f %f %f\n",
		//	   i,  funk(aSimp->p[i]), aSimp->p[i][0], aSimp->p[i][1], aSimp->p[i][2]);
		aSimp->funkEvaluationCount++;
	}

#if SIMPLEX_VERBOSE
	printf("\nnewSimplex: initial p values, followed by y values..\n");
	for(i = 0; i < dim + 1; i++) {
		for(j = 0; j < dim; j++) {
			printf("%7.3f", aSimp->p[i][j]);
		}
		printf(":  %12.5f", aSimp->y[i]);
		printf("\n");
	}

	//printf("newSimplex: initial y values are\n");
	//for(i = 0; i < dim + 1; i++) {
	//	printf("%7.3f", aSimp->y[i]);
	//}
	printf("\n\n");
#endif	

	return aSimp;
}

void freeSimplex(simplex *aSimp)
{
	if(aSimp->p) free_pdmatrix(aSimp->p);
	aSimp->p = NULL;
	if(aSimp->y) free(aSimp->y);
	aSimp->y = NULL;
	if(aSimp->pSum) free(aSimp->pSum);
	aSimp->pSum = NULL;
	if(aSimp->pTry) free(aSimp->pTry);
	aSimp->pTry = NULL;
	free(aSimp);
	aSimp = NULL;
}



double amoeba(simplex *aSimp)
{
	int i, ihi, ilo, inhi, j;
	int mpts = aSimp->dim + 1;
	double rtol, sum, ysave, ytry;
	//double swap;

	for (j = 0; j < aSimp->dim; j++) {
		sum = 0.0;
		for (i = 0; i < mpts; i++) {
			sum = sum +  aSimp->p[i][j];
		}
		aSimp->pSum[j] = sum;
		//printf("pSum[%i] = %f\n", j, aSimp->pSum[j]);
	}
	
	// ihi is index of highest, inhi is index of next highest, ilo is index of lowest.
	for (;;) {  
		ilo = 0;
		if(aSimp->y[0] > aSimp->y[1]) {
			ihi = 0;
			inhi = 1;
		} else {
			ihi = 1;
			inhi = 0;
		}
		for (i = 0; i < mpts; i++) {
			if (aSimp->y[i] <= aSimp->y[ilo]) ilo=i;
			if (aSimp->y[i] > aSimp->y[ihi]) {
				inhi = ihi;
				ihi = i;
			} else if (aSimp->y[i] > aSimp->y[inhi] && i != ihi) inhi = i;
		}

#if SIMPLEX_VERBOSE
		//printf("ihi = %i, inhi = %i, ilo = %i\n", ihi, inhi, ilo);
		//printf("y's hi to low, are %4.3f, %4.3f, %4.3f     ", 
		//							aSimp->y[ihi], aSimp->y[inhi], aSimp->y[ilo]);
		for(i = 0; i < aSimp->dim; i++) {
			printf("%6.2f", aSimp->p[ilo][i]);
		}
		printf(", y = %6.2f    ", aSimp->y[ilo]);  // follow this with rtol, below
#endif
		
		// This is no good, cuz it depends on how close it is to zero.
		//rtol = 2.0 * fabs(y[ihi] - y[ilo])/(fabs(y[ihi]) + fabs(y[ilo]));
		rtol = fabs(aSimp->y[ihi] - aSimp->y[ilo]);
#if SIMPLEX_VERBOSE
		printf("rtol is %f \n", rtol);
#endif
		if (rtol < aSimp->yDiffRequested || 
			aSimp->funkEvaluationCount >= aSimp->nMaxEvaluationsAllowed) {  // finished!
			if(aSimp->funkEvaluationCount >= aSimp->nMaxEvaluationsAllowed) {
				printf("simplex funkEvaluationCount maxed out.\n");
			}
			for(i = 0; i < aSimp->dim; i++) {
				aSimp->prams[i] = aSimp->p[ilo][i];
			}
			//SWAP(aSimp->y[0],aSimp->y[ilo])
			//for (i = 0; i < aSimp->dim; i++) SWAP(aSimp->p[0][i],aSimp->p[ilo][i])
#if SIMPLEX_VERBOSE
	        printf("\nsimplex finished: p values are\n");
			for(j = 0; j < aSimp->dim; j++) {
				printf("%7.3f", aSimp->p[ilo][j]);
			}
			printf(", \n    and y = %7.3f\n", aSimp->y[ilo]);
			printf("    Did %i function evaluations\n", aSimp->funkEvaluationCount);
#endif	
			return aSimp->y[ilo];
		}

		ytry = amotry(aSimp, ihi, -1.0);
		//printf("amotry returned ytry = %f\n", ytry);
		if (ytry <= aSimp->y[ilo]) {
			//printf(" ytry (%f) is less than y[ilo] = y[%i] = %f\n", ytry, ilo, aSimp->y[ilo]);
			// its the right direction, so go twice as far...
			ytry = amotry(aSimp, ihi, 2.0); 
		}
		else if (ytry >= aSimp->y[inhi]) {  // is it higher than the second highest?
			//printf(" ytry (%f) is greater than y[ilo] = y[%i] = %f\n", ytry, ilo, aSimp->y[ilo]);
			ysave = aSimp->y[ihi];  // the worst point
			ytry = amotry(aSimp, ihi, 0.5);
			if (ytry >= ysave) {   // did it just get worse than the worst point?
				for (i = 0; i < mpts; i++) {                               // for each point
					if (i != ilo) {                                        // except the best point
						for (j = 0; j < aSimp->dim; j++)                   // replace the old point with
							aSimp->p[i][j] = aSimp->pSum[j] =              // the average of the old
								0.5 * (aSimp->p[i][j] + aSimp->p[ilo][j]); // point and the best point
						// I don't think that we need to worry about bounds or composition
						// issues here, because each old point and the best point must all have
						// been ok, so the new points, being averages, would also be ok.  Right?
						aSimp->y[i] = aSimp->funk(aSimp->pSum);
						aSimp->funkEvaluationCount++;
					}
				}
				for (j = 0; j < aSimp->dim; j++) {
					sum = 0.0;
					for (i = 0; i < mpts; i++) {
						sum = sum +  aSimp->p[i][j];
					}
					aSimp->pSum[j] = sum;
				}
			}
		}
	}
}


double amotry(simplex *aSimp, int ihi, double fac)
{
	int     i, j;
	double  fac1, fac2, ytry;
	double  compSum;
	int speak = SIMPLEX_VERBOSE;
	int compLen = 0;  // composition dimension - 1, ie nSymbols - 1

	fac1 = (1.0 - fac) / aSimp->dim;
	fac2 = fac1 - fac;
	for (j = 0; j < aSimp->dim; j++) {
		aSimp->pTry[j] = aSimp->pSum[j] * fac1 - aSimp->p[ihi][j] * fac2;
		// at this point we can check whether aSimp->pTry[j] is out of bounds
		if(aSimp->pTry[j] < aSimp->lowerBounds[j]) {
			aSimp->pTry[j] = aSimp->lowerBounds[j] + (JIGGLE_FACTOR * aSimp->pSum[j]);
			if(speak) printf("      try less than lower bounds, setting to %g instead\n", aSimp->pTry[j]);
		} else if(aSimp->pTry[j] > aSimp->upperBounds[j]) {
			aSimp->pTry[j] = aSimp->upperBounds[j] - (JIGGLE_FACTOR * aSimp->pSum[j]);
			if(speak) printf("      try greater than upper bounds, setting to %g instead\n", aSimp->pTry[j]);
		}
	}

	// now deal with composition issues
	
	compSum = 0.0;
	for (i = 0; i < aSimp->dim; i++) {
		if(aSimp->compStarts[i]) {
			compLen = aSimp->compStarts[i];  // nSymbols - 1, or dim (sensu model) - 1.
			compSum = 0.0;
			for(j = i; j < i + compLen; j++) {
				compSum = compSum + aSimp->pTry[j];
			}
			if(compSum > 1.0) {
				if(speak) printf("    comp too large (total %f), scaling back\n", compSum);
				for(j = i; j < i + compLen; j++) {
					aSimp->pTry[j] = aSimp->pTry[j] / compSum;
				}
			}
		}
	}
	
	
	ytry = aSimp->funk(aSimp->pTry);
	aSimp->funkEvaluationCount++;
	if (ytry < aSimp->y[ihi]) {
		aSimp->y[ihi] = ytry;
		for (j = 0; j < aSimp->dim; j++) {
			aSimp->pSum[j] += aSimp->pTry[j] - aSimp->p[ihi][j];
			aSimp->p[ihi][j] = aSimp->pTry[j];
		}
	}
	return ytry;
}
