#include <Python.h>
#include "pftypes.h"
#include "util.h"
#include <stdio.h>
//#include <libc.h>   // for random and srandom, but on linux srandom is in stdlib.h
#include <stdlib.h>
#include "pmatrices.h"


void checkpoint(char *location)
{
    FILE *dfile;

    if ((dfile = fopen("debugFile", "a+")) == NULL) {
        printf("checkpoint error: Couldn't open debugFile for appending \n");
        return;
    }

    fprintf(dfile, "%s\n", location);
    fclose(dfile);

}

void checkpointOneInt(char *location, int theInt)
{
    FILE *dfile;

    if ((dfile = fopen("debugFile", "a+")) == NULL) {
        printf("checkpoint error: Couldn't open debugFile for appending \n");
        return;
    }

    fprintf(dfile, "%s %i\n", location, theInt);
    fclose(dfile);

}

double ranDoubleUpToOne(void)
{
	return ((double)random()) / ((double)((long)RAND_MAX)) ;
}

void setBigQFromRMatrixDotCharFreq(double **theBigQ, double **theRMatrix, double *charFreq, int dim)
{
    int	row, col;
    double	sum;

 	
	for(col = 0; col < dim; col++){
		for(row = 0; row < dim; row++){
			theBigQ[row][col] = theRMatrix[row][col] * charFreq[col];
		}
	}

    // set theBigQ diagonals
    for(row = 0; row < dim; row++) {
        sum = 0.0;
        for(col = 0; col < dim; col++) {
            if(row != col) sum = sum + theBigQ[row][col];
        }
        theBigQ[row][row] = -sum;
    }
#if 0
    printf("sbq: charFreq\n");
    dump_pdvector(charFreq, dim);
    printf("sbq: rMatrix\n");
    dump_psdmatrix(theRMatrix, dim);
    printf("sbq: bigQ\n");
    dump_psdmatrix(theBigQ, dim);
#endif


}	

void normalizeBigQ(double **theBigQ, double *charFreq, int dim)
{
	double  sumOfCharFreqElements;
	double	sumODE;	// sum of off-diagonal elements of bigPi . bigQ
	int	row, col;

	// first check whether charFreq has been set, by summing them all.  Equals 1?
		sumOfCharFreqElements = 0.0;
		for(row = 0; row < dim; row++) {
			sumOfCharFreqElements = sumOfCharFreqElements + charFreq[row];
		}
		if(sumOfCharFreqElements < 0.999 || sumOfCharFreqElements > 1.001) {
			printf("Model: normalizeBigQ: Something wrong with the charFreq\n");
			printf("    sumOfCharFreqElements is %f\n", sumOfCharFreqElements);
			exit(1);
		}
	
	// get sum of off-diag elements of bigPi . bigQ 
		sumODE = 0.0;
		for(row = 0; row < dim; row++) {
			for(col = 0; col < dim; col++) {
				if(row != col) {
					sumODE = sumODE + (charFreq[row] * theBigQ[row][col]);
				}
			}
		}
	
	// now multiply each element of theBigQ by the inverse of sumODE
		sumODE = 1.0 / sumODE;
		for(row = 0; row < dim; row++) {
			for(col = 0; col < dim; col++) {
				theBigQ[row][col] = theBigQ[row][col] * sumODE;
			}
		}
#if 0
    printf("nbq: bigQ\n");
    dump_psdmatrix(theBigQ, dim);
#endif

}


int indexOfIntInArray(int theInt, int *theIntArray, int arrayLength)
{
	int i;

	for(i = 0; i < arrayLength; i++) {
		if(theInt == theIntArray[i]) {
			return i;
		}
	}
	return -1;
}

#if 0
long int pdbmalloc(size_t size)
{
    FILE *dfile;
	void *ptr;

    if ((dfile = fopen("pdbmallocFile", "a+")) == NULL) {
        printf("pdbmalloc error: Couldn't open pdbmallocFile for appending \n");
        return;
    }

	ptr = malloc(size);
    fprintf(dfile, "%li\n", (long int)ptr);
    fclose(dfile);
	return (long int)ptr;
}


void pdbfree(void *ptr)
{

    FILE *dfile;

    if ((dfile = fopen("pdbfreeFile", "a+")) == NULL) {
        printf("pdbfree error: Couldn't open pdbfreeFile for appending \n");
        return;
    }

    fprintf(dfile, "%li\n", ptr);
    fclose(dfile);
	free(ptr);
}
#endif


double logOfSum(double *theLogs, int len)
{
	double theMax;
	double diff;
	double mySum;
	double x;
	int i;

	theMax = theLogs[0];
	for(i = 0; i < len; i++) {
		if(theLogs[i] > theMax) {
			theMax = theLogs[i];
		}
	}
	diff = 600. - theMax;
	mySum = 0.0;
	for(i = 0; i < len; i++) {
		x = theLogs[i] + diff;
		mySum += SAFE_EXP(x);
	}
	
	return log(mySum) - diff;
}
double newtonRaftery94_eqn16(double *logLikes, int len, double harmMean, double delta, int verbose)
{
	double smallenough = 0.001;
	double logP4, oldLogP4;
	int count, i;
	double *work, *doublet;
	double numRatSOL, denomRatSOL, numerator, denominator;
	double logDelta, logOneMinusDelta, logM;

	logDelta = log(delta);
	logOneMinusDelta = log(1 - delta);
	logM = log(len);
	//printf("logDelta %f, logOneMinusDelta %f, logM %f\n", logDelta, logOneMinusDelta, logM);

	work = pdvector(len);
	doublet = pdvector(2);

	logP4 = harmMean; // to start
	oldLogP4 = logP4;

	if(verbose) {
		printf("Starting the iteration with logP4= %f\n", logP4);
	}
	count = 0;
	while(1) {
		// numerator summation
		for(i = 0; i < len; i++) {
			doublet[0] = logDelta + logP4;
			doublet[1] = logOneMinusDelta + logLikes[i];
			work[i] = logLikes[i] - logOfSum(doublet, 2);
		}
		numRatSOL = logOfSum(work, len);

		// denominator summation
		for(i = 0; i < len; i++) {
			doublet[0] = logDelta + logP4;
			doublet[1] = logOneMinusDelta + logLikes[i];
			work[i] = - logOfSum(doublet, 2);
		}
		denomRatSOL = logOfSum(work, len);

		doublet[0] = logDelta + logM - logOneMinusDelta;
		doublet[1] = numRatSOL;
		numerator = logOfSum(doublet, 2);

		doublet[0] = logDelta + logM - (logOneMinusDelta + logP4);
		doublet[1] = denomRatSOL;
		denominator = logOfSum(doublet, 2);
		//printf("count %3i  numRatSOL %f, denomRatSOL %f, numerator %f, denominator %f \n", 
		//	   count, numRatSOL, denomRatSOL, numerator, denominator);

		count += 1;
		oldLogP4 = logP4;
		logP4 = numerator - denominator;
		if(fabs(logP4 - oldLogP4) < smallenough) {
			break;
		}
		if(verbose && (count % 100 == 0)) {
			printf("looped %i times. current logP4 is %f, hasnt converged yet.\n", count, logP4);
		}

		//break;
	}
	if(verbose) {
		printf("Looped %i times.  converged on logP4= %f\n", count, logP4);
	}

	free(work);
	free(doublet);
	return logP4;
}

