#include <Python.h>
#include "pftypes.h"
#include "pmatrices.h"
#include "data.h"
#include "part.h"

data *newData(int nTax, int nParts)
{
	data *theData;
	int  i;

	theData = malloc(sizeof(data));
	theData->nTax = nTax;
	theData->nParts = nParts;
	theData->parts = malloc(nParts * sizeof(part *));
	for(i = 0; i < nParts; i++) {
		theData->parts[i] = NULL;
	}
	theData->unconstrainedLogLike = 0.0;
	return theData;
}

void freeData(data *theData)
{
	// note that Python parts created cParts,  
	// and so are responsible for free-ing cParts
	free(theData->parts);
	theData->parts = NULL;
	//printf("data.c: free-ing data %li\n", (long int)theData);
	free(theData);
	theData = NULL;
	//printf("data.c: finished freeData\n");
}

void dumpData(data *theData)
{
	int i;

	printf("                Data dump (%li)\n", (long int)theData);
	printf("                    nTax = %i\n", theData->nTax);
	printf("                    nParts = %i\n", theData->nParts);
	printf("                    parts = %li\n", (long int)theData->parts);
	if(theData->nParts < 25) {
		for(i = 0; i < theData->nParts; i++) {
			printf("                            parts[%i] = %li\n", i, (long int)theData->parts[i]);
		}
	}
}


PyObject *rell(int nBoots, rellStuff *rStuff)
{
	PyObject *thePyList;
	int i,j,k;
	int *winners;
	double *sums;
	int  *indxs;
	float theMax;
	int   theIndx;
	
	winners = pivector(rStuff->nTrees);
	indxs = pivector(rStuff->nChar);
	sums = pdvector(rStuff->nTrees);

	for(j = 0; j < rStuff->nTrees; j++) {
		winners[j] = 0;
	}	

	for(i = 0; i < nBoots; i++) {
		
		for(j = 0; j < rStuff->nChar; j++) {
			indxs[j] = (int)gsl_rng_uniform_int(rStuff->gsl_rng, (unsigned long int)(rStuff->nChar));
		}
		for(j = 0; j < rStuff->nTrees; j++) {
			sums[j] = 0.0;
		}
		for(j = 0; j < rStuff->nTrees; j++) {
			for(k = 0; k < rStuff->nChar; k++) {
				sums[j] = sums[j] + rStuff->mat[j][indxs[k]];
			}
		}
		theMax = sums[0];
		theIndx = 0;
		for(j = 0; j < rStuff->nTrees; j++) {
			if(sums[j] > theMax) {
				theMax = sums[j];
				theIndx = j;
			}
		}
		
		winners[theIndx] = winners[theIndx] + 1;
	}
		
	

	   
	thePyList = PyList_New(rStuff->nTrees);
	for(i = 0; i < rStuff->nTrees; i++) {
		PyList_SetItem(thePyList, i, PyInt_FromLong((long int)winners[i]));
	}
	free(winners);
	free(indxs);
	free(sums);
	return thePyList;
}

void bootstrapData(data *reference, data *toFill, const gsl_rng *gsl_rng)
{
	int     partNum, seqNum, pos, ran;
	part   *referencePart = NULL;
	part   *toFillPart = NULL;

#if 0
	// test that we are getting random numbers
	for(partNum = 0; partNum < 10; partNum++) {
		ran = (int)gsl_rng_uniform_int(gsl_rng, (unsigned long int)(10));
		printf("%i ", ran);
	}
	printf("\n");
	//return;
#endif	

	for(partNum = 0; partNum < toFill->nParts; partNum++) {
		referencePart = reference->parts[partNum];
		toFillPart = toFill->parts[partNum];
		for(pos = 0; pos < toFillPart->nChar; pos++) {
			ran = (int)gsl_rng_uniform_int(gsl_rng, (unsigned long int)(toFillPart->nChar));
			//printf("pos=%i, \tran=%i\n", pos,ran);
			//if(pos < 5) {
			//	printf("pos=%i, \tran=%i\n", pos,ran);
			//}
			for(seqNum = 0; seqNum < toFillPart->nTax; seqNum++) {
				toFillPart->sequences[seqNum][pos] = referencePart->sequences[seqNum][ran];
			}
		}
		makePatterns(toFillPart);
	}

}
