#include <stdio.h>  // for printf
#include <stdlib.h> // for exit
#include "defines.h" // for N_LIKE

void logDetFillFxy(int *nUnambig, int *nAmbig, int *nDoubleGap, char *seq, int sNum1, int sNum2, int *refUnambigCountMatrix, int *refAmbigCountMatrix, int *allSymbolNums, int nTax, int nChar, int dim, int bigDim, double *normUnambig, double *bigFxy, int *equatesArray)
{
	int i, j, i1, j1;
	int firstChar, secondChar, firstIndx, secondIndx;
	int na, nSlots;
	double nUnambigD, fsum, oneOverNSlots;
	int dimSquared = dim * dim;

	//printf("logDetFillFxy here.\n");

	nUnambig[0] = 0;
	nAmbig[0] = 0;
	nDoubleGap[0] = 0;


	for(i = 0; i < nChar; i++) {
		firstChar = (int)seq[(sNum1 * nChar) + i];
		secondChar = (int)seq[(sNum2 * nChar) + i];
		if(firstChar >= 0 && secondChar >= 0) {
			refUnambigCountMatrix[(firstChar * dim) +secondChar] += 1;
			nUnambig[0] += 1;
		} 
		else {
			if(firstChar == N_LIKE && secondChar == N_LIKE) {
				nDoubleGap[0] += 1;
			} 
			else {
				firstIndx = 0;
				secondIndx = 0;
				for(j1 = 0; j1 < bigDim; j1++) {
					if(allSymbolNums[j1] == firstChar) {
						firstIndx = j1;
						break;
					}
				}
				for(j1 = 0; j1 < bigDim; j1++) {
					if(allSymbolNums[j1] == secondChar) {
						secondIndx = j1;
						break;
					}
				}
				refAmbigCountMatrix[(firstIndx * bigDim) + secondIndx] += 1;
				nAmbig[0] += 1;
			}
		}
	}
    //printf("nAmbig=%i, nUnambig=%i, nDoubleGap=%i, nChar=%i\n", nAmbig[0], nUnambig[0], nDoubleGap[0], nChar);
	if((nAmbig[0] + nUnambig[0] + nDoubleGap[0]) != nChar) {
		printf("problem with Fxy\n");
		exit(0);
	}
	if(!nUnambig[0]) {
		printf("logDetFillFxy().  No unambiguous chars.\n");
		exit(0);
	}
	nUnambigD = (double)nUnambig[0];
	for(i = 0; i < dimSquared; i++) {
		normUnambig[i] = refUnambigCountMatrix[i] / nUnambigD;
	}
	if(0) {
		printf("refUnambigCountMatrix = \n");
		for(i = 0; i < dim; i++) {
			for(j = 0; j < dim; j++) {
				printf("%6i", refUnambigCountMatrix[(dim * i) + j]);
			}
			printf("\n");
		}
	}
	if(0) {
		printf("normUnambig = \n");
		for(i = 0; i < dim; i++) {
			for(j = 0; j < dim; j++) {
				printf("%6.2f", normUnambig[(dim * i) + j]);
			}
			printf("\n");
		}
	}
	if(0) {
		printf("refAmbigCountMatrix = \n");
		for(i = 0; i < bigDim; i++) {
			for(j = 0; j < bigDim; j++) {
				printf("%6i", refAmbigCountMatrix[(bigDim * i) + j]);
			}
			printf("\n");
		}
		exit(0);
	}
	for(i = 0; i < dimSquared; i++) {
		bigFxy[i] = (double)refUnambigCountMatrix[i];
	}
	
	// ---------- if ambigs --------------
	if(nAmbig[0]) {
		for(i = 0; i < bigDim; i++) {
			for(j = 0; j < bigDim; j++) {
				na = refAmbigCountMatrix[(bigDim * i) + j];
				if(na) {
					fsum = 0.0;
					nSlots = 0;
					if(i < dim) {  // first char is a symbol
						if(j == bigDim - 1) {  // second char is N_LIKE
							for(j1 = 0; j1 < dim; j1++) {
								fsum += normUnambig[(dim * i) + j1];
								nSlots += 1;
							}
						}
						else if(j >= dim) {  // second char is an equate
							for(j1 = 0; j1 < dim; j1++) {
								if(equatesArray[(dim * (j - dim)) + j1]) {
									fsum += normUnambig[(dim * i) + j1];
									nSlots += 1;
								}
							}
						} 
						else {
							printf("error in logDetFillFxy\n");
							exit(0);
						}
					}
					else if(i == bigDim - 1) {  // first char is N_LIKE
						if(j < dim) {  // second char is a symbol
							for(i1 = 0; i1 < dim; i1++) {
								fsum += normUnambig[(dim * i1) + j];
								nSlots += 1;
							}
						} 
						else if(j == bigDim - 1) {  // second char is N_LIKE
							printf("error in logDetFillFxy\n");
							exit(0);
						}
						else {  // second char is an equate
							for(i1 = 0; i1 < dim; i1++) {
								for(j1 = 0; j1 < dim; j1++) {
									if(equatesArray[(dim * (j - dim)) + j1]) {
										fsum += normUnambig[(dim * i1) + j1];
										nSlots += 1;
									}
								}
							}
						}
					}
					else {  // first char is an equate
						if(j < dim) {  // second char is a symbol
							for(i1 = 0; i1 < dim; i1++) {
								if(equatesArray[(dim * (i - dim)) + i1]) {
									fsum += normUnambig[(dim * i1) + j];
									nSlots += 1;
								}
							}
						} 
						else if(j == (bigDim - 1)) { // second char is N_LIKE
							for(i1 = 0; i1 < dim; i1++) {
								if(equatesArray[(dim * (i - dim)) + i1]) {
									for(j1 = 0; j1 < dim; j1++) {
										fsum += normUnambig[(dim * i1) + j1];
										nSlots += 1;
									}
								}
							}
						}
						else {  // second char is an equate
							for(i1 = 0; i1 < dim; i1++) {
								if(equatesArray[(dim * (i - dim)) + i1]) {
									for(j1 = 0; j1 < dim; j1++) {
										if(equatesArray[(dim * (j - dim)) + j1]) {
											fsum += normUnambig[(dim * i1) + j1];
											nSlots += 1;
										}
									}
								}
							}
						}
					}

					//printf("fsum = %f, nSlots=%i\n", fsum, nSlots);

					if(fsum < 1.e-15) { // ie equals zero
						oneOverNSlots = 1.0/(double)nSlots;
						if(i < dim) {  // first char is a symbol
							if(j == (bigDim - 1)) {  // second char is N_LIKE
								for(j1 = 0; j1 < dim; j1++) {
									bigFxy[(dim * i) + j1] += na * oneOverNSlots;
								}
							}
							else if(j >= dim) {  // second char is an equate
								for(j1 = 0; j1 < dim; j1++) {
									if(equatesArray[(dim * (j - dim)) + j1]) {
										bigFxy[(dim * i) + j1] += na * oneOverNSlots;
									}
								}
							}
							else { // second char is a symbol
								printf("error in logDetFillFxy\n");
								exit(0);
							}
						}
						else if(i == (bigDim - 1)) {  // first char is N_LIKE
							if(j < dim) {   // second char is a symbol
								for(i1 = 0; i1 < dim; i1++) {
									bigFxy[(dim * i1) + j] += na * oneOverNSlots;
								}
							}
							else if(j == (bigDim - 1)) {  // second char is N_LIKE
								printf("error in logDetFillFxy\n");
								exit(0);
							}
							else {   // second char is an equate
								for(i1 = 0; i1 < dim; i1++) {
									for(j1 = 0; j1 < dim; j1++) {
										if(equatesArray[(dim * (j - dim)) + j1]) {
											bigFxy[(dim * i1) + j1] += na * oneOverNSlots;
										}
									}
								}
							}
						}
						else {  // first char is an equate
							if(j < dim) {  // second char is a symbol
								for(i1 = 0; i1 < dim; i1++) {
									if(equatesArray[(dim * (i - dim)) + i1]) {
										bigFxy[(dim * i1) + j] += na * oneOverNSlots;
									}
								}
							}
							else if(j == (bigDim - 1)) { // second char is N_LIKE
								for(i1 = 0; i1 < dim; i1++) {
									if(equatesArray[(dim * (i - dim)) + i1]) {
										for(j1 = 0; j1 < dim; j1++) {
											bigFxy[(dim * i1) + j1] += na * oneOverNSlots;
										}
									}
								}
							}
							else {  // second char is an equate
								for(i1 = 0; i1 < dim; i1++) {
									if(equatesArray[(dim * (i - dim)) + i1]) {
										for(j1 = 0; j1 < dim; j1++) {
											if(equatesArray[(dim * (j - dim)) + j1]) {
												bigFxy[(dim * i1) + j1] += na * oneOverNSlots;
											}
										}
									}
								}
							}
						}
					}
					else { // fsum is not zero
						if(i < dim) {  // first char is a symbol
							if(j == (bigDim - 1)) {  // second char is N_LIKE
								for(j1 = 0; j1 < dim; j1++) {
									bigFxy[(dim * i) + j1] += na * normUnambig[(dim * i) + j1] / fsum;
								}
							}
							else if(j >= dim) {  // second char is an equate
								for(j1 = 0; j1 < dim; j1++) {
									if(equatesArray[(dim * (j - dim)) + j1]) {
										bigFxy[(dim * i) + j1] += na * normUnambig[(dim * i) + j1] / fsum;
									}
								}
							}
							else { // second char is a symbol
								printf("error in logDetFillFxy\n");
								exit(0);
							}
						}
						else if(i == (bigDim - 1)) {  // first char is N_LIKE
							if(j < dim) {  // second char is a symbol
								for(i1 = 0; i1 < dim; i1++) {
									bigFxy[(dim * i1) + j] += na * normUnambig[(dim * i1) + j] / fsum;
								}
							}
							else if(j == (bigDim - 1)) { // second char is N_LIKE
								printf("error in logDetFillFxy\n");
								exit(0);
							}
							else {  // second char is an equate
								for(i1 = 0; i1 < dim; i1++) {
									for(j1 = 0; j1 < dim; j1++) {
										if(equatesArray[(dim * (j - dim)) + j1]) {
											bigFxy[(dim * i1) + j1] += na * normUnambig[(dim * i1) + j1] / fsum;
										}
									}
								}
							}
						}
						else {  // first char is an equate
							if(j < dim) {  // second char is a symbol 
								for(i1 = 0; i1 < dim; i1++) {
									if(equatesArray[(dim * (i - dim)) + i1]) {
										bigFxy[(dim * i1) + j] += na * normUnambig[(dim * i1) + j] / fsum;
									}
								}
							}
							else if(j == (bigDim - 1)) {  // second char is N_LIKE
								for(i1 = 0; i1 < dim; i1++) {
									if(equatesArray[(dim * (i - dim)) + i1]) {
										for(j1 = 0; j1 < dim; j1++) {
											bigFxy[(dim * i1) + j1] += na * normUnambig[(dim * i1) + j1] / fsum;
										}
									}
								}
							}
							else { // second char is an equate
								for(i1 = 0; i1 < dim; i1++) {
									if(equatesArray[(dim * (i - dim)) + i1]) {
										for(j1 = 0; j1 < dim; j1++) {
											if(equatesArray[(dim * (j - dim)) + j1]) {
												bigFxy[(dim * i1) + j1] += na * normUnambig[(dim * i1) + j1] / fsum;
											}
										}
									}
								}
							}
						}
					} // end else fsum is not zero
					if(0) {
						printf("bigFxy, after partial ambig resolution\n");
						for(i1 = 0; i1 < dim; i1++) {
							for(j1 = 0; j1 < dim; j1++) {
								printf("%7.3f", bigFxy[(dim * i1) + j1]);
							}
							printf("\n");
						}
					}
				} // end if(na)
			}
		}
	} // ---------- end of if(ambigs) ------------
}
