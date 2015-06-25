#include "pftypes.h"
#include "eig.h"
#include <stdlib.h>
#include "pmatrices.h"
#include "pcomplex.h"
#include "linalg.h"
#include <stdio.h>
#include "util.h"
#include <math.h>

eig *allocEig(int theDim, double **theMat)
{
	eig	*anEig;
	
	anEig = malloc(sizeof(eig));
	anEig->complexEig = 0;
	anEig->dim = theDim;
	anEig->mat = theMat;
	anEig->mwork = psdmatrix(theDim);
	anEig->eigvecs = psdmatrix(theDim);
	anEig->inverseEigvecs = psdmatrix(theDim);
	anEig->eigvals = pdvector(theDim);
	anEig->eigvalsImag = pdvector(theDim);
	anEig->dwork = pdvector(theDim);
	anEig->iwork = pivector(theDim);
	
	anEig->Ceigvecs = NULL;
	anEig->CinverseEigvecs = NULL;
	anEig->Cwork = NULL;
	anEig->Ccol = NULL;
	anEig->complexStuffHasBeenAllocated = 0;
	//(void) eigensystem(anEig);
	return anEig;
}

void allocComplexStuff(eig *anEig)
{
	anEig->Ceigvecs = pscmatrix(anEig->dim);
	anEig->CinverseEigvecs = pscmatrix(anEig->dim);
	anEig->Ccol = (pcomplex *)malloc(anEig->dim * sizeof(pcomplex));
	anEig->Cwork = pscmatrix(anEig->dim);
	anEig->complexStuffHasBeenAllocated = 1;
}

void freeEig(eig *anEig)
{
	//checkpoint("freeEig: begin");
	if(anEig->mwork) free_psdmatrix(anEig->mwork);
	if(anEig->eigvecs) free_psdmatrix(anEig->eigvecs);
	if(anEig->inverseEigvecs) free_psdmatrix(anEig->inverseEigvecs);
	if(anEig->eigvals) free(anEig->eigvals);
	if(anEig->eigvalsImag) free(anEig->eigvalsImag);
	if(anEig->dwork) free(anEig->dwork);
	if(anEig->iwork) free(anEig->iwork);
	if(anEig->Ceigvecs) free_pscmatrix(anEig->Ceigvecs);
	if(anEig->CinverseEigvecs) free_pscmatrix(anEig->CinverseEigvecs);
	if(anEig->Ccol) free(anEig->Ccol);
	if(anEig->Cwork) free_pscmatrix(anEig->Cwork);
	free(anEig);
	//checkpoint("freeEig: end");
}	

int eigensystem(eig *anEig)
{
	int	i, j, retval;
	
	//printf("eig eigensystem: starting\n");
	copy_psdmatrix(anEig->mat, anEig->mwork, anEig->dim);
	//checkpoint("about to do EigenRealGeneral");
	retval = EigenRealGeneral(anEig->dim, anEig->mwork, anEig->eigvals, 
							  anEig->eigvalsImag, anEig->eigvecs, anEig->iwork, anEig->dwork);
	//checkpoint("eig eigensystem: afterEigenRealGeneral");
	// zero is no error, 1 is error
	if(retval) {
		printf("EigenRealGeneral function returned a %i\n", retval);
		return retval;
	}
	
	// check for non-zeros in eigvalsImag
	for(i = 0; i < anEig->dim; i++) {
		if(anEig->eigvalsImag[i] != 0.0) {
			//printf(
			//		"anEig->eigvalsImag[%i] is %f\n", i, anEig->eigvalsImag[i]);
			anEig->complexEig = 1;
		}
		if(anEig->complexEig) break;	// don't bother looking any more, 
		// its complex
	}
	
	if(anEig->complexEig == 1 && anEig->complexStuffHasBeenAllocated == 0) {
		printf("eig eigensystem: alloc-ing complexStuff\n");
		allocComplexStuff(anEig);
	}
	
	// EigenRealGeneral eigvectors are columns of the matrix eigvecs,
	// which is what I want.
	
	// Now get the inverse of the eigvecs, stick it in inverseEigvecs
	// or CinverseEigvecs
	// But beware- the input matrix is mucked!  So copy it first
		
	if(!anEig->complexEig) {
		//printf("eig eigensystem: inverting non-complex\n");
		copy_psdmatrix(anEig->eigvecs, anEig->mwork, anEig->dim);
		//checkpoint("about to InvertMatrix");
		InvertMatrix(anEig->mwork, anEig->dim, anEig->dwork, anEig->iwork, anEig->inverseEigvecs);
	}
		
	if(anEig->complexEig) {
		// printf("eig eigensystem: inverting complex\n");
		// unconvolve the eigvecs into Ceigvecs.  
		// unusual: i is col, j is row
		for(i = 0; i < anEig->dim; i++) {
			if(anEig->eigvalsImag[i] == 0.0) { 
				for(j = 0; j < anEig->dim; j++) {
					anEig->Ceigvecs[j][i].re = anEig->eigvecs[j][i];
					anEig->Ceigvecs[j][i].im = 0.0;
				}
			} else if(anEig->eigvalsImag[i] > 0.0) { 
				for(j = 0; j < anEig->dim; j++) {
					anEig->Ceigvecs[j][i].re = anEig->eigvecs[j][i];
					anEig->Ceigvecs[j][i].im = anEig->eigvecs[j][i + 1];
				}
			} else if(anEig->eigvalsImag[i] < 0.0) { 
				for(j = 0; j < anEig->dim; j++) {
					anEig->Ceigvecs[j][i].re = anEig->eigvecs[j][i - 1];
					anEig->Ceigvecs[j][i].im = -anEig->eigvecs[j][i];
				}
			}
		}
			
		copy_pscmatrix(anEig->Ceigvecs, anEig->Cwork, anEig->dim);
		//checkpoint("about to ComplexInvertMatrix");
		ComplexInvertMatrix(anEig->Cwork, anEig->dim, anEig->dwork, anEig->iwork, 
							anEig->CinverseEigvecs, anEig->Ccol);
	}
   
#if 0
	printf("anEig->eigvals(real) ");
	dump_pdvector(anEig->eigvals, anEig->dim);
		
	//printf("anEig->eigvals(imag) ");
	//dump_pdvector(anEig->eigvalsImag, anEig->dim);
		
	printf("eigenvectors (as columns)\n");
	dump_psdmatrix(anEig->eigvecs, anEig->dim);
	if(anEig->complexEig && anEig->Ceigvecs){
		dump_pscmatrix(anEig->Ceigvecs, anEig->dim);
	}
			
	printf("Inverse of eigenvectors\n");
	if(!anEig->complexEig) {
		dump_psdmatrix(anEig->inverseEigvecs, anEig->dim);
	} else {
		dump_pscmatrix(anEig->CinverseEigvecs, anEig->dim);
	}
#endif

	return retval;

}

void matrixExpTimesBranchLength(eig *anEig, double branchLength, double **result)
{
	int	i, j, k;
	//double prod;
	//printf("      --------- here\n"); fflush(stdout);
	if(!anEig->complexEig) {
		//printf("not complex\n"); fflush(stdout);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				//printf("i=%i, j=%i\n", i, j); fflush(stdout);
				result[i][j] = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					// Mar 15_00, I found this test necessary on the dec alpha, at least
					/*
					  prod = anEig->eigvals[k] * branchLength;
					  if(prod > -700.0) {
					  result[i][j] = result[i][j] + (anEig->eigvecs[i][k] * 
					  anEig->inverseEigvecs[k][j] * 
					  exp(anEig->eigvals[k] * branchLength));
					  }
					*/
					// Aug8 replaced the above hack with this, stolen from bambe
					result[i][j] = result[i][j] + (anEig->eigvecs[i][k] *
												   anEig->inverseEigvecs[k][j] *
												   SAFE_EXP(anEig->eigvals[k] * branchLength));
				}
			}
		}
	}

	if(anEig->complexEig) {
		//printf("is complex\n"); fflush(stdout);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				anEig->Cwork[i][j].re = 0.0;
				anEig->Cwork[i][j].im = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					anEig->Cwork[i][j] = Cadd(anEig->Cwork[i][j], Cmul(anEig->Ceigvecs[i][k], 
																	   Cmul(anEig->CinverseEigvecs[k][j], 
																			Cexp(RCmul(branchLength, 
																					   Complex(anEig->eigvals[k], anEig->eigvalsImag[k]))))));
					
								
				}
			}
		}
		//printf("Eig: matrixExpTimes: ...\n"); 
		//dump_pscmatrix(anEig->Cwork, anEig->dim);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = anEig->Cwork[i][j].re;
				if(fabs(anEig->Cwork[i][j].im) > 0.000000000001) {
					printf("Eig: matrixExpTimes: imaginary numbers\n");
				}
			}
		}

	}
	

}

//- (void) matrixLogInto: (double **) result;
// this could be speeded up lots by taking the logs of the eigsvals first and re-using them
void matrixLog(eig *anEig, double **result)
{
	int	i, j, k;
	//complex temp;
	
	if(!anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					result[i][j] = result[i][j] + (anEig->eigvecs[i][k] * 
												   anEig->inverseEigvecs[k][j] * 
												   log(anEig->eigvals[k]));
				}
			}
		}
	}
	
	if(anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				anEig->Cwork[i][j].re = 0.0;
				anEig->Cwork[i][j].im = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					//temp = Cmul(Ceigvecs[i][k], 
					//			Cmul(CinverseEigvecs[k][j], 
					//			Clog(Complex(eigvals[k], eigvalsImag[k]))));
					//printf("k = %i, addend = %f + %f I\n",
					//					k, temp.re, temp.im);
					anEig->Cwork[i][j] = Cadd(anEig->Cwork[i][j], Cmul(anEig->Ceigvecs[i][k], 
																	   Cmul(anEig->CinverseEigvecs[k][j], 
																			Clog(Complex(anEig->eigvals[k], anEig->eigvalsImag[k])))));
				}
			}
		}
		//printf("Eig: matrixLogInto: log of matrix is...\n"); 
		//dump_pscmatrix(Cwork, dim);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = anEig->Cwork[i][j].re;
				if(fabs(anEig->Cwork[i][j].im) > 0.00000000001) {
					printf("Eig: matrixLogInto: imaginary numbers\n");
				}
			}
		}

	}
	
}

//- (void) matrixPower:(double) pow into: (double **) result;
void matrixPower(eig *anEig, double pow, double **result)
{
	int	i, j, k;
	
	for(i = 0; i < anEig->dim; i++) {
		for(j = 0; j < anEig->dim; j++) {
			result[i][j] = 0.0;
			for(k = 0; k < anEig->dim; k++) {
				result[i][j] = result[i][j] + (anEig->eigvecs[i][k] * 
											   anEig->inverseEigvecs[k][j] * 
											   SAFE_EXP(pow * log(anEig->eigvals[k])));
			}
		}
	}
}

void firstDerivativeOfMatrixExpTimesBranchLength(eig *anEig, double branchLength, double **result, double rate)
{
	int	i, j, k;
	pcomplex	ceval;
	
	// note that the branchLength already has the rate built-in
	if(!anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					result[i][j] = result[i][j] + (anEig->eigvecs[i][k] * 
												   anEig->inverseEigvecs[k][j] * 
												   anEig->eigvals[k] *
												   rate *
												   SAFE_EXP(anEig->eigvals[k] * branchLength));
				}
			}
		}
	}
	
	if(anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				anEig->Cwork[i][j].re = 0.0;
				anEig->Cwork[i][j].im = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					ceval = Complex(anEig->eigvals[k], anEig->eigvalsImag[k]);
					anEig->Cwork[i][j] = Cadd(anEig->Cwork[i][j], Cmul(anEig->Ceigvecs[i][k], 
																	   Cmul(anEig->CinverseEigvecs[k][j], Cmul(RCmul(rate, ceval),
																											   Cexp(RCmul(branchLength, ceval))))));
				
				}
			}
		}
		//printf("Eig: firstDerivativeOfMatrixExpTimes: "
		//			"first derivative of matrix is...\n"); 
		//dump_pscmatrix(Cwork, dim);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = anEig->Cwork[i][j].re;
				if(fabs(anEig->Cwork[i][j].im) > 0.00000000001) {
					printf("Eig: firstDerivative: imaginary numbers\n");
				}
			}
		}

	}
	
}

void secondDerivativeOfMatrixExpTimesBranchLength(eig *anEig, double branchLength, double **result, double rate)
{
	int	i, j, k;
	pcomplex	ceval;
	
    // note that the branchLength already has the rate built-in
	if(!anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					result[i][j] = result[i][j] + (anEig->eigvecs[i][k] * 
												   anEig->inverseEigvecs[k][j] * 
												   anEig->eigvals[k] * anEig->eigvals[k] *
												   rate * rate *
												   SAFE_EXP(anEig->eigvals[k] * branchLength));
				}
			}
		}
	}
	
	if(anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				anEig->Cwork[i][j].re = 0.0;
				anEig->Cwork[i][j].im = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					ceval = Complex(anEig->eigvals[k], anEig->eigvalsImag[k]);
					anEig->Cwork[i][j] = Cadd(anEig->Cwork[i][j], Cmul(anEig->Ceigvecs[i][k], 
																	   Cmul(anEig->CinverseEigvecs[k][j], 
																			Cmul(ceval, Cmul(RCmul(rate * rate,ceval),
																							 Cexp(RCmul(branchLength, ceval)))))));
				
				}
			}
		}
		//printf("Eig: secondDerivativeOfMatrixExpTimes: "
		//			"second derivative of matrix is...\n"); 
		//dump_pscmatrix(Cwork, dim);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = anEig->Cwork[i][j].re;
				if(fabs(anEig->Cwork[i][j].im) > 0.00000000001) {
					printf("Eig: secondDerivative: imaginary numbers\n");
				}
			}
		}

	}

}

#if 0

void firstDerivativeOfMatrixExpTimesBranchLengthL(eig *anEig, double branchLength, double **result, double rate)
{
	int	i, j, k;
	pcomplex	ceval;
	
	// note that the branchLength already has the rate built-in
	if(!anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					result[i][j] = result[i][j] + (anEig->eigvecs[i][k] * 
												   anEig->inverseEigvecs[k][j] * 
												   anEig->eigvals[k] *
												   rate * branchLength *  // change is here
												   SAFE_EXP(anEig->eigvals[k] * branchLength));
				}
			}
		}
	}
	
	if(anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				anEig->Cwork[i][j].re = 0.0;
				anEig->Cwork[i][j].im = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					ceval = Complex(anEig->eigvals[k], anEig->eigvalsImag[k]);
					anEig->Cwork[i][j] = Cadd(anEig->Cwork[i][j], 
											  Cmul(anEig->Ceigvecs[i][k], 
												   Cmul(anEig->CinverseEigvecs[k][j], 
														Cmul(RCmul(rate, ceval), 
															 Cexp(RCmul(branchLength, ceval))))));
				
				}
			}
		}
		//printf("Eig: firstDerivativeOfMatrixExpTimes: "
		//			"first derivative of matrix is...\n"); 
		//dump_pscmatrix(Cwork, dim);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = anEig->Cwork[i][j].re;
				if(fabs(anEig->Cwork[i][j].im) > 0.00000000001) {
					printf("Eig: firstDerivative: imaginary numbers\n");
				}
			}
		}

	}
	
}

void secondDerivativeOfMatrixExpTimesBranchLengthL(eig *anEig, double branchLength, double **result, double rate)
{
	int	i, j, k;
	pcomplex	ceval;
	
    // note that the branchLength already has the rate built-in
	if(!anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					result[i][j] = result[i][j] + (anEig->eigvecs[i][k] * 
												   anEig->inverseEigvecs[k][j] * 
												   (anEig->eigvals[k] * anEig->eigvals[k] *
													rate * rate * branchLength * branchLength +
													(anEig->eigvals[k] * branchLength)) *
												   SAFE_EXP(anEig->eigvals[k] * branchLength));
				}
			}
		}
	}
	
	if(anEig->complexEig) {
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				anEig->Cwork[i][j].re = 0.0;
				anEig->Cwork[i][j].im = 0.0;
				for(k = 0; k < anEig->dim; k++) {
					ceval = Complex(anEig->eigvals[k], anEig->eigvalsImag[k]);
					anEig->Cwork[i][j] = Cadd(anEig->Cwork[i][j], Cmul(anEig->Ceigvecs[i][k], 
																	   Cmul(anEig->CinverseEigvecs[k][j], 
																			Cmul(ceval, Cmul(RCmul(rate * rate,ceval),
																							 Cexp(RCmul(branchLength, ceval)))))));
				
				}
			}
		}
		//printf("Eig: secondDerivativeOfMatrixExpTimes: "
		//			"second derivative of matrix is...\n"); 
		//dump_pscmatrix(Cwork, dim);
		for(i = 0; i < anEig->dim; i++) {
			for(j = 0; j < anEig->dim; j++) {
				result[i][j] = anEig->Cwork[i][j].re;
				if(fabs(anEig->Cwork[i][j].im) > 0.00000000001) {
					printf("Eig: secondDerivative: imaginary numbers\n");
				}
			}
		}

	}

}

#endif
