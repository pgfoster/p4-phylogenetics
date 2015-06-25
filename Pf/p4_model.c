#include "pftypes.h"
#include <Python.h>
//#include <Numeric/arrayobject.h>
#include <numpy/arrayobject.h>
#include "p4_model.h"
#include "pmatrices.h"
#include "eig.h"
#include "proteinModels.h"
#include "defines.h"
#include "util.h"


p4_model *p4_newModel(int nParts, int doRelRates, int relRatesAreFree, int nFreePrams, int isHet, int *rMatrixNormalizeTo1)
{
	p4_model  *aModel = NULL;

	//printf("p4_newModel here.  nParts=%i\n", nParts);
	aModel = malloc(sizeof(p4_model));
	if(!aModel) {
		printf("Failed to alloc memory for model\n");
		exit(1);
	}
		
	aModel->nParts = nParts;
	aModel->parts = malloc(nParts * sizeof(p4_modelPart *));
	if(!aModel->parts) {
		printf("Failed to alloc memory for model->parts\n");
		exit(1);
	}

	aModel->doRelRates = doRelRates;
	aModel->relRatesAreFree = relRatesAreFree;
	aModel->nFreePrams = nFreePrams;
	aModel->isHet = isHet;
	aModel->rMatrixNormalizeTo1 = rMatrixNormalizeTo1;

	return aModel;
}

void p4_freeModel(p4_model *aModel)
{
	int i;

	//printf("p4_freeModel here.\n");
	if(aModel->parts) {
		for(i = 0; i < aModel->nParts; i++) {
			if(aModel->parts[i]) {
				p4_freeModelPart(aModel->parts[i]);
				aModel->parts[i] = NULL;
			}
		}
		free(aModel->parts);
		aModel->parts = NULL;
	}
	aModel->rMatrixNormalizeTo1 = NULL; // a numpy array
	free(aModel);
	aModel = NULL;
	//printf("freed p4_model\n");
}

void p4_dumpModel(p4_model *aModel)
{
	int i;

	printf("p4_dumpModel()\n");
	printf("  nParts %2i\n", aModel->nParts);
	for(i = 0; i < aModel->nParts; i++) {
		printf("part\n");
	}
}


//==========
// modelPart
//==========

void p4_newModelPart(p4_model *aModel, int pNum, int dim, int nComps, int nRMatrices, int nGdasrvs, int nCat, int isMixture, int mixtureIsFree, int pInvarFree, int doTSCovarion, int tSCovIsFree, int *bQETneedsReset)
{
	p4_modelPart  *mp;
	int i,j;


	//printf("..p4_newModelPart here.  pNum=%i, nComps=%i\n", pNum, nComps);
	mp = malloc(sizeof(p4_modelPart));
	if(!mp) {
		printf("Failed to alloc memory for modelPart\n");
		exit(1);
	}

	mp->doTSCovarion = doTSCovarion;
	if(doTSCovarion) {
		dim += dim;
	}
	mp->dim = dim;
	//mp->clNeedsUpdating = 0;  // init

	// comp
	mp->nComps = nComps;
	mp->comps = malloc(nComps * sizeof(p4_comp *));
	if(!mp->comps) {
		printf("Failed to alloc memory for modelPart->comps\n");
		exit(1);
	}
	for(i = 0; i < nComps; i++) {
		mp->comps[i] = NULL;
	}

	// rMatrices
	mp->nRMatrices = nRMatrices;
	mp->rMatrices = malloc(nRMatrices * sizeof(p4_rMatrix *));
	if(!mp->rMatrices) {
		printf("Failed to alloc memory for modelPart->rMatrices\n");
		exit(1);
	}
	for(i = 0; i < nRMatrices; i++) {
		mp->rMatrices[i] = NULL;
	}

	// gdasrvs
	mp->nGdasrvs = nGdasrvs;
	// it is possible that nGdasrvs is 0, if nCat is 1
	if(mp->nGdasrvs) {
		mp->gdasrvs = malloc(nGdasrvs * sizeof(p4_gdasrv *));
		if(!mp->gdasrvs) {
			printf("Failed to alloc memory for modelPart->gdasrvs\n");
			exit(1);
		}
	} else {
		mp->gdasrvs = NULL;
	}
	for(i = 0; i < nGdasrvs; i++) {
		mp->gdasrvs[i] = NULL;
	}
	

	// set isMixture and nCat
	mp->isMixture = isMixture;
	mp->nCat = nCat;
	
	// pInvar
	mp->pInvar = (p4_pInvar *)malloc(sizeof(p4_pInvar));
	if(!mp->pInvar) {
		printf("Failed to alloc memory for modelPart->pInvar\n");
		exit(1);
	}
	mp->pInvar->val = (double *)malloc(sizeof(double));
	if(!mp->pInvar->val) {
		printf("Failed to alloc memory for modelPart->pInvar->val\n");
		exit(1);
	}
	mp->pInvar->free = pInvarFree;
	mp->pInvar->val[0]  = -1.0;

	// mixture
	mp->mixture = NULL;
	if(mp->isMixture) {
		mp->mixture = (p4_mixture *)malloc(sizeof(p4_mixture));
		if(!mp->mixture) {
			printf("Failed to alloc memory for modelPart->mixture\n");
			exit(1);
		}
		mp->mixture->freqs = malloc(mp->nCat * sizeof(double));
		mp->mixture->rates = malloc(mp->nCat * sizeof(double));
		if(!mp->mixture->freqs || !mp->mixture->rates) {
			printf("Failed to alloc memory for mp->mixture->freqs/rates\n");
			exit(1);
		}
		mp->mixture->free = mixtureIsFree;
	}


	// relRate
	mp->relRate = (double *)malloc(sizeof(double));
	if(!mp->relRate) {
		printf("Failed to alloc memory for modelPart->relRate\n");
		exit(1);
	}
	mp->relRate[0] = -1;


	// bigQAndEigThing
	mp->bigQAndEigThing = (p4_bigQAndEig ***)malloc(nComps * sizeof(p4_bigQAndEig **));
	if(!mp->bigQAndEigThing) {
		printf("Failed to allocate memory for mp->bigQAndEigThing.\n");
		exit(1);
	}
	for(i = 0; i < nComps; i++) {
		mp->bigQAndEigThing[i] = (p4_bigQAndEig **)malloc(nRMatrices * sizeof(p4_bigQAndEig *));
		if(!mp->bigQAndEigThing[i]) {
			printf("Failed to allocate memory for mp->bigQAndEigThing[i].\n");
			exit(1);
		}
		for(j = 0; j < nRMatrices; j++) {
			//mp->bigQAndEigThing[i][j] = (p4_bigQAndEig *)malloc(nGdasrvs * sizeof(p4_bigQAndEig));
			mp->bigQAndEigThing[i][j] = (p4_bigQAndEig *)malloc(sizeof(p4_bigQAndEig));
			if(!mp->bigQAndEigThing[i][j]) {
				printf("Failed to allocate memory for mp->bigQAndEigThing[i][j].\n");
				exit(1);
			}
		}
	}
	for(i = 0; i < nComps; i++) {
		for(j = 0; j < nRMatrices; j++) {
			mp->bigQAndEigThing[i][j]->bigQ = NULL;
			mp->bigQAndEigThing[i][j]->qEig = NULL;
			//mp->bigQAndEigThing[i][j]->needsReset = bQETneedsReset;
		}
	}

	// bQETneedsReset is a numpy 2-D array, but it is 1-D in C, so it is indexed as [cNum * nRMatrices + rNum]
	mp->bQETneedsReset = bQETneedsReset;

	// TSCovarion
	mp->tSCov = malloc(sizeof(p4_tSCovarion));
	if(!mp->tSCov) {
		printf("Failed to alloc memory for tSCovarion.\n");
		exit(1);
	}
	mp->tSCov->free = tSCovIsFree;
	mp->tSCov->s1 = malloc(sizeof(double));
	mp->tSCov->s2 = malloc(sizeof(double));
	if(!mp->tSCov->s1 || !mp->tSCov->s1) {
		printf("Failed to allocate memory for tSCov->sN\n");
		exit(1);
	}
	//mp->tSCov->halfDim = dim / 2;
	mp->tSCov->halfComp = malloc((dim / 2) * sizeof(double));
	if(!mp->tSCov->halfComp) {
		printf("Failed to allocate memory for tSCov->halfComp\n");
		exit(1);
	}

	


	mp->freqsTimesOneMinusPInvar = NULL;
		
	//printf("finished p4_newModelPart()\n");

	// install
	aModel->parts[pNum] = mp;
}

void p4_freeModelPart(p4_modelPart *mp)
{
	int i,j;

	//printf("..p4_freeModelPart here.\n");
	// free comps
	if(mp->comps) {
		for(i = 0; i < mp->nComps; i++) {
			if(mp->comps[i]) {
				if(mp->comps[i]->val) {
					free(mp->comps[i]->val);
					mp->comps[i]->val = NULL;
				}
				free(mp->comps[i]);
				mp->comps[i] = NULL;
			}
		}
		free(mp->comps);
		mp->comps = NULL;
	}

	// free rMatrices
	if(mp->rMatrices) {
		for(i = 0; i < mp->nRMatrices; i++) {
			if(mp->rMatrices[i]) {
				if(mp->rMatrices[i]->bigR) {
					free_psdmatrix(mp->rMatrices[i]->bigR);
					mp->rMatrices[i]->bigR = NULL;
				}
				if(mp->rMatrices[i]->kappa) {
					free(mp->rMatrices[i]->kappa);
					mp->rMatrices[i]->kappa = NULL;
				}
				free(mp->rMatrices[i]);
				mp->rMatrices[i] = NULL;
			}
		}
		free(mp->rMatrices);
		mp->rMatrices = NULL;
	}

	// free gdasrvs
	if(mp->gdasrvs) {
		for(i = 0; i < mp->nGdasrvs; i++) {
			if(mp->gdasrvs[i]) {
				//free(mp->gdasrvs[i]->freqs);
				mp->gdasrvs[i]->freqs = NULL;
				//free(mp->gdasrvs[i]->rates);
				mp->gdasrvs[i]->rates = NULL;
				//free(mp->gdasrvs[i]->val);
				mp->gdasrvs[i]->val = NULL;
				free(mp->gdasrvs[i]);
				mp->gdasrvs[i] = NULL;
			}
		}
		free(mp->gdasrvs);
		mp->gdasrvs = NULL;
	}

	// free pInvar
	if(mp->pInvar) {
		if(mp->pInvar->val) {
			free(mp->pInvar->val);
			mp->pInvar->val = NULL;
		}
		free(mp->pInvar);
		mp->pInvar = NULL;
	}

	// mixture
	if(mp->mixture) {
		free(mp->mixture->freqs);
		mp->mixture->freqs = NULL;
		free(mp->mixture->rates);
		mp->mixture->rates = NULL;
		free(mp->mixture);
		mp->mixture = NULL;
	}

	// relRate
	free(mp->relRate);
	mp->relRate = NULL;


	// bigQAndEigThing
	if(mp->bigQAndEigThing) {
		for(i = 0; i < mp->nComps; i++) {
			for(j = 0; j < mp->nRMatrices; j++) {
				if(mp->bigQAndEigThing[i][j]->bigQ) {
					free_psdmatrix(mp->bigQAndEigThing[i][j]->bigQ);
					mp->bigQAndEigThing[i][j]->bigQ = NULL;
					freeEig(mp->bigQAndEigThing[i][j]->qEig);
					mp->bigQAndEigThing[i][j]->qEig = NULL;
				}
				free(mp->bigQAndEigThing[i][j]);
				mp->bigQAndEigThing[i][j] = NULL;
			}
			free(mp->bigQAndEigThing[i]);
			mp->bigQAndEigThing[i] = NULL;
		}
		free(mp->bigQAndEigThing);
		mp->bigQAndEigThing = NULL;
	}

	if(mp->bQETneedsReset) {
		mp->bQETneedsReset = NULL;
	}

	// tSCovarion
	if(mp->doTSCovarion) {
		free(mp->tSCov->s1);
		mp->tSCov->s1 = NULL;
		free(mp->tSCov->s2);
		mp->tSCov->s2 = NULL;
		free(mp->tSCov->halfComp);
		mp->tSCov->halfComp = NULL;

		free(mp->tSCov);
		mp->tSCov = NULL;
	}

	if(mp->freqsTimesOneMinusPInvar) {
		free(mp->freqsTimesOneMinusPInvar);
	}
	mp->freqsTimesOneMinusPInvar = NULL;


	free(mp);
	mp = NULL;
	//printf("..freed p4_modelPart\n");
}



void p4_newComp(p4_model *aModel, int pNum, int mNum, int free)
{
	p4_comp  *aComp;

	//printf("....p4_newComp here.  pNum=%i, mNum=%i\n", pNum, mNum);
	aComp = malloc(sizeof(p4_comp));
	if(!aComp) {
		printf("Failed to alloc memory for comp\n");
		exit(1);
	}
	aComp->free = free;
	aComp->val = malloc(aModel->parts[pNum]->dim * sizeof(double));
	if(!aComp->val) {
		printf("Failed to alloc memory for comp->val\n");
		exit(1);
	}
	aModel->parts[pNum]->comps[mNum] = aComp;
	
}


void p4_newRMatrix(p4_model *aModel, int pNum, int mNum, int free, int spec)
{
	p4_rMatrix  *aRMatrix;
	int i, j;
	int dim;

	//printf("....p4_newRMatrix here.  pNum=%i, mNum=%i\n", pNum, mNum);
	aRMatrix = malloc(sizeof(p4_rMatrix));
	if(!aRMatrix) {
		printf("Failed to alloc memory for rMatrix\n");
		exit(1);
	}
	aRMatrix->free = free;
	aRMatrix->spec = spec;
	aRMatrix->kappa = NULL;

	// malloc bigR
	dim = aModel->parts[pNum]->dim;
	aRMatrix->bigR = psdmatrix(dim);

	// fill the bigR
	if(spec == RMATRIX_CPREV) {
	  cpREVRMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_D78) {
	  d78RMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_JTT) {
	  jttRMatrix(aRMatrix->bigR);
	}
	else if (spec == RMATRIX_MTREV24) {
	  mtREV24RMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_MTMAM) {
	  mtmamRMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_WAG) {
	  wagRMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_RTREV) {
	  rtRevRMatrix(aRMatrix->bigR);
	}

	else if(spec == RMATRIX_TMJTT94) {
	  tmjtt94RMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_TMLG99) {
	  tmlg99RMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_LG) {
	  lgRMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_BLOSUM62_A) {
	  blosum62RMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_HIVB) {
	  hivbRMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_MTART) {
	  mtartRMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_MTZOA) {
	  mtzoaRMatrix(aRMatrix->bigR);
	}
	else if(spec == RMATRIX_GCPREV) {
	  gcpREVRMatrix(aRMatrix->bigR);
	}
	//else if(spec == RMATRIX_BLOSUM62_B) {
	//	printf("blosum62b is not ready yet.\n");
	//	exit(1);
	//}
	//else if(spec == RMATRIX_PHAT70) {
	//	printf("phat70 is not ready yet.\n");
	//	exit(1);
	//}
	else if(spec == RMATRIX_2P) {
		aRMatrix->kappa = malloc(sizeof(double *));
		if(!aRMatrix->kappa) {
			printf("Failed to alloc memory for kappa.\n");
			exit(1);
		}
		aRMatrix->kappa[0] = 2.0;
	}
	else {
		// poisson etc, all ones, the default.  If it is a specified
		// matrix, ie RMATRIX_SPECIFIED, then values will be poked
		// into aRMatrix->bigR later.
		for(i = 0; i < dim; i++) {
			for(j = 0; j < dim; j++) {
				aRMatrix->bigR[i][j] = 1.0;
			}
		}
	}
	
	// install
	aModel->parts[pNum]->rMatrices[mNum] = aRMatrix;
	//printf("allocated rMatrix=%li, rMatrix->bigR=%li\n", (long int)aRMatrix, (long int)(aRMatrix->bigR));
}


#if 0

p4_gdasrv *p4_newGdasrv(p4_model *aModel, int pNum, int mNum, int nCat, int free) //, PyArrayObject *val, PyArrayObject *freqs, PyArrayObject *rates)
{
	p4_gdasrv  *aGdasrv;

	//printf("=p4_newGdasrv here.  pNum=%i, mNum=%i", pNum, mNum);
	aGdasrv = malloc(sizeof(p4_gdasrv));
	if(!aGdasrv) {
		printf("Failed to alloc memory for gdasrv\n");
		exit(1);
	}
	//printf(", new aGdasrv=%li\n", (long int)aGdasrv);
	aGdasrv->free = free;

	aGdasrv->val = malloc(sizeof(double));
	if(!aGdasrv->val) {
		printf("Failed to alloc memory for gdasrv->val\n");
		exit(1);
	}
	aGdasrv->val[0] = -1.0;

	aGdasrv->freqs = malloc(aModel->parts[pNum]->nCat * sizeof(double));
	aGdasrv->rates = malloc(aModel->parts[pNum]->nCat * sizeof(double));
	if(!aGdasrv->freqs || !aGdasrv->rates) {
		printf("Failed to alloc memory for gdasrv->freqs/rates\n");
		exit(1);
	}

	aModel->parts[pNum]->gdasrvs[mNum] = aGdasrv;
	return aGdasrv;
}
#endif


p4_gdasrv *p4_newGdasrv(p4_model *aModel, int pNum, int mNum, int nCat, int free, PyArrayObject *val, PyArrayObject *freqs, PyArrayObject *rates)
{
	p4_gdasrv  *aGdasrv;

	//printf("=p4_newGdasrv here.  pNum=%i, mNum=%i", pNum, mNum);
	aGdasrv = malloc(sizeof(p4_gdasrv));
	if(!aGdasrv) {
		printf("Failed to alloc memory for gdasrv\n");
		exit(1);
	}
	//printf(", new aGdasrv=%li\n", (long int)aGdasrv);
	aGdasrv->free = free;
	aGdasrv->nCat = nCat;

	aGdasrv->val = (double *)val->data;
	aGdasrv->freqs = (double *)freqs->data;
	aGdasrv->rates = (double *)rates->data;

	aModel->parts[pNum]->gdasrvs[mNum] = aGdasrv;
	//printf("p4_newGdasrv returning %li\n", (long int)aGdasrv);
	return aGdasrv;
}

void p4_resetBQET(p4_model *aModel, int pNum, int compNum, int rMatrixNum)
{
	int ret;
	p4_modelPart   *mp;
	p4_bigQAndEig  *aQE;
	p4_rMatrix     *r;
	p4_comp        *c;

	mp = aModel->parts[pNum];
	aQE = mp->bigQAndEigThing[compNum][rMatrixNum];
	c = mp->comps[compNum];
	r = mp->rMatrices[rMatrixNum];

	if(!aQE->bigQ) {
		aQE->bigQ = psdmatrix(mp->dim);
		setBigQFromRMatrixDotCharFreq(aQE->bigQ, r->bigR, c->val, mp->dim);
		normalizeBigQ(aQE->bigQ, c->val, mp->dim);
		aQE->qEig = allocEig(mp->dim, aQE->bigQ);
	}
	else {
		setBigQFromRMatrixDotCharFreq(aQE->bigQ, r->bigR, c->val, mp->dim);
		normalizeBigQ(aQE->bigQ, c->val, mp->dim);
	}
	ret = eigensystem(aQE->qEig);
	if(ret) {
		printf("p4_resetBQET()  There is a problem with the eigensystem.\n");
		exit(1);
	}
	mp->bQETneedsReset[(compNum * mp->nRMatrices) + rMatrixNum] = 0;
	
}
