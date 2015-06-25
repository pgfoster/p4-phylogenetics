#include "pftypes.h"
#include "p4_node.h"
#include "p4_model.h"
#include "pmatrices.h"
#include "eig.h"
#include "defines.h"
#include <math.h> //For fabs


p4_node *p4_newNode(int nodeNum, p4_tree *aTree, int seqNum, int isLeaf, int inTree)
{
	p4_node  *aNode;
	int i, j;

	aNode = malloc(sizeof(p4_node));
	if(!aNode) {
		printf("Failed to allocate memory for p4_node.\n");
		exit(1);
	}
	aNode->nodeNum = nodeNum;
	//printf("p4_newNode %i\n", aNode->nodeNum);
	aNode->tree = aTree;
	aNode->parent = NULL;
	aNode->leftChild = NULL;
	aNode->sibling = NULL;
	aNode->seqNum = seqNum;
	aNode->isLeaf = isLeaf;
	aNode->nParts = aTree->nParts;

	aNode->brLen = malloc(sizeof(double));
	if(!aNode->brLen) {
		printf("Failed to allocate memory for brLen\n");
		exit(1);
	}
	aNode->brLen[0] = -1.0;

	// compNums
	aNode->compNums = malloc(aNode->nParts * sizeof(int));
	if(!aNode->compNums) {
		printf("Failed to allocate memory for compNums\n");
		exit(0);
	}
	for(i = 0; i < aNode->nParts; i++) {
		aNode->compNums[i] = 0;
	}

	// rMatrixNums
	aNode->rMatrixNums = malloc(aNode->nParts * sizeof(int));
	if(!aNode->rMatrixNums) {
		printf("Failed to allocate memory for rMatrixNums\n");
		exit(1);
	}
	for(i = 0; i < aNode->nParts; i++) {
		aNode->rMatrixNums[i] = 0;
	}

	// gdasrvNums
	aNode->gdasrvNums = malloc(aNode->nParts * sizeof(int));
	if(!aNode->gdasrvNums) {
		printf("Failed to allocate memory for gdasrvNums\n");
		exit(0);
	}
	for(i = 0; i < aNode->nParts; i++) {
		aNode->gdasrvNums[i] = 0;
	}

		
	// bigPDecks
	aNode->bigPDecks = (double ****)malloc(aNode->nParts * sizeof(double ***));
	if(!aNode->bigPDecks) {
		printf("Failed to allocate memory for bigPDecks.\n");
		exit(1);
	}
	for(i = 0; i < aNode->nParts; i++) {
		aNode->bigPDecks[i] = (double ***)malloc(aNode->tree->model->parts[i]->nCat * sizeof(double **));
		if(!aNode->bigPDecks[i]) {
			printf("Failed to allocate memory for bigPDecks[i].\n");
			exit(1);
		}
		for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
			aNode->bigPDecks[i][j] = psdmatrix(aNode->tree->model->parts[i]->dim);
		}
	}


	aNode->bigPDecks_1stD = NULL;
	aNode->bigPDecks_2ndD = NULL;


	// cl  conditionalLikelihoods
	if(!aNode->isLeaf) {
		aNode->cl = (double ****)malloc(aNode->nParts * sizeof(double ***));
		if(!aNode->cl) {
			printf("Failed to allocate memory for cl.\n");
			exit(1);
		}
		for(i = 0; i < aNode->nParts; i++) {
			aNode->cl[i] = (double ***)malloc(aNode->tree->model->parts[i]->nCat * sizeof(double **));
			if(!aNode->cl[i]) {
				printf("Failed to allocate memory for cl[i].\n");
				exit(1);
			}
			for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
				aNode->cl[i][j] = pdmatrix(aNode->tree->model->parts[i]->dim, aNode->tree->data->parts[i]->nChar);
			}
		}

	} else {
		aNode->cl = NULL;
	}

	aNode->cl2 = NULL;
	aNode->pickerDecks = NULL;
	if(aNode->isLeaf) {
	  aNode->clNeedsUpdating = 0; // init
	} else {
	  aNode->clNeedsUpdating = 1;  // init
		if(!inTree) {
			aNode->clNeedsUpdating = 0;  // the node is not part of the tree
		}
	}
	aNode->cl2NeedsUpdating = 1;

	aNode->expectedComp = NULL;
	aNode->observedComp = NULL;

	// install
	aTree->nodes[nodeNum] = aNode;
	//printf("finished alloc bigQ\n");

	//if(aNode->nodeNum == 0) {
	//	printf("newNode node 0.\n");
	//}

	return aNode;
}

void p4_freeNode(p4_node *aNode)
{
	int i, j;

	//printf("    x1 freeing node %i.\n", aNode->nodeNum);
	//if(aNode->nodeNum == 0) {
	//	printf("x1 freeing node %i. ", aNode->nodeNum);
	//	printf(" aNode->tree->model = %li  ", (long int)aNode->tree->model);
	//}

	if(aNode->brLen) {
		free(aNode->brLen);
		aNode->brLen = NULL;
	}
	if(aNode->compNums) {
		free(aNode->compNums);
		aNode->compNums = NULL;
	}
	if(aNode->rMatrixNums) {
		free(aNode->rMatrixNums);
		aNode->rMatrixNums = NULL;
	}
	if(aNode->gdasrvNums) {
		free(aNode->gdasrvNums);
		aNode->gdasrvNums = NULL;
	}
	//printf("    x2 freeing node %i.\n", aNode->nodeNum);

	// bigPDecks
	if(aNode->bigPDecks) {
		//printf("    x2.1 nParts %i.\n", aNode->nParts);
		//printf("    x2.2 tree %li.\n", (long int)aNode->tree);
		//printf("    x2.3 model %li.\n", (long int)aNode->tree->model);
		//printf("    x2.4 parts %li.\n", (long int)aNode->tree->model->parts);
		//printf("    x2.5 nCat %i.\n", aNode->tree->model->parts[0]->nCat);

		for(i = 0; i < aNode->nParts; i++) {
			for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
				free_psdmatrix(aNode->bigPDecks[i][j]);
				aNode->bigPDecks[i][j] = NULL;
			}
			free(aNode->bigPDecks[i]);
			aNode->bigPDecks[i] = NULL;
		}
		free(aNode->bigPDecks);
		aNode->bigPDecks = NULL;
	}

	//printf("    x3 freeing node %i.\n", aNode->nodeNum);

	// bigPDecks_1stD
	if(aNode->bigPDecks_1stD) {
		for(i = 0; i < aNode->nParts; i++) {
			for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
				free_psdmatrix(aNode->bigPDecks_1stD[i][j]);
				aNode->bigPDecks_1stD[i][j] = NULL;
			}
			free(aNode->bigPDecks_1stD[i]);
			aNode->bigPDecks_1stD[i] = NULL;
		}
		free(aNode->bigPDecks_1stD);
		aNode->bigPDecks_1stD = NULL;
	}

	// bigPDecks_2ndD
	if(aNode->bigPDecks_2ndD) {
		for(i = 0; i < aNode->nParts; i++) {
			for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
				free_psdmatrix(aNode->bigPDecks_2ndD[i][j]);
				aNode->bigPDecks_2ndD[i][j] = NULL;
			}
			free(aNode->bigPDecks_2ndD[i]);
			aNode->bigPDecks_2ndD[i] = NULL;
		}
		free(aNode->bigPDecks_2ndD);
		aNode->bigPDecks_2ndD = NULL;
	}

	// cl
	if(aNode->cl) {
		for(i = 0; i < aNode->nParts; i++) {
			for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
				free_pdmatrix(aNode->cl[i][j]);
				aNode->cl[i][j] = NULL;
			}
			free(aNode->cl[i]);
			aNode->cl[i] = NULL;
		}
		free(aNode->cl);
		aNode->cl = NULL;
	}
		
	// cl2
	if(aNode->cl2) {
		for(i = 0; i < aNode->nParts; i++) {
			for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
				free_pdmatrix(aNode->cl2[i][j]);
				aNode->cl2[i][j] = NULL;
			}
			free(aNode->cl2[i]);
			aNode->cl2[i] = NULL;
		}
		free(aNode->cl2);
		aNode->cl2 = NULL;
	}

	// pickerDecks, malloc'd in p4_calculatePickerDecks()
	if(aNode->pickerDecks) {
		for(i = 0; i < aNode->nParts; i++) {
			for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
				free_psdmatrix(aNode->pickerDecks[i][j]);
				aNode->pickerDecks[i][j] = NULL;
			}
			free(aNode->pickerDecks[i]);
			aNode->pickerDecks[i] = NULL;
		}
		free(aNode->pickerDecks);
		aNode->pickerDecks = NULL;
	}


	if(aNode->expectedComp) {  // malloc'd in p4_calculateExpectedComp, below
		for(i = 0; i < aNode->nParts; i++) {
			free_pdmatrix(aNode->expectedComp[i]);
			aNode->expectedComp[i] = NULL;
		}
		free(aNode->expectedComp);
		aNode->expectedComp = NULL;
	}

	if(aNode->observedComp) {
		for(i = 0; i < aNode->nParts; i++) {
			free(aNode->observedComp[i]);
			aNode->observedComp[i] = NULL;
		}
		free(aNode->observedComp);
		aNode->observedComp = NULL;
	}
	
	//if(aNode->nodeNum == 0) {
	//	printf("...done\n");
	//}
	//printf("    p4_node.   finished freeing node %i\n", aNode->nodeNum);
	free(aNode);
	aNode = NULL;
}


void p4_calculateBigPDecks(p4_node *aNode)
{
	int pNum;
	
	for(pNum = 0; pNum < aNode->tree->model->nParts; pNum++) {
		p4_calculateBigPDecksPart(aNode, pNum);
	}
}


void p4_calculateBigPDecksPart(p4_node *aNode, int pNum)
{
	int cNum, rNum, rate;
	p4_modelPart  *mp;
	p4_bigQAndEig *aQE;
	p4_gdasrv     *aGdasrv=NULL;
	
	mp = aNode->tree->model->parts[pNum];
	if(mp->isMixture) {
		rate = 0;
		for(cNum = 0; cNum < mp->nComps; cNum++) {
			for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
				aQE = aNode->tree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];
				if(mp->pInvar->val[0] == 0.0) {
					matrixExpTimesBranchLength(aQE->qEig,
											   (aNode->brLen[0] * mp->mixture->rates[rate] * mp->relRate[0]),
											   aNode->bigPDecks[pNum][rate]);
				}
				else {
					matrixExpTimesBranchLength(aQE->qEig,
											   (aNode->brLen[0] * mp->mixture->rates[rate] * mp->relRate[0]) / 
											   (1.0 - mp->pInvar->val[0]),
											   aNode->bigPDecks[pNum][rate]);
				}
				rate++;
			}
		}
	}
	else {
		cNum = aNode->compNums[pNum];
		rNum = aNode->rMatrixNums[pNum];
		//printf("p4_calculateBigPDecksPart()  b  pNum=%i, cNum=%i, rNum=%i\n", pNum, cNum, rNum);
		aQE = mp->bigQAndEigThing[cNum][rNum];

		// This should not happen, but it might.
		if(mp->bQETneedsReset[(cNum * mp->nRMatrices) + rNum]) {
			printf("p4_calculateBigPDecksPart() pNum=%i, compNum=%i, rMatrixNum=%i, needsReset. Fix me.\n",
				   pNum, cNum, rNum);
			p4_resetBQET(aNode->tree->model, pNum, cNum, rNum);
			//exit(1);
		}
			
		//printf("p4_calculateBigPDecksPart()  c\n");
		if(mp->nGdasrvs) {
			aGdasrv = mp->gdasrvs[aNode->gdasrvNums[pNum]];
		}
		for(rate = 0; rate < mp->nCat; rate++) {
			if(mp->pInvar->val[0] == 0.0) {
				if(mp->nGdasrvs) {
					matrixExpTimesBranchLength(aQE->qEig,
											   (aNode->brLen[0] * aGdasrv->rates[rate] * mp->relRate[0]),
											   aNode->bigPDecks[pNum][rate]);
				} else { // no gdasrvs, nCat must be 1
					matrixExpTimesBranchLength(aQE->qEig,
											   (aNode->brLen[0] * mp->relRate[0]),
											   aNode->bigPDecks[pNum][rate]);
				}
				//dump_psdmatrix(aNode->bigPDecks[pNum][rate], 4);
				
			} else { // so pInvar is more than zero
				if(mp->nGdasrvs) {
					matrixExpTimesBranchLength(aQE->qEig,
											   (aNode->brLen[0] * aGdasrv->rates[rate] * mp->relRate[0]) / 
											   (1.0 - mp->pInvar->val[0]),
											   aNode->bigPDecks[pNum][rate]);
				} else { // no gdasrvs, nCat must be 1
					matrixExpTimesBranchLength(aQE->qEig,
											   (aNode->brLen[0] * mp->relRate[0]) / (1.0 - mp->pInvar->val[0]),
											   aNode->bigPDecks[pNum][rate]);
				}
			}
		}
	}
}


void p4_calculateBigPDecks_1stD(p4_node *aNode)
{
	int pNum, cNum, rNum, rate;
	p4_modelPart  *mp;
	p4_bigQAndEig *aQE;
	p4_gdasrv     *aGdasrv=NULL;
	double  temp;
	
	for(pNum = 0; pNum < aNode->tree->model->nParts; pNum++) {
		mp = aNode->tree->model->parts[pNum];
		if(mp->isMixture) {
			rate = 0;
			for(cNum = 0; cNum < mp->nComps; cNum++) {
				for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
					aQE = aNode->tree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];
					if(mp->pInvar->val[0] == 0.0) {
						firstDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
												   (aNode->brLen[0] * mp->mixture->rates[rate] * mp->relRate[0]),
																	aNode->bigPDecks_1stD[pNum][rate],
																	mp->mixture->rates[rate] * mp->relRate[0]);
					}
					else {
						firstDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
												   (aNode->brLen[0] * mp->mixture->rates[rate] * mp->relRate[0]) / 
												   (1.0 - mp->pInvar->val[0]),
												   aNode->bigPDecks_1stD[pNum][rate], 
																	(mp->mixture->rates[rate] * mp->relRate[0]) / 
												   (1.0 - mp->pInvar->val[0]));
					}
					rate++;
				}
			}

		}
		else {
			cNum = aNode->compNums[pNum];
			rNum = aNode->rMatrixNums[pNum];
			aQE = aNode->tree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];				
			if(mp->nGdasrvs) {
				aGdasrv = mp->gdasrvs[aNode->gdasrvNums[pNum]];
			}
			for(rate = 0; rate < mp->nCat; rate++) {
				if(mp->pInvar->val[0] == 0.0) {
					if(mp->nGdasrvs) {
						temp = aGdasrv->rates[rate] * mp->relRate[0];
						firstDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
													  (aNode->brLen[0] * temp),
													  aNode->bigPDecks_1stD[pNum][rate], 
													  temp);
					} else { // no gdasrvs, nCat must be 1
						firstDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
																	(aNode->brLen[0] * mp->relRate[0]),
																	aNode->bigPDecks_1stD[pNum][rate], mp->relRate[0]);
					}
				
				} else { // so pInvar is more than zero
					if(mp->nGdasrvs) {
						temp = aGdasrv->rates[rate] * mp->relRate[0];
						firstDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
												 (aNode->brLen[0] * temp) / (1.0 - mp->pInvar->val[0]),
												 aNode->bigPDecks_1stD[pNum][rate],
												 temp /(1.0 - mp->pInvar->val[0]));
					} else { // no gdasrvs, nCat must be 1
						firstDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
												(aNode->brLen[0] * mp->relRate[0]) / (1.0 - mp->pInvar->val[0]),
												aNode->bigPDecks_1stD[pNum][rate],
												mp->relRate[0] / (1.0 - mp->pInvar->val[0]));
					}
				}
			}
		}
	}
}



void p4_calculateBigPDecks_2ndD(p4_node *aNode)
{
	int pNum, cNum, rNum, rate;
	p4_modelPart  *mp;
	p4_bigQAndEig *aQE;
	p4_gdasrv     *aGdasrv=NULL;
	double  temp;
	
	for(pNum = 0; pNum < aNode->tree->model->nParts; pNum++) {
		mp = aNode->tree->model->parts[pNum];
		if(mp->isMixture) {
			rate = 0;
			for(cNum = 0; cNum < mp->nComps; cNum++) {
				for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
					aQE = aNode->tree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];
					if(mp->pInvar->val[0] == 0.0) {
						secondDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
												   (aNode->brLen[0] * mp->mixture->rates[rate] * mp->relRate[0]),
																	aNode->bigPDecks_2ndD[pNum][rate],
																	mp->mixture->rates[rate] * mp->relRate[0]);
					}
					else {
						secondDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
												   (aNode->brLen[0] * mp->mixture->rates[rate] * mp->relRate[0]) / 
												   (1.0 - mp->pInvar->val[0]),
												   aNode->bigPDecks_2ndD[pNum][rate], 
																	(mp->mixture->rates[rate] * mp->relRate[0]) / 
												   (1.0 - mp->pInvar->val[0]));
					}
					rate++;
				}
			}

		}
		else {
			cNum = aNode->compNums[pNum];
			rNum = aNode->rMatrixNums[pNum];
			aQE = aNode->tree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];				
			if(mp->nGdasrvs) {
				aGdasrv = mp->gdasrvs[aNode->gdasrvNums[pNum]];
			}
			for(rate = 0; rate < mp->nCat; rate++) {
				if(mp->pInvar->val[0] == 0.0) {
					if(mp->nGdasrvs) {
						temp = aGdasrv->rates[rate] * mp->relRate[0];
						secondDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
												 (aNode->brLen[0] * temp),
																	 aNode->bigPDecks_2ndD[pNum][rate], 
																	 temp * mp->relRate[0]);
					} else { // no gdasrvs, nCat must be 1
						secondDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
																	 (aNode->brLen[0] * mp->relRate[0]),
																	 aNode->bigPDecks_2ndD[pNum][rate], mp->relRate[0]);
					}
				
				} else { // so pInvar is more than zero
					if(mp->nGdasrvs) {
						temp = aGdasrv->rates[rate] * mp->relRate[0];
						secondDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
																	 (aNode->brLen[0] * temp) / (1.0 - mp->pInvar->val[0]),
																	 aNode->bigPDecks_2ndD[pNum][rate],
																	 temp /(1.0 - mp->pInvar->val[0]));
					} else { // no gdasrvs, nCat must be 1
						secondDerivativeOfMatrixExpTimesBranchLength(aQE->qEig,
													 (aNode->brLen[0] * mp->relRate[0]) / (1.0 - mp->pInvar->val[0]),
																	 aNode->bigPDecks_2ndD[pNum][rate],
																	 mp->relRate[0] / (1.0 - mp->pInvar->val[0]));
					}
				}
			}
		}
	}
}




void p4_calculatePickerDecks(p4_node *aNode)
{
	int i, j, k, l;
	p4_modelPart *mp;
	
	//printf("Node %i calculatePickerDecks\n", aNode->nodeNum);

	if(!aNode->pickerDecks) {
		aNode->pickerDecks = (double ****)malloc(sizeof(double ***) * aNode->nParts);
		if(!aNode->pickerDecks) {
			printf("Failed to allocate memory for pickerDecks.\n");
			exit(0);
		}
		for(i = 0; i < aNode->nParts; i++) {
			mp = aNode->tree->model->parts[i];
			//printf("nCat = %i, dim = %i\n", mp->nCat, mp->dim);
			aNode->pickerDecks[i] = (double ***)malloc(mp->nCat * sizeof(double **));
			if(!aNode->pickerDecks[i]) {
				printf("Failed to allocate memory for pickerDecks[i].\n");
				exit(0);
			}
			for(j = 0; j < mp->nCat; j++) {
				aNode->pickerDecks[i][j] = psdmatrix(mp->dim);
			}
		}
	}


	// fill
	for(i = 0; i < aNode->nParts; i++) {
		mp = aNode->tree->model->parts[i];
		for(j = 0; j < mp->nCat; j++) {
			for(k = 0; k < mp->dim; k++) {
				aNode->pickerDecks[i][j][k][0] = aNode->bigPDecks[i][j][k][0];
				for(l = 1; l < (mp->dim - 1); l++) {
					aNode->pickerDecks[i][j][k][l] = aNode->pickerDecks[i][j][k][l - 1] + 
						aNode->bigPDecks[i][j][k][l];
				}
				aNode->pickerDecks[i][j][k][mp->dim - 1] = 1.0;
			}
		}
	}

#if 0
	if(aNode->nodeNum == 1) {
		printf("pickerdeck for node 1\n");
		for(i = 0; i < 4; i++) {
			for(j = 0; j < 4; j++) {
				printf(" %f", aNode->pickerDecks[0][0][i][j]);
			}
			printf("\n");
		}
	}
#endif
}

void p4_setConditionalLikelihoodsOfInteriorNode(p4_node *aNode)
{
	int pNum;

	for(pNum = 0; pNum < aNode->nParts; pNum++) {	
		p4_setConditionalLikelihoodsOfInteriorNodePart(aNode, pNum);
	}
	aNode->clNeedsUpdating = 0;  // for treeNewt, if nothing else.

}

void p4_setConditionalLikelihoodsOfInteriorNodePart(p4_node *aNode, int pNum)
{
	int	seqPos, symb, chSymb;
	p4_node *child;
	double	sum = 0.0;
	int rate;
	//int j;
	int charCode;
	int isN = 0;
	int dim;
	part   *dp;         // a data part
	p4_modelPart  *mp;  // a modelPart
	
	dp = aNode->tree->data->parts[pNum];
	mp = aNode->tree->model->parts[pNum];
	dim = mp->dim;
	for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {			
		for(rate = 0; rate < mp->nCat; rate++){
			for(symb = 0; symb < dim; symb++) {
				//printf("part %i, seqPos %i, rate %i, symb %i\n", pNum, seqPos, rate, symb);
				//printf("symbolPos %i\n", symb);
					
				// first, do the leftChild...

				if(aNode->leftChild->isLeaf) {   // If its a leaf, don't use conditionalLikelihoods of leftChild

					// The charCode is one of:
					// - a regular unambiguous character, from zero to dim-1
					// - a gap or a question mark, -1 or -2
					// - an equate, ie well less than zero.
					//     - an n, or n-like character, which would be treated like a gap
					//     - a less general equate, eg "y".

					charCode = dp->patterns[aNode->leftChild->seqNum][seqPos];
					//printf("  Got charCode = %i\n", charCode);
					if(charCode >= 0) {
#if 0
						if(charCode >= dp->dim) {
							printf("0 Got charCode = %i\n", charCode);
							printf("part=%i, seqPos=%i, rate=%i, symb=%i\n", pNum,seqPos,rate,symb);
							printf("nPatterns=%i, nCat=%i, data dim=%i, sequenceNum=%i\n", 
								   dp->nPatterns, 
								   mp->nCat, 
								   dp->dim,
								   aNode->leftChild->seqNum);
							exit(1);
						}
#endif
						if(mp->doTSCovarion) {
							aNode->cl[pNum][rate][symb][seqPos] = 
								aNode->leftChild->bigPDecks[pNum][rate][symb][charCode] + 
								aNode->leftChild->bigPDecks[pNum][rate][symb][charCode + dp->dim];
						}
						else {
							//printf("  Setting cl[%i] to %f\n", 
							//symb, aNode->leftChild->bigPDecks[pNum][rate][symb][charCode]);
							aNode->cl[pNum][rate][symb][seqPos] = 
								aNode->leftChild->bigPDecks[pNum][rate][symb][charCode];
						}
					}
					else if (charCode == GAP_CODE || charCode == QMARK_CODE) {
						aNode->cl[pNum][rate][symb][seqPos] = 1.0;
					}
					else if(charCode >= EQUATES_BASE && 
							charCode < EQUATES_BASE + dp->nEquates) {
						// If we get this to this point, most
						// usually, it will be an "n" or n-like
						// character (eg x).  Assume it is so unless
						// datapart->equates says otherwise.
						isN = 1;
						for(chSymb = 0; chSymb < dp->dim; chSymb++) {
							if(!dp->equates[charCode - EQUATES_BASE][chSymb]) {
								isN = 0;
								break;
							}
						}
						//printf("isN = %i\n", isN);
						if(isN) {
							aNode->cl[pNum][rate][symb][seqPos] = 1.0;
						}
						else {
							sum = 0.0;
							for(chSymb = 0; chSymb < dp->dim; chSymb++) {
								if(dp->equates[charCode - EQUATES_BASE][chSymb]) {
									sum += aNode->leftChild->bigPDecks[pNum][rate][symb][chSymb];
								}
							}
							aNode->cl[pNum][rate][symb][seqPos] = sum;
						}
					}

					else {
						printf("node %i: setConditionalLikelihoodsOfInteriorNodes:\n", aNode->nodeNum);
						printf("   Programming error.  This shouldn't happen\n");
						exit(0);
					}
						
				}
				else {  // The leftChild is not terminal, and so use conditionalLikelihoods
					sum = 0.0;
					for(chSymb = 0; chSymb < dim; chSymb++) {

#if 0
						printf("   %f = %f + (prob[%i]=%f * condLike[%i]=%f)\n",
							   sum + (aNode->leftChild->bigPDecks[pNum][rate][symb][chSymb] * 
									  aNode->leftChild->cl[pNum][rate][chSymb][seqPos]),
							   sum,
							   chSymb,
							   aNode->leftChild->bigPDecks[pNum][rate][symb][chSymb],
							   chSymb,
							   aNode->leftChild->cl[pNum][rate][chSymb][seqPos]);
#endif
						sum = sum + (aNode->leftChild->bigPDecks[pNum][rate][symb][chSymb] * 
									 aNode->leftChild->cl[pNum][rate][chSymb][seqPos]);
					}		
					aNode->cl[pNum][rate][symb][seqPos] = sum;
				}
				//printf("                                  after leftChild, cl[%i] = %f\n", 
				//	   symb, aNode->cl[pNum][rate][symb][seqPos]);



				// then do the rest of the children until they run out
				child = aNode->leftChild->sibling;
				while(child != NULL) {
					if(child->isLeaf) {   // don't use conditionalLikelihoods of child
						//if(0) {
						charCode = dp->patterns[child->seqNum][seqPos];
						//printf("leaf sibs Got charCode = %i\n", charCode);
						if(charCode >= 0) {
							if(mp->doTSCovarion) {
								aNode->cl[pNum][rate][symb][seqPos] *= (child->bigPDecks[pNum][rate][symb][charCode] + 
																		child->bigPDecks[pNum][rate][symb][charCode + dp->dim]);
							}
							else {
#if 0
								printf("node %i, adding cl for leaf node %i, charCode=%i, bigP=%f\n",
									   aNode->nodeNum,
									   child->nodeNum,
									   charCode,
									   child->bigPDecks[pNum][rate][symb][charCode]);
								printf("   cl(before)=%f   %f * %f = ", 
									   aNode->cl[pNum][rate][symb][seqPos],
									   aNode->cl[pNum][rate][symb][seqPos],
									   child->bigPDecks[pNum][rate][symb][charCode]);
#endif
								aNode->cl[pNum][rate][symb][seqPos] *= child->bigPDecks[pNum][rate][symb][charCode];
								//printf("%f\n", aNode->cl[pNum][rate][symb][seqPos]);
							}
						}
						else if (charCode == GAP_CODE || charCode == QMARK_CODE) {
							//aNode->cl[pNum][rate][symb][seqPos] *= 1.0;
						}
						else if(charCode >= EQUATES_BASE && 
								charCode < EQUATES_BASE + dp->nEquates) {
							isN = 1;
							for(chSymb = 0; chSymb < dp->dim; chSymb++) {
								if(!dp->equates[charCode - EQUATES_BASE][chSymb]) {
									isN = 0;
									break;
								}
							}
							if(isN) {
								//printf("x isN = %i\n", isN);
								//aNode->cl[pNum][rate][symb][seqPos] = 1.0;
							}
							else {
								//printf("x isN = %i\n", isN);
								sum = 0.0;
								for(chSymb = 0; chSymb < dp->dim; chSymb++) {
									if(dp->equates[charCode - EQUATES_BASE][chSymb]) {
										sum = sum + child->bigPDecks[pNum][rate][symb][chSymb];
									}
								}
								aNode->cl[pNum][rate][symb][seqPos] *= sum;
							}
						}

						else {
							printf("node %i: setConditionalLikelihoodsOfInteriorNodes:\n", aNode->nodeNum);
							printf("x   Programming error.  This shouldn't happen\n");
							exit(0);
						}
						//printf("                                  after sib, cl[%i] = %f\n", 
						//	   symb, aNode->cl[pNum][rate][symb][seqPos]);
							
						
					}
					else {  // child is not terminal, so use conditionalLikelihoods of child
						sum = 0.0;
						for(chSymb = 0; chSymb < dim; chSymb++) {

#if 0
							printf("   %f = %f + (prob=%f * cl=%f)\n",
								   sum + (child->bigPDecks[pNum][rate][symb][chSymb] * 
										  child->cl[pNum][rate][chSymb][seqPos]),
								   sum,
								   child->bigPDecks[pNum][rate][symb][chSymb],
								   child->cl[pNum][rate][chSymb][seqPos]);
#endif

							sum = sum + (child->bigPDecks[pNum][rate][symb][chSymb] * 
										 child->cl[pNum][rate][chSymb][seqPos]);
						}		
#if 0					
						printf("                                after sibling, cl[%i] = (%f * %f) = %f\n",
							   symb,
							   aNode->cl[pNum][rate][symb][seqPos],
							   sum,
							   aNode->cl[pNum][rate][symb][seqPos] * sum);
#endif
						aNode->cl[pNum][rate][symb][seqPos] *= sum;
					}
					child = child->sibling;
				} // while(child != NULL)
			} // for symb
		} // for rate
	}
#if 0
	printf("   node %i condLikes: \n", aNode->nodeNum);
	for(seqPos = 0; seqPos < 1; seqPos++) {			
		for(rate = 0; rate < 2; rate++){
			printf("seqPos %i, rate %i: ", seqPos, rate);
			for(symb = 0; symb < mp->dim; symb++) {  // if covarion, full dim
				printf(" %6g", aNode->cl[pNum][rate][symb][seqPos]);
			}
			printf("\n");
		}
	}
#endif		

}


void p4_initializeCL2ToRootComp(p4_node *aNode)
{
	int	pNum, seqPos, symb, rate, dim, rootCompNum;
	part   *dp;         // a data part
	p4_modelPart  *mp;  // a modelPart
	int cNum, rNum;

	//printf("p4_initializeCL2(), node %i\n", aNode->nodeNum);
	for(pNum = 0; pNum < aNode->nParts; pNum++) {
		dp = aNode->tree->data->parts[pNum];
		mp = aNode->tree->model->parts[pNum];
		dim = mp->dim; // in case it is covarion, take dim from model, not the data
		if(mp->isMixture) {
			for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
				rate = 0;
				for(cNum = 0; cNum < mp->nComps; cNum++) {
					for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
						for(symb = 0; symb < dim; symb++) {
							aNode->cl2[pNum][rate][symb][seqPos] = mp->comps[cNum]->val[symb];
						}
						rate++;
					}
				}
			}
		}
		else {
			rootCompNum = aNode->tree->root->compNums[pNum];
			for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {			
				for(rate = 0; rate < mp->nCat; rate++){
					for(symb = 0; symb < dim; symb++) {
						aNode->cl2[pNum][rate][symb][seqPos] = mp->comps[rootCompNum]->val[symb];
					}
				}
			}
		}
	}
}


void p4_setCL2Up(p4_node *cl2Node)
{
	p4_node  *aNode = cl2Node->parent;
	int	pNum, seqPos, symb, rate, dim, from;
	//int charCode, isN;
	double sum;
	part   *dp;         // a data part
	p4_modelPart  *mp;  // a modelPart

	//printf("p4_setCL2Up(), node %i\n", aNode->nodeNum);

	/*
                                    +--------4:A
                           +--------cl2Node
                   +-------aNode    +--------5:B
                   |       |
        x:0--------1       +--------6:C
                   |
        .          +-------7:D

		We are doing cl2Node->cl2.  Node aNode->cl2 is already done,
		getting info from branches on 1 and 7.  Now we want the info
		to come from aNode branch.  So its the aNode->cl2 convolved
		with the aNode->bigP to make the cl2Node->cl2.

	*/

	//printf("node %i cl2\n", aNode->nodeNum);
	//dumpCl(aNode->cl2);
	//printf("node %i bigP\n", aNode->nodeNum);
	//dump_psdmatrix(aNode->bigPDecks[0][0], 4);

	if(aNode->isLeaf) {
		printf("p4_setCL2Up().  parent is a leaf.  Not implemented yet.\n");
		exit(1);
	}

	for(pNum = 0; pNum < aNode->nParts; pNum++) {
		dp = aNode->tree->data->parts[pNum];
		mp = aNode->tree->model->parts[pNum];
		dim = mp->dim;
		for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {			
			for(rate = 0; rate < mp->nCat; rate++){
				for(symb = 0; symb < dim; symb++) {
					// aNode is not aLeaf.  So use cond likes
					sum = 0.0;
					for(from = 0; from < dim; from++) {
						sum += (aNode->bigPDecks[pNum][rate][from][symb] * 
								aNode->cl2[pNum][rate][from][seqPos]);
					}
					cl2Node->cl2[pNum][rate][symb][seqPos] = sum;
					
				}
			}
		}
	}
	//printf("node %i cl2\n", cl2Node->nodeNum);
	//dumpCl(cl2Node->cl2);
	
}




void p4_setCL2Down(p4_node *cl2Node, p4_node *aNode)
{
	int	pNum, seqPos, symb, rate, dim, chSymb, charCode, isN;
	double sum;
	part   *dp;         // a data part
	p4_modelPart  *mp;  // a modelPart
	
	//if(cl2Node->nodeNum == 3) {
	//	printf("node %i cl2 before p4_setCL2Down()\n", cl2Node->nodeNum);
	//	p4_dumpCL(cl2Node->cl2);
	//}

	//printf("p4_setCL2Down(), cl2Node=%i, node %i\n", cl2Node->nodeNum, aNode->nodeNum);
	for(pNum = 0; pNum < aNode->nParts; pNum++) {
		dp = aNode->tree->data->parts[pNum];
		mp = aNode->tree->model->parts[pNum];
		dim = mp->dim;
		for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {			
			for(rate = 0; rate < mp->nCat; rate++){
				for(symb = 0; symb < dim; symb++) {
					if(aNode->isLeaf) {
						charCode = dp->patterns[aNode->seqNum][seqPos];
						if(charCode >= 0) {
							if(mp->doTSCovarion) {
								cl2Node->cl2[pNum][rate][symb][seqPos] *=
									aNode->bigPDecks[pNum][rate][symb][charCode] +
									aNode->bigPDecks[pNum][rate][symb][charCode + dp->dim];
							}
							else {
								cl2Node->cl2[pNum][rate][symb][seqPos] *= 
									aNode->bigPDecks[pNum][rate][symb][charCode];
							}
						}
						else if(charCode == GAP_CODE || charCode == QMARK_CODE) {
							//cl2Node->cl2[pNum][rate][symb][seqPos] *= 1.0;
						}
						else if(charCode >= EQUATES_BASE && charCode < EQUATES_BASE + dp->nEquates) {
							isN = 1;
							for(chSymb = 0; chSymb < dp->dim; chSymb++) {
								if(!dp->equates[charCode - EQUATES_BASE][chSymb]) {
									isN = 0;
									break;
								}
							}
							if(isN) {
								//cl2Node->cl2[pNum][rate][symb][seqPos] *= 1.0;
							}
							else {
								sum = 0.0;
								for(chSymb = 0; chSymb < dp->dim; chSymb++) {
									if(dp->equates[charCode - EQUATES_BASE][chSymb]) {
										sum += aNode->bigPDecks[pNum][rate][symb][chSymb];
									}
								}
								cl2Node->cl2[pNum][rate][symb][seqPos] *= sum;
							}
						}
						else {
							printf("p4_setCL2down().  Programming error.\n");
							exit(1);
						}
					}
					else { // aNode is not aLeaf.  So use cond likes
						sum = 0.0;
						for(chSymb = 0; chSymb < dim; chSymb++) {
							sum += (aNode->bigPDecks[pNum][rate][symb][chSymb] * 
									aNode->cl[pNum][rate][chSymb][seqPos]);
						}
						cl2Node->cl2[pNum][rate][symb][seqPos] *= sum;
					}
				}
			}
		}
	}
	//if(cl2Node->nodeNum == 3) {
	//	printf("node %i cl2\n", cl2Node->nodeNum);
	//	p4_dumpCL(cl2Node->cl2);
	//}
}


void p4_dumpCL(double ****theCL)
{
	int i, j;

	for(i = 0; i < 1; i++) {
		printf("    pos %i: ", i);
		for(j = 0; j < 4; j++) {
			printf(" %8.6f", theCL[0][0][j][i]);
		}
		printf("\n");
	}
}


//==================================================


void p4_calculateExpectedComp(p4_node *aNode)
{
	int i, j;
	int pNum, catNum;
	double factor = 0.0;
	double   **invarCompo = NULL;  // no malloc.  Rather just re-use the observedComp.
	p4_modelPart  *mp;
	int debug = 0;

	//if(aNode != aNode->tree->root && aNode->bigPDecks[0][0][0][0] < 0.001) {
	//	printf("p4_calculateExpectedComp node %i bad bigP\n", aNode->nodeNum);
	//	exit(0);
	//}
	
	// Note that the expectedComp for a given part is not a vector,
	// its a matrix, one row for each nCat.  We need to do this
	// because the sites that are in a particular rate category remain
	// in that rate category in the entire tree.  And there are nParts
	// of these expectedComp matrices.  So its
	// expectedComp[nParts][nCat][dim]

	if(aNode->expectedComp == NULL) {
		aNode->expectedComp = malloc(aNode->nParts * sizeof(double **));
		if(!aNode->expectedComp) {
			printf("Problem allocating memory for aNode->expectedComp\n");
			exit(0);
		}
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			aNode->expectedComp[pNum] = pdmatrix(mp->nCat, mp->dim);
		}
	}

	// If its the root node, the expected comp is just the model comp
	if(aNode == aNode->tree->root) {
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			for(catNum = 0; catNum < mp->nCat; catNum++) {
				for(i = 0; i < mp->dim; i++) {
					aNode->expectedComp[pNum][catNum][i] = mp->comps[aNode->compNums[pNum]]->val[i];
				}
			}
		}
		return;
	}

	// If its not the root node, we need the expected char freq 
	// of the parent, and the bigPDecks
	else {
		
		// first, initialize the entire expectedComp to zero
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			for(catNum = 0; catNum < mp->nCat; catNum++) {
				for(i = 0; i < mp->dim; i++) {
					aNode->expectedComp[pNum][catNum][i] = 0.0;
				}
			}
		}

		// We re-use aNode->pickerDecks to do the calcs.  If it does not
		// exist, we can use p4_calculatePickerDecks() to alloc.
		if(!aNode->pickerDecks) {
			p4_calculatePickerDecks(aNode);
		}
		// We need more space to hold the invarCompo.  aNode->observedComp is the right size, so use it.
		if(!aNode->observedComp) {
			aNode->observedComp = (double **)malloc(aNode->nParts * sizeof(double *));
			if(!aNode->observedComp) {
				printf("failed to malloc aNode->observedComp\n");
				exit(1);
			}
			for(pNum = 0; pNum < aNode->nParts; pNum++) {
				aNode->observedComp[pNum] = (double *)malloc(aNode->tree->model->parts[pNum]->dim * sizeof(double));
				if(!aNode->observedComp[pNum]) {
					printf("failed to malloc aNode->observedComp[pNum]\n");
					exit(1);
				}
			}
		}
		invarCompo = aNode->observedComp;			

		// If pInvar[0] > 0.0, then part of the root compo does not change 
		// over the tree.  Get it from the root.
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			if(mp->pInvar->val[0] > 0.0) { 		// mp->pInvar->val[0] is initialized to -1.0;
				for(i = 0; i < mp->dim; i++) {
					invarCompo[pNum][i] =  mp->pInvar->val[0] * 
						aNode->tree->root->expectedComp[pNum][0][i]; // Here catNum=0
				}
			}
			else {
				for(i = 0; i < mp->dim; i++) {
					invarCompo[pNum][i] = 0.0;  
				}
			}
		}
				
		// copy the bigPDeck into pickerDecks
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			for(catNum = 0; catNum < mp->nCat; catNum++) {
				for(i = 0; i < mp->dim; i++) {
					for(j = 0; j < mp->dim; j++) {
						aNode->pickerDecks[pNum][catNum][i][j] = aNode->bigPDecks[pNum][catNum][i][j];
					}
				}
			}
		}

		if(debug) {
			if(aNode->nodeNum == 7) {
				printf("part 0, cat 0, node %i, bigP\n",aNode->nodeNum);
				dump_psdmatrix(aNode->bigPDecks[0][0], 4);
				//dump_psdmatrix(aNode->bigPDecks[0][1], 4);
				printf("\n");
				printf("part 0, cat 0, node %i, pickerDecks\n",aNode->nodeNum);
				dump_psdmatrix(aNode->pickerDecks[0][0], 4);
				//dump_psdmatrix(aNode->pickerDecks[0][1], 4);
				printf("\n");
				printf("node %i, parent expectedComp\n",aNode->nodeNum);
				dump_pdvector(aNode->parent->expectedComp[0][0], 4);
				//dump_pdvector(aNode->parent->expectedComp[0][1], 4);
				printf("\n");
				printf("node %i, invarCompo[pNum=0]\n", aNode->nodeNum);
				dump_pdvector(invarCompo[0], 4);
			}
			printf("\n");
		}
		
		// multiply the rows by the parent expectedComp of variable sites
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			for(catNum = 0; catNum < mp->nCat; catNum++) {
				for(i = 0; i < mp->dim; i++) {
					for(j = 0; j < mp->dim; j++) {
						aNode->pickerDecks[pNum][catNum][i][j] *= 
							aNode->parent->expectedComp[pNum][catNum][i] - invarCompo[pNum][i];
					}
				}
			}
		}

		if(debug) {
			if(aNode->nodeNum == 7) {
				printf("node %i, parent charFreq times bigP\n",aNode->nodeNum);
				dump_psdmatrix(aNode->pickerDecks[0][0], 4);
				//dump_psdmatrix(aNode->pickerDecks[0][1], 4);
			}
			printf("\n");
		}

		// add up the columns of aNode->pickerDecks
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			for(catNum = 0; catNum < mp->nCat; catNum++) {
				for(i = 0; i < mp->dim; i++) {
					for(j = 0; j < mp->dim; j++) {
						aNode->expectedComp[pNum][catNum][j] += aNode->pickerDecks[pNum][catNum][i][j];
					}
				}
			}
		}
		if(debug) {
			if(aNode->nodeNum == 7) {
				printf("node %i, sum of columns = expectedComp\n",aNode->nodeNum);
				dump_pdvector(aNode->expectedComp[0][0], 4);
				//dump_pdvector(aNode->expectedComp[0][1], 4);
			}
			printf("\n");
		}

		// Make each expectedComp sum to 1.0
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			for(catNum = 0; catNum < mp->nCat; catNum++) {
				factor = 0.0;
				for(i = 0; i < mp->dim; i++) {
					factor += aNode->expectedComp[pNum][catNum][i];
				}
				for(i = 0; i < mp->dim; i++) {
					aNode->expectedComp[pNum][catNum][i] /= factor;
				}
			}
		}

		// If there is pInvar[0], make room for it, to be added in in the next part
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			if(mp->pInvar->val[0] > 0.0) {
				factor = 1.0 - mp->pInvar->val[0];
				for(catNum = 0; catNum < mp->nCat; catNum++) {
					for(i = 0; i < mp->dim; i++) {
						aNode->expectedComp[pNum][catNum][i] *= factor;
					}
				}
			}
		}

		// Finally, add in the contribution of invariant sites
		for(pNum = 0; pNum < aNode->nParts; pNum++) {
			mp = aNode->tree->model->parts[pNum];
			if(mp->pInvar->val[0] > 0.0) {
				for(catNum = 0; catNum < mp->nCat; catNum++) {
					for(i = 0; i < mp->dim; i++) {
						aNode->expectedComp[pNum][catNum][i] += invarCompo[pNum][i];  
					}
				}
			}
		}

		if(debug) {
			if(aNode->nodeNum == 7) {
				printf("node %i, final expectedComp\n",aNode->nodeNum);
				dump_pdvector(aNode->expectedComp[0][0], 4);
			}
		}
	}
}




#ifdef NEWSTUFF


void calculateObservedCharFreq(p4_node *aNode, int verbose)
{
	int i, j;
	int counts = 0;
	
	if(!aNode->isTerminal) {
		printf("node: calculateObservedCharFreq: node %i is not a terminal node.\n", aNode->nodeNum);
		exit(1);
	}

	if(aNode->observedCharFreq == NULL) {
		aNode->observedCharFreq = malloc(aNode->nParts * sizeof(double *));
		for(i = 0; i < aNode->nParts; i++) {
			aNode->observedCharFreq[i] = malloc(aNode->models[i]->dim * sizeof(double));
		}
	}

	for(i = 0; i < aNode->nParts; i++) {
		for(j = 0; j < aNode->models[i]->dim; j++) {
			aNode->observedCharFreq[i][j] = 0.0;
		}
	}

	// tally up base counts
	for(i = 0; i < aNode->nParts; i++) {
		if(aNode->data->parts[0]->sequences) {
			//printf("node %i: calculateObservedCharFreq, using sequences, sequenceNum %i\n", 
			//	   aNode->nodeNum, aNode->sequenceNum);
			counts = 0;
			for(j = 0; j < aNode->data->parts[i]->nChar; j++) {
				if(aNode->data->parts[i]->sequences[aNode->sequenceNum][j] >= 0) {
					aNode->observedCharFreq[i][aNode->data->parts[i]->sequences[aNode->sequenceNum][j]] = 
						aNode->observedCharFreq[i][aNode->data->parts[i]->sequences[aNode->sequenceNum][j]] + 1.0;
					counts++;
				}
				else if(aNode->data->parts[i]->sequences[aNode->sequenceNum][j] == GAP_CODE ||
						aNode->data->parts[i]->sequences[aNode->sequenceNum][j] == QMARK_CODE) {
					// pass
				}
				else {
					printf("node: calculateObservedCharFreq (from sequences):\n"); 
					printf("    ambiguity code not implemented yet\n");
					exit(1);
				}
			}
			// divide by counts to get charFreq
			for(j = 0; j < aNode->data->parts[i]->dim; j++) {
				aNode->observedCharFreq[i][j] = aNode->observedCharFreq[i][j] / (double)counts;
			}
		}
		else {
			printf("node: calculateObservedCharFreq: something wrong here.\n");
			exit(1);
		}
	}
#if 0
	else if(aNode->parts[i]->patternCount) {
		//printf("node %i: calculateObservedCharFreq, using patterns, sequenceNum %i\n", 
		//					aNode->nodeNum, aNode->sequenceNum);
		counts = 0;
		for(j = 0; j < aNode->parts[i]->patternCount; j++) {
			if(aNode->parts[i]->patterns[aNode->sequenceNum][j] >= 0) {
				aNode->observedCharFreq[aNode->parts[i]->patterns[aNode->sequenceNum][j]] = 
						aNode->observedCharFreq[aNode->parts[i]->patterns[aNode->sequenceNum][j]] + 
									(double)aNode->parts[i]->patternCounts[j];
				counts = counts + aNode->parts[i]->patternCounts[j];
			}
			else if(aNode->parts[i]->patterns[aNode->sequenceNum][j] >= -4) {
				// pass, its one of -, ?, n, or x
			}
			else {
				printf("node: calculateObservedCharFreq: ambiguity code not implemented yet\n");
				exit(1);
			}
		}
		// divide by counts to get charFreq
		for(j = 0; j < aNode->models[i]->dim; j++) {
			aNode->observedCharFreq[j] = aNode->observedCharFreq[j] / (double)counts;
		}
	} 
#endif


    if(verbose) {
		printf("node %2i, sequence %i ", aNode->nodeNum, aNode->sequenceNum);
		for(i = 0; i < aNode->nParts; i++) {
			printf("  part %i: ", i);
			for(j = 0; j < aNode->data->parts[i]->dim; j++) {
				printf("%9.6f", aNode->observedCharFreq[i][j]);
			}
			printf("\n");
		}
		printf("\n");
	}
}

#endif // NEWSTUFF


////===========================================
////
////  Copy and verify cond likes and bigP
////
////===========================================

void p4_copyCondLikesFromNodeToNode(p4_node *nA, p4_node *nB)
{
	int partNum, catNum, stateNum, patNum;
	
	/*
	  cl:
	    aNode->nParts
	    aNode->model->parts[partNum]->nCat
	    aNode->model->parts[partNum]->dim, 
	    aNode->tree->data->parts[partNum]->nChar

		but for the latter, we can use nPatterns
	*/

	//printf("copyCondLikesFromNodeToNode.  nParts=%i\n", nA->nParts);
	for(partNum = 0; partNum < nA->nParts; partNum++) {
		for(catNum = 0; catNum < nA->tree->model->parts[partNum]->nCat; catNum++) {
			for(stateNum = 0; stateNum < nA->tree->model->parts[partNum]->dim; stateNum++) {
				for(patNum = 0; patNum < nA->tree->data->parts[partNum]->nPatterns; patNum++) {
					nB->cl[partNum][catNum][stateNum][patNum] = 
						nA->cl[partNum][catNum][stateNum][patNum];
				}
			}
		}
		nB->tree->data->parts[partNum]->nPatterns = nA->tree->data->parts[partNum]->nPatterns;
	}
	

}


int p4_verifyCondLikesFromNodeToNode(p4_node *nA, p4_node *nB)
{
	int partNum, catNum, stateNum, patNum;
	double epsilon, diff;

        epsilon = 1.e-15;
	
	/*
	  cl:
	    aNode->nParts
	    aNode->tree->model->parts[partNum]->nCat
	    aNode->tree->model-parts[partNum]->dim, 
	    aNode->tree->data->parts[partNum]->nChar

		but for the latter, we can use nPatterns
	*/

	for(partNum = 0; partNum < nA->nParts; partNum++) {
		for(catNum = 0; catNum < nA->tree->model->parts[partNum]->nCat; catNum++) {
			for(stateNum = 0; stateNum < nA->tree->model->parts[partNum]->dim; stateNum++) {
				for(patNum = 0; patNum < nA->tree->data->parts[partNum]->nPatterns; patNum++) {
					if(fabs(nB->cl[partNum][catNum][stateNum][patNum] - nA->cl[partNum][catNum][stateNum][patNum]) > epsilon) {
						printf("  p4_verifyCondLikesFromNodeToNode().  part=%i, category=%i, charState=%i, patNum=%i\n",
							   partNum, catNum, stateNum, patNum);
						diff = fabs(nB->cl[partNum][catNum][stateNum][patNum] - nA->cl[partNum][catNum][stateNum][patNum]);
						printf("  Nodes %i and %i: %g and %g.  diff = %f (%g)\n", 
							   nA->nodeNum, 
							   nB->nodeNum, 
							   nB->cl[partNum][catNum][stateNum][patNum], 
							   nA->cl[partNum][catNum][stateNum][patNum], diff, diff);
						return DIFFERENT;
					}
				}
			}
		}
		if(nB->tree->data->parts[partNum]->nPatterns != nA->tree->data->parts[partNum]->nPatterns) {
			return DIFFERENT;
		}
	}
	
	return SAME;
}


void p4_copyBigPDecksFromNodeToNode(p4_node *nA, p4_node *nB)
{
	int partNum, catNum, i, j;
	
	/*
	  bigPDecks:
	    aNode->nParts
	    aNode->tree->model->parts[partNum]->nCat
	    aNode->tree->model->parts[partNum]->dim, 
	    aNode->tree->model->parts[partNum]->dim 
	*/

	for(partNum = 0; partNum < nA->nParts; partNum++) {
		for(catNum = 0; catNum < nA->tree->model->parts[partNum]->nCat; catNum++) {
			for(i = 0; i < nA->tree->model->parts[partNum]->dim; i++) {
				for(j = 0; j < nA->tree->model->parts[partNum]->dim; j++) {
					nB->bigPDecks[partNum][catNum][i][j] = nA->bigPDecks[partNum][catNum][i][j];
				}
			}
		}
	}	
}

int p4_verifyBigPDecksFromNodeToNode(p4_node *nA, p4_node *nB)
{
	int partNum, catNum, i, j;
	double epsilon, diff;
        
        epsilon = 1.e-15;
	
	/*
	  bigPDecks:
	    aNode->nParts
	    aNode->tree->model->parts[partNum]->nCat
	    aNode->tree->model->parts[partNum]->dim, 
	    aNode->tree->model->parts[partNum]->dim 
	*/

	for(partNum = 0; partNum < nA->nParts; partNum++) {
		for(catNum = 0; catNum < nA->tree->model->parts[partNum]->nCat; catNum++) {
			for(i = 0; i < nA->tree->model->parts[partNum]->dim; i++) {
				for(j = 0; j < nA->tree->model->parts[partNum]->dim; j++) {
					if(fabs(nB->bigPDecks[partNum][catNum][i][j] - 
					   nA->bigPDecks[partNum][catNum][i][j]) > epsilon) {
						printf(" p4_verifyBigPDecksFromNodeToNode()  part=%i, category=%i, i=%i, j=%i\n", partNum, catNum, i, j);
						diff = fabs(nB->bigPDecks[partNum][catNum][i][j] - nA->bigPDecks[partNum][catNum][i][j]);
						printf(" Nodes %i and %i,  (%g, %g) diff= %f (%g)\n", 
							   nA->nodeNum, nB->nodeNum, 
							   nA->bigPDecks[partNum][catNum][i][j], 
							   nB->bigPDecks[partNum][catNum][i][j], diff, diff);
						return DIFFERENT;
					}
				}
			}
		}
	}	
	return SAME;
}



