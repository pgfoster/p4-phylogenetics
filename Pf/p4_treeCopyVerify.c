#include "pftypes.h"
#include "p4_tree.h"
#include "defines.h"
#include "p4_treeCopyVerify.h"
#include "p4_node.h"

void p4_copyCondLikes(p4_tree *treeA, p4_tree *treeB, int doAll)
{
  int   i, j;
  p4_node *nA, *nB;

  for(j = 0; j < treeA->nNodes; j++) {
    i = treeA->preOrder[j];
    if(i != NO_ORDER) {
      nA = treeA->nodes[i];
      if(!(nA->isLeaf)) {
	nB = treeB->nodes[i];
	if(doAll) {
	  p4_copyCondLikesFromNodeToNode(nA, nB);
	}
	else {
	  if(nA->clNeedsUpdating || nB->clNeedsUpdating) {
	    p4_copyCondLikesFromNodeToNode(nA, nB);
	    nA->clNeedsUpdating = 0;
	    nB->clNeedsUpdating = 0;
	  }
	}
      }
    }
  }
}

void p4_copyBigPDecks(p4_tree *treeA, p4_tree *treeB, int doAll)
{
  int   i, j;
  p4_node *nA, *nB;

  if(!doAll) {
    printf("p4_copyBigPDecks() doAll is not set. Programming error?");
    exit(0);
  }
  for(j = 0; j < treeA->nNodes; j++) {
    i = treeA->preOrder[j];
    if(i != NO_ORDER) {
      nA = treeA->nodes[i];
      if(nA != treeA->root) {
	nB = treeB->nodes[i];
	//if(doAll) {
	p4_copyBigPDecksFromNodeToNode(nA, nB);
	//}
	//else {
	//if(nA->brLenChanged || nB->brLenChanged) {
	//	p4_copyBigPDecksFromNodeToNode(nA, nB);
	//	nA->brLenChanged = 0;
	//	nB->brLenChanged = 0;
	//}
	//}
      }
    }
  }
}

void p4_copyModelPrams(p4_tree *treeA, p4_tree *treeB)
{
  int pNum, mtNum, rNum, catNum, dim, i, j;
  p4_modelPart  *a, *b;
  double **aM, **bM, *aV, *bV;
  p4_bigQAndEig *aBQE, *bBQE;

  for(pNum = 0; pNum < treeA->model->nParts; pNum++) {
    a = treeA->model->parts[pNum];
    b = treeB->model->parts[pNum];
    dim = treeA->model->parts[pNum]->dim;

    // comps
    for(mtNum = 0; mtNum < a->nComps; mtNum++) {
      if(a->comps[mtNum]->free) {
	for(i = 0; i < dim; i++) {
	  b->comps[mtNum]->val[i] = a->comps[mtNum]->val[i];
	}
      }
    }
			
	
    // rMatrices
    for(mtNum = 0; mtNum < a->nRMatrices; mtNum++) {
      if(a->rMatrices[mtNum]->free) {
	for(i = 0; i < dim - 2; i++) {
	  for(j = i + 1; j < dim; j++) {
	    b->rMatrices[mtNum]->bigR[i][j] = a->rMatrices[mtNum]->bigR[i][j];
	    b->rMatrices[mtNum]->bigR[j][i] = a->rMatrices[mtNum]->bigR[j][i];
	  }
	}
	if(a->rMatrices[mtNum]->spec == RMATRIX_2P) {
	  b->rMatrices[mtNum]->kappa[0] = a->rMatrices[mtNum]->kappa[0];
	}
      }
    }
			
    // bigQ and eig
    for(mtNum = 0; mtNum < a->nComps; mtNum++) {
      for(rNum = 0; rNum < a->nRMatrices; rNum++) {
	aBQE = a->bigQAndEigThing[mtNum][rNum];
	bBQE = b->bigQAndEigThing[mtNum][rNum];
	if(!aBQE || !bBQE) {
	  printf("----------- aBQE=%li, bBQE=%li\n", (long int)aBQE, (long int)bBQE);
	}
	else{
	  // first do bigQ
	  aM = aBQE->bigQ;
	  bM = bBQE->bigQ;
	  if(bM && aM) {
	    for(i = 0; i < dim; i++) {
	      for(j = 0; j < dim; j++) {
		bM[i][j] = aM[i][j];
	      }
	    }
	  }

	  if(aBQE->qEig && bBQE->qEig) {
	    // do the eigvecs
	    aM = aBQE->qEig->eigvecs;
	    bM = bBQE->qEig->eigvecs;

	    if(bM && aM) {
	      for(i = 0; i < dim; i++) {
		for(j = 0; j < dim; j++) {
		  bM[i][j] = aM[i][j];
		}
	      }
	    }
	    // do the inverseEigvecs
	    aM = aBQE->qEig->inverseEigvecs;
	    bM = bBQE->qEig->inverseEigvecs;
	    if(bM && aM) {
	      for(i = 0; i < dim; i++) {
		for(j = 0; j < dim; j++) {
		  bM[i][j] = aM[i][j];
		}
	      }
	    }
	    // do the eigvals
	    aV = aBQE->qEig->eigvals;
	    bV = bBQE->qEig->eigvals;
	    if(bV && aV) {
	      for(i = 0; i < dim; i++) {
		bV[i] = aV[i];
	      }
	    }
	  }
	}
      }
    }



    // gdasrvs
    for(mtNum = 0; mtNum < a->nGdasrvs; mtNum++) {
      if(a->gdasrvs[mtNum]->free) {
	b->gdasrvs[mtNum]->val[0] = a->gdasrvs[mtNum]->val[0];
	for(catNum = 0; catNum < a->nCat; catNum++) {
	  b->gdasrvs[mtNum]->rates[catNum] = a->gdasrvs[mtNum]->rates[catNum];
	  //b->gdasrvs[mtNum]->freqs[catNum] = a->gdasrvs[mtNum]->freqs[catNum];
	}
      }
    }

    // pInvar
    if(a->pInvar->free) {
      b->pInvar->val[0] = a->pInvar->val[0];
    }

    // mixture
    if(a->isMixture) {
      if(a->mixture->free) {
	for(mtNum = 0; mtNum < a->nCat; mtNum++) {
	  b->mixture->freqs[mtNum] = a->mixture->freqs[mtNum];
	  b->mixture->rates[mtNum] = a->mixture->rates[mtNum];
	}
      }
    }

    // relRate
    if(treeA->model->doRelRates && treeA->model->relRatesAreFree) {
      b->relRate[0] = a->relRate[0];
    }

    // TSCovarion
    if(a->doTSCovarion) {
      printf("p4_copyModelPrams().  Fix me!\n");
    }

  }
  //printf("p4_copyModelPrams().  Done.\n");
}






int p4_verifyIdentityOfTwoTrees(p4_tree *treeA, p4_tree *treeB)
{
  int ret;
  int i, j;
  int isBad = 0;

  //printf("x "); fflush(stdout);
  //printf("%li   %li\n", (long int)treeA, (long int)treeB);
  //p4_dumpTree(treeA);
  //p4_dumpTree(treeB);

  for(j = 0; j < treeA->nNodes; j++) {
    i = treeA->preOrder[j];
    if(i != NO_ORDER) {
      if(treeA->nodes[i]->clNeedsUpdating) {
	printf("Verify: aTree node %i clNeedsUpdating is set.  Bad.\n", i);
	isBad = 1;
      }
      if(treeB->nodes[i]->clNeedsUpdating) {
	printf("Verify: bTree node %i clNeedsUpdating is set.  Bad.\n", i);
	isBad = 1;
      }

      //if(treeA->nodes[i]->brLenChanged) {
      //	printf("Verify: aTree node %i brLenChanged is set.  Bad.\n", i);
      //	isBad = 1;
      //}
      //if(treeB->nodes[i]->brLenChanged) {
      //	printf("Verify: bTree node %i brLenChanged is set.  Bad.\n", i);
      //	isBad = 1;
      //}
    }
  }
  //printf("y "); fflush(stdout);

  ret = p4_verifyModelPrams(treeA, treeB);
  if(ret == DIFFERENT) {
    printf("Verify: model prams are different.  Bad.\n");
    isBad = 1;
  }

  //printf("y2 "); fflush(stdout);
  ret = p4_verifyNodeRelations(treeA, treeB);
  if(ret == DIFFERENT) {
    printf("Verify: nodes relations are different.  Bad.\n");
    isBad = 1;
  }
	
  //printf("y3 "); fflush(stdout);
  ret = p4_verifyNodeInfo(treeA, treeB);
  if(ret == DIFFERENT) {
    printf("Verify: node brLen, or clNeedsUpdating different.  Bad.\n");
    isBad = 1;
  }
	
  //printf("y4 "); fflush(stdout);
  ret = p4_verifyNodeModelUsage(treeA, treeB);
  if(ret == DIFFERENT) {
    printf("Verify: model arrangements are different.  Bad.\n");
    isBad = 1;
  }

  //printf("y5 "); fflush(stdout);
  ret = p4_verifyCondLikes(treeA, treeB);
  if(ret == DIFFERENT) {
    printf("Verify: cond likes are different.  Bad.\n");
    isBad = 1;
  }

  //printf("y6 "); fflush(stdout);
  ret = p4_verifyBigPDecks(treeA, treeB);
  if(ret == DIFFERENT) {
    printf("Verify: bigPDecks are different.  Bad.\n");
    isBad = 1;
  }
  //printf("y7 "); fflush(stdout);

  if(isBad) {
    return DIFFERENT;
  }
  return SAME;
}

int p4_verifyModelPrams(p4_tree *treeA, p4_tree *treeB)
{
  int pNum, mtNum, rNum, catNum, dim, i, j;
  p4_modelPart  *a, *b;
  double diff;
  double epsilon;
  double **aM, **bM, *aV, *bV; // M for matrix, V for vector.
  p4_bigQAndEig *aBQE, *bBQE;
	

  epsilon = 1.e-15;
  //printf("verifyModelPrams()\n");

  for(pNum = 0; pNum < treeA->model->nParts; pNum++) {
    a = treeA->model->parts[pNum];
    b = treeB->model->parts[pNum];
    dim = treeA->model->parts[pNum]->dim;

    // comps
    for(mtNum = 0; mtNum < a->nComps; mtNum++) {
      if(a->comps[mtNum]->free) {
	for(i = 0; i < dim; i++) {
	  if(fabs(a->comps[mtNum]->val[i] - b->comps[mtNum]->val[i]) > epsilon) {
	    printf("verifyModelPrams(): comps differ. part=%i, mtNum=%i, symb=%i, %f and %f\n", 
		   pNum, mtNum, i, a->comps[mtNum]->val[i],b->comps[mtNum]->val[i]);
	    diff = fabs(a->comps[mtNum]->val[i] - b->comps[mtNum]->val[i]);
	    printf("diff= %f (%g)\n", diff, diff);
	    return DIFFERENT;
	  }
	}
      }
    }
			
	
    // rMatrices
    for(mtNum = 0; mtNum < a->nRMatrices; mtNum++) {
      if(a->rMatrices[mtNum]->free) {
	for(i = 0; i < dim - 2; i++) {
	  for(j = i + 1; j < dim; j++) {
	    if(fabs(a->rMatrices[mtNum]->bigR[i][j] - b->rMatrices[mtNum]->bigR[i][j]) > epsilon) {
	      printf("verifyModelPrams(): rMatrices differ.\n");
	      printf("  part=%i, mtNum=%i, bigR[%i][%i] = %f and %f\n", pNum, mtNum, i, j,
		     a->rMatrices[mtNum]->bigR[i][j], b->rMatrices[mtNum]->bigR[i][j]);
	      diff = fabs(a->rMatrices[mtNum]->bigR[i][j] - b->rMatrices[mtNum]->bigR[i][j]);
	      printf("  diff= %f (%g)\n", diff, diff);
	      return DIFFERENT;
	    }
	    else if(fabs(a->rMatrices[mtNum]->bigR[j][i] - b->rMatrices[mtNum]->bigR[j][i]) > epsilon) {
	      printf("verifyModelPrams(): rMatrices differ.\n");
	      printf("  part=%i, mtNum=%i, bigR[%i][%i] = %f and %f\n", pNum, mtNum, j, i,
		     a->rMatrices[mtNum]->bigR[j][i], b->rMatrices[mtNum]->bigR[j][i]);
	      diff = fabs(a->rMatrices[mtNum]->bigR[j][i] - b->rMatrices[mtNum]->bigR[j][i]);
	      printf("  diff= %f (%g)\n", diff, diff);
	      return DIFFERENT;
	    }
	  }
	}
	if(a->rMatrices[mtNum]->spec == RMATRIX_2P) {
	  if(fabs(a->rMatrices[mtNum]->kappa[0] - b->rMatrices[mtNum]->kappa[0]) > epsilon) {
	    printf("verifyModelPrams(): rMatrices kappa differ.\n");
	    printf("  part=%i, mtNum=%i, kappa = %f and %f\n", pNum, mtNum,
		   a->rMatrices[mtNum]->kappa[0], b->rMatrices[mtNum]->kappa[0]);
	    diff = fabs(a->rMatrices[mtNum]->kappa[0] - b->rMatrices[mtNum]->kappa[0]);
	    printf("  diff= %f (%g)\n", diff, diff);
	    return DIFFERENT;
	  }
	}
      }
    }

    // bigQ
    for(mtNum = 0; mtNum < a->nComps; mtNum++) {
      for(rNum = 0; rNum < a->nRMatrices; rNum++) {
	aBQE = a->bigQAndEigThing[mtNum][rNum];
	bBQE = b->bigQAndEigThing[mtNum][rNum];
	if(!aBQE || !bBQE) {
	  printf("-------------------- verifyModelPrams() aBQE=%li, bBQE=%li\n", (long int)aBQE, (long int)bBQE);
	}
	else {
	  // first do the bigQ
	  aM = aBQE->bigQ;
	  bM = bBQE->bigQ;
	  if(bM && aM) {
	    for(i = 0; i < dim; i++) {
	      for(j = 0; j < dim; j++) {
							
		if(fabs(aM[i][j] - bM[i][j]) > epsilon) {
		  printf("verifyModelPrams(): bigQs differ.\n");
		  printf("  part=%i, compNum=%i, rNum=%i, bigQ[%i][%i] = %f and %f\n", 
			 pNum, mtNum, rNum, i, j, aM[i][j], bM[i][j]);
		  diff = fabs(aM[i][j] - bM[i][j]);
		  printf("  diff= %f (%g)\n", diff, diff);
		  return DIFFERENT;
		}
	      }
	    }
	  }
	  if(aBQE->qEig && bBQE->qEig) {
	    // do the eigvecs
	    aM = aBQE->qEig->eigvecs;
	    bM = bBQE->qEig->eigvecs;
	    if(bM && aM) {
	      for(i = 0; i < dim; i++) {
		for(j = 0; j < dim; j++) {
		  if(fabs(aM[i][j] - bM[i][j]) > epsilon) {
		    printf("verifyModelPrams(): eigvecs differ.\n");
		    printf("  part=%i, compNum=%i, rNum=%i, eigvecs[%i][%i] = %f and %f\n", 
			   pNum, mtNum, rNum, i, j, aM[i][j], bM[i][j]);
		    diff = fabs(aM[i][j] - bM[i][j]);
		    printf("  diff= %f (%g)\n", diff, diff);
		    return DIFFERENT;
		  }
		}
	      }
	    }
	    // do the inverseEigvecs
	    aM = aBQE->qEig->inverseEigvecs;
	    bM = bBQE->qEig->inverseEigvecs;
	    if(bM && aM) {
	      for(i = 0; i < dim; i++) {
		for(j = 0; j < dim; j++) {
		  if(fabs(aM[i][j] - bM[i][j]) > epsilon) {
		    printf("verifyModelPrams(): inverseEigvecs differ.\n");
		    printf("  part=%i, compNum=%i, rNum=%i, inverseEigvecs[%i][%i] = %f and %f\n", 
			   pNum, mtNum, rNum, i, j, aM[i][j], bM[i][j]);
		    diff = fabs(aM[i][j] - bM[i][j]);
		    printf("  diff= %f (%g)\n", diff, diff);
		    return DIFFERENT;
		  }
		}
	      }
	    }
	    // do the eigvals
	    aV = aBQE->qEig->eigvals;
	    bV = bBQE->qEig->eigvals;
	    if(bV && aV) {
	      for(i = 0; i < dim; i++) {
		if(fabs(aV[i] - bV[i]) > epsilon) {
		  printf("verifyModelPrams(): eigvals differ.\n");
		  printf("  part=%i, compNum=%i, rNum=%i, eigvals[%i] = %f and %f\n", pNum, mtNum, rNum, i,
			 aV[i], bV[i]);
		  diff = fabs(aV[i] - bV[i]);
		  printf("  diff= %f (%g)\n", diff, diff);
		  return DIFFERENT;
		}
	      }
	    }
	  }
	}
      }
    }
				
			
    // gdasrvs
    for(mtNum = 0; mtNum < a->nGdasrvs; mtNum++) {
      if(a->gdasrvs[mtNum]->free) {
	if(fabs(a->gdasrvs[mtNum]->val[0] - b->gdasrvs[mtNum]->val[0]) > epsilon) {
	  printf("verifyModelPrams(): gdasrvs differ.\n");
	  printf("  part=%i, mtNum=%i, %f and %f\n", pNum, mtNum, 
		 a->gdasrvs[mtNum]->val[0], b->gdasrvs[mtNum]->val[0]);
	  diff = fabs(a->gdasrvs[mtNum]->val[0] - b->gdasrvs[mtNum]->val[0]);
	  printf("  diff= %f (%g)\n", diff, diff);
	  return DIFFERENT;
	}
	for(catNum = 0; catNum < a->nCat; catNum++) {
	  if(fabs(a->gdasrvs[mtNum]->rates[catNum] - b->gdasrvs[mtNum]->rates[catNum]) > epsilon) {
	    printf("verifyModelPrams(): gdasrvs rates differ.\n");
	    printf("  part=%i, mtNum=%i, catNum=%i, %f and %f\n", pNum, mtNum, catNum,
		   a->gdasrvs[mtNum]->rates[catNum], b->gdasrvs[mtNum]->rates[catNum]);
	    diff = fabs(a->gdasrvs[mtNum]->rates[catNum] - b->gdasrvs[mtNum]->rates[catNum]);
	    printf("  diff= %f (%g)\n", diff, diff);
	    return DIFFERENT;
	  }
	  //if(fabs(a->gdasrvs[mtNum]->freqs[catNum] - b->gdasrvs[mtNum]->freqs[catNum]) > epsilon) {
	  //  printf("verifyModelPrams(): gdasrvs freqs differ.\n");
	  //  printf("  part=%i, mtNum=%i, catNum=%i, %f and %f\n", pNum, mtNum, catNum,
	  //   a->gdasrvs[mtNum]->freqs[catNum], b->gdasrvs[mtNum]->freqs[catNum]);
	  //  diff = fabs(a->gdasrvs[mtNum]->freqs[catNum] - b->gdasrvs[mtNum]->freqs[catNum]);
	  //  printf("  diff= %f (%g)\n", diff, diff);
	  //  return DIFFERENT;
	  //}
	}
	  
      }
    }

    // pInvar
    if(a->pInvar->free) {
      if(fabs(a->pInvar->val[0] - b->pInvar->val[0]) > epsilon) {
	printf("verifyModelPrams(): pInvars differ.\n");
	printf("  part=%i, %f and %f\n" , pNum, a->pInvar->val[0], b->pInvar->val[0]);
	diff = fabs(a->pInvar->val[0] - b->pInvar->val[0]);
	printf("  diff= %f (%g)\n", diff, diff);
	return DIFFERENT;
      }
    }


    // mixture
    if(a->isMixture) {
      if(a->mixture->free) {
	for(mtNum = 0; mtNum < a->nCat; mtNum++) {
	  if((fabs(b->mixture->freqs[mtNum] - a->mixture->freqs[mtNum]) > epsilon) || 
	     (fabs(b->mixture->rates[mtNum] - a->mixture->rates[mtNum]) > epsilon)) {
	    printf("verifyModelPrams(): mixtures differ.\n");
	    return DIFFERENT;
	  }
	}
      }
    }

    // relRate
    if(treeA->model->doRelRates && treeA->model->relRatesAreFree) {
      if(fabs(b->relRate[0] - a->relRate[0]) > epsilon) {
	printf("verifyModelPrams(): relRates differ.\n");
	printf("  %f and %f\n", b->relRate[0], a->relRate[0]);
	diff = fabs(b->relRate[0] - a->relRate[0]);
	printf("  diff= %f (%g)\n", diff, diff);
	return DIFFERENT;
      }
    }

    // TSCovarion
    if(a->doTSCovarion) {
      printf("p4_verifyModelPrams().  doTSCovarion. Fix me!\n");
    }
  }
  //printf("p4_verifyModelPrams().  Done\n");
  return SAME;

}
int p4_verifyNodeRelations(p4_tree *treeA, p4_tree *treeB)
{
	int i, j;
	p4_tree *a, *b;

	a = treeA;
	b = treeB;

	for(j = 0; j < a->nNodes; j++) {
		i = a->preOrder[j];
		if(i != NO_ORDER) {
			if(b->nodes[i]->parent) {
				if(!(a->nodes[i]->parent)) {
					return DIFFERENT;
				}
				if(b->nodes[i]->parent != b->nodes[a->nodes[i]->parent->nodeNum]) {
					return DIFFERENT;
				}
			}
			if(a->nodes[i]->parent) {
				if(!(b->nodes[i]->parent)) {
					return DIFFERENT;
				}
			}		

			if(b->nodes[i]->leftChild) {
				if(!(a->nodes[i]->leftChild)) {
					return DIFFERENT;
				}
				if(b->nodes[i]->leftChild != b->nodes[a->nodes[i]->leftChild->nodeNum]) {
					return DIFFERENT;
				}
			}
			if(a->nodes[i]->leftChild) {
				if(!(b->nodes[i]->leftChild)) {
					return DIFFERENT;
				}
			}

			if(b->nodes[i]->sibling) {
				if(!(a->nodes[i]->sibling)) {
					return DIFFERENT;
				}
				if(b->nodes[i]->sibling != b->nodes[a->nodes[i]->sibling->nodeNum]) {
					return DIFFERENT;
				}
			}
			if(a->nodes[i]->sibling) {
				if(!(b->nodes[i]->sibling)) {
					return DIFFERENT;
				}
			}
		}

	}
	if(b->root != b->nodes[a->root->nodeNum]) {
		return DIFFERENT;
	}
	
	return SAME;
}




int p4_verifyNodeInfo(p4_tree *treeA, p4_tree *treeB)
{
	int i, j;
	p4_tree *a, *b;
	double epsilon;
	double diff;

	a = treeA;
	b = treeB;
	epsilon = 1.e-15;

	for(j = 0; j < a->nNodes; j++) {
		i = a->preOrder[j];
		if(i != NO_ORDER) {
			if(a->nodes[i] != a->root) {
				if(fabs(b->nodes[i]->brLen[0] - a->nodes[i]->brLen[0]) > epsilon) {
					printf("p4_verifyNodeInfo() brLens differ.  ");
					printf("  node %i,  %g   %g\n", i, a->nodes[i]->brLen[0], b->nodes[i]->brLen[0]);
					diff = fabs(b->nodes[i]->brLen[0] - a->nodes[i]->brLen[0]);
					printf("  diff = %f (%g)\n", diff, diff);
					return DIFFERENT;
				}
			}
		}
	}	
	for(i = 0; i < a->nNodes; i++) {
		if(b->postOrder[i] != a->postOrder[i]) {
			printf("postOrder different.\n");
			return DIFFERENT;
		}
	}	
	for(i = 0; i < a->nNodes; i++) {
		if(b->preOrder[i] != a->preOrder[i]) {
			printf("preOrder different.\n");
			return DIFFERENT;
		}
	}	
	return SAME;
}




int p4_verifyNodeModelUsage(p4_tree *treeA, p4_tree *treeB)
{
	p4_node *a, *b;
	int  pNum, nNum, i;

	for(i = 0; i < treeA->nNodes; i++) {
		nNum = treeA->preOrder[i];
		if(nNum != NO_ORDER) {
			if(treeA->nodes[nNum] != treeA->root) {
				a = treeA->nodes[nNum];
				b = treeB->nodes[nNum];
				for(pNum = 0; pNum < treeA->nParts; pNum++) {
					if(a->compNums[pNum] != b->compNums[pNum]) {
						printf("compNums differ.  node %i, part %i, compNumA=%i, compNumB=%i\n", 
							   nNum, pNum, a->compNums[pNum], b->compNums[pNum]); 
						return DIFFERENT;
					}
					else if(a->rMatrixNums[pNum] != b->rMatrixNums[pNum]) {
						printf("rMatrixNums differ.  node %i, part %i, compNumA=%i, compNumB=%i\n", 
							   nNum, pNum, a->rMatrixNums[pNum], b->rMatrixNums[pNum]); 
						return DIFFERENT;
					}
					else if(treeA->model->parts[pNum]->nCat > 1) {
						if(a->gdasrvNums[pNum] != b->gdasrvNums[pNum]) {
							printf("gdasrvNums differ.  node %i, part %i, compNumA=%i, compNumB=%i\n", 
								   nNum, pNum, a->gdasrvNums[pNum], b->gdasrvNums[pNum]); 
							return DIFFERENT;
						}
					}
				}
			}
		}
	}
	
	return SAME;
}


int p4_verifyCondLikes(p4_tree *treeA, p4_tree *treeB)
{	
	int   i, j, ret;
	p4_node *nA, *nB;

	for(j = 0; j < treeA->nNodes; j++) {
		i = treeA->preOrder[j];
		if(i != NO_ORDER) {
			nA = treeA->nodes[i];
			if(!(nA->isLeaf)) {
				nB = treeB->nodes[i];
				//printf("n%i ", i);
				ret = p4_verifyCondLikesFromNodeToNode(nA, nB);
				if(ret == DIFFERENT) {
					return DIFFERENT;
				}
			}
		}
	}
	return SAME;
}


int p4_verifyBigPDecks(p4_tree *treeA, p4_tree *treeB)
{	
	int   i,j, ret;
	p4_node *nA, *nB;

	for(j = 0; j < treeA->nNodes; j++) {
		i = treeA->preOrder[j];
		if(i != NO_ORDER) {
			nA = treeA->nodes[i];
			if(nA != treeA->root) {
				nB = treeB->nodes[i];
				//printf("n%i ", i);
				ret = p4_verifyBigPDecksFromNodeToNode(nA, nB);
				if(ret == DIFFERENT) {
					return DIFFERENT;
				}
			}
		}
	}
	return SAME;
}




