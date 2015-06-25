#include "pftypes.h"
#include "p4_tree.h"
#include "p4_node.h"
#include "gamma.h"    // Yang funcs, DiscreteGamma et al.
#include "defines.h"
#include "util.h"
#include "eig.h"
/* #include <stdio.h> */
/* #include <stdlib.h> */
#include <math.h>    // log
#include "pmatrices.h"

//static  int anInt = 0;  // for a recursion, below.

p4_tree *p4_newTree(int nNodes, int nLeaves, int *preOrder, int *postOrder, double *partLikes, data *aData, p4_model *aModel)
{
	p4_tree	*aTree;
	int		i;

	//printf("p4_newTree here.\n");

	aTree = malloc(sizeof(p4_tree));
	if(!aTree) {
		printf("Failed to alloc memory for tree.\n");
		exit(1);
	}
	
	aTree->nNodes = nNodes;
	aTree->nLeaves = nLeaves;
	aTree->root = NULL;
	aTree->data = aData;
	aTree->model = aModel;
	if(aData) {
		aTree->nParts = aData->nParts;
	} else {  // should never happen.
		aTree->nParts = 0;
	}
		

	aTree->nodes = malloc(nNodes * sizeof(p4_node *));
	if(!aTree->nodes) {
		printf("Failed to alloc memory for tree->nodes.\n");
		exit(1);
	}

	aTree->preOrder = preOrder;  // numeric array
	aTree->postOrder = postOrder;
	aTree->partLikes = partLikes;
	//for(i = 0; i < nNodes; i++) {
	//	printf("node %i, ", i);
	//	printf("preOrder = %i\n", preOrder[i]);
	//}

	aTree->ints = malloc(nNodes * sizeof(int));
	if(!aTree->ints) {
		printf("Failed to alloc memory for tree->ints.\n");
		exit(1);
	}

	aTree->stack = malloc(nNodes * sizeof(p4_node *));
	if(!aTree->stack) {
		printf("Failed to alloc memory for tree->stack.\n");
		exit(1);
	}

	for(i = 0; i < nNodes; i++) {
		aTree->nodes[i] = NULL;
		aTree->ints[i] = -1;
		aTree->stack[i] = NULL;
	}
	
	aTree->logLike = 0.0;

	aTree->simSequences = NULL;
	aTree->internalSequences = NULL;
	


#if 0 
	aTree->expPi = NULL;
	aTree->path = malloc(nNodes * sizeof(int));
	if(!aTree->path) {
		printf("Failed to alloc memory for tree path\n");
		exit(1);
	}
	for(i = 0; i < nNodes; i++) {
		aTree->path[i] = -1;
	}
#endif
	//printf("newTree: returning %li\n", (long int) aTree);
	return aTree;

}

void p4_freeTree(p4_tree *aTree)
{
			
	//printf("starting p4_freeTree\n");
	int i, j;
	int nInternals = 0;

	// I need to free the sequences used for simulation before I free
	// the nodes.  The simSequences are composed of sequences from the
	// parts, and internalSequences from the tree.  The former is
	// freed by the part, and the latter is freed here.  Then the
	// remaining skeleton of simSequences is freed.

	//printf("p4_tree.c  about to freeTree.\n");
	if(aTree->internalSequences) {
		nInternals = aTree->nNodes - aTree->nLeaves;
		for(i = 0; i < aTree->nParts; i++) {
			//k = 0;
			//for(j = 0; j < aTree->nNodes; j++) {
			//	if(aTree->nodes[j]->isLeaf) {
			//		// pass
			//	} else {
			//		free(aTree->internalSequences[i][k]);
			//		aTree->internalSequences[i][k] = NULL;
			//		k++;
			//	}
			for(j = 0; j < nInternals; j++) {
				free(aTree->internalSequences[i][j]);
				aTree->internalSequences[i][j] = NULL;
			}
			free(aTree->internalSequences[i]);
			aTree->internalSequences[i] = NULL;
		}
		free(aTree->internalSequences);
		aTree->internalSequences = NULL;
	}

	// Just free the skeleton
	if(aTree->simSequences) {
		for(i = 0; i < aTree->root->nParts; i++) {
			for(j = 0; j < aTree->nNodes; j++) {
				aTree->simSequences[i][j] = NULL;
			}
			free(aTree->simSequences[i]);
			aTree->simSequences[i] = NULL;
		}
		free(aTree->simSequences);
		aTree->simSequences = NULL;
	}

#if 0
	
	if(aTree->expPi) {
		for(i = 0; i < aTree->root->nParts; i++) {
			free_pdmatrix(aTree->expPi[i]);
			aTree->expPi[i] = NULL;
		}
	free(aTree->expPi);
	aTree->expPi = NULL;
	}

#endif

	// We don't want to free the nodes, but do want to free the nodes array.
	free(aTree->nodes);
	aTree->nodes = NULL;
	aTree->root = NULL;
	//free(aTree->preOrder);
	aTree->preOrder = NULL;
	//free(aTree->postOrder);
	aTree->postOrder = NULL;
	free(aTree->ints);
	aTree->ints = NULL;
	free(aTree->stack);
	aTree->stack = NULL;
	free(aTree);
	aTree = NULL;
	//printf("p4_tree.c: finished free-ing tree\n");
}

void p4_dumpTree(p4_tree *aTree)
{
	//int i;

	printf("\ndumpTree:\n");
	printf("    %15s %i\n", "nNodes", aTree->nNodes);
	printf("    %15s %li\n", "data", (long int)aTree->data);
	printf("    %15s %li\n", "model", (long int)aTree->model);
	printf("    %15s %i\n", "root", aTree->root->nodeNum);
#if 0
	printf("    %15s ", "preOrder");
	for(i = 0; i < aTree->nNodes; i++) {
		printf("%i ", aTree->preOrder[i]);
	}
	printf("\n");
	printf("    %15s ", "postOrder");
	for(i = 0; i < aTree->nNodes; i++) {
		printf("%i ", aTree->postOrder[i]);
	}
	printf("\n");
#endif
}


void p4_setPrams(p4_tree *aTree)
{
	int pNum;

	for(pNum = 0; pNum < aTree->model->nParts; pNum++) {
		p4_setPramsPart(aTree, pNum);
	}
}



void p4_setPramsPart(p4_tree *aTree, int pNum)
{

  int mNum, nNum, cNum, rNum;
  int dim;
  int i, j;
  int ret;
  double alpha, beta;
  p4_modelPart   *mp;
  part           *dp;
  p4_gdasrv      *g;
  p4_rMatrix     *r;
  p4_node        *n;
  p4_comp        *c;
  p4_bigQAndEig  *aQE;

  //printf("p4_setPramsPart here.\n");
  mp = aTree->model->parts[pNum];
  if(mp->nCat > 1) {
    // Set gdasrv->rates using Yang's DiscreteGamma function.
    for(mNum = 0; mNum < mp->nGdasrvs; mNum++) {
      g = mp->gdasrvs[mNum];
      // All rates were set when the gdasrv was given its
      // val, in pf_p4_setGdasrvVal() in pfmodule.c.  We
      // need only re-do the ones with a free shape.
      if(g->free) {
	//printf("part %i, gdasrv %i, is free, with val %f\n", pNum, mNum, g->val[0]);
	// int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median)
	// Here alfa = beta.  The freqK array is used as a
	// working space for internal calcs by
	// DiscreteGamma, but boringly ends up as all 1.0/K
	DiscreteGamma(g->freqs, g->rates, g->val[0], g->val[0], mp->nCat, 0);
	if(0) {
	  printf("gdasrv[%i] = %f, free=%i\n", mNum, g->val[0], g->free);
	  for(i = 0; i < mp->nCat; i++) {
	    printf("  %i  %f   %f\n", i, g->freqs[i], g->rates[i]);
	  }
	}
      }
    }
  }
	


#if 0  // I'll need to do this when I have a free kappa.  And doTSCovarion.
  // rMatrix
  //printf("\np4_tree setPramsPart(): about to set rMatrix\n");
  mp = aTree->model->parts[pNum];
  dim = mp->dim;
  for(mNum = 0; mNum < mp->nRMatrices; mNum++) {
    r = mp->rMatrices[mNum];
    if(r->free) {
      //printf("part %i, rMatrix %i, is free.  spec = %i\n", pNum, mNum, r->spec);
      if(r->spec == RMATRIX_SPECIFIED ||
	 r->spec == RMATRIX_POISSON) {
	k = 0;
	for(i = 0; i < dim - 2; i++) {
	  for(j = i + 1; j < dim; j++) {
	    //printf("i=%i, j=%i, k=%i\n", i,j,k);
	    r->bigR[i][j] = r->val[k];
	    r->bigR[j][i] = r->val[k];
	    k++;
	  }
	}
	r->bigR[dim - 2][dim - 1] = 1.0;
	r->bigR[dim - 1][dim - 2] = 1.0;
	//dump_psdmatrix(r->bigR, dim);
      }
      else {
	printf("I don't recognize the rMatrix spec\n");
	exit(1);
      }
    }
  }
	
#endif


  // Set rMatrix for 2p (ie nst=2) models (k2p and hky)
  mp = aTree->model->parts[pNum];
  dim = mp->dim;
  for(mNum = 0; mNum < mp->nRMatrices; mNum++) {
    r = mp->rMatrices[mNum];
    if(r->free && (r->spec == RMATRIX_2P)) {
      //printf("part %i, rMatrix %i, is free.  spec = %i\n", pNum, mNum, r->spec);
      alpha = 1.0 / 3.0;
      beta = alpha * r->kappa[0];
      r->bigR[0][0] = 0.0;
      r->bigR[0][1] = alpha;
      r->bigR[0][2] = beta;
      r->bigR[0][3] = alpha;
      r->bigR[1][0] = alpha;
      r->bigR[1][1] = 0.0;
      r->bigR[1][2] = alpha;
      r->bigR[1][3] = beta;
      r->bigR[2][0] = beta;
      r->bigR[2][1] = alpha;
      r->bigR[2][2] = 0.0;
      r->bigR[2][3] = alpha;
      r->bigR[3][0] = alpha;
      r->bigR[3][1] = beta;
      r->bigR[3][2] = alpha;
      r->bigR[3][3] = 0.0;
    }
  }
	


#if 1
  // doTSCovarion
  mp = aTree->model->parts[pNum];
  dp = aTree->data->parts[pNum];
  if(mp->doTSCovarion) {
    //printf("p4_setPrams().  s1=%f, s2=%f, pOn=%f, pOff=%f\n", 
    //	   mp->tSCov->s1[0], mp->tSCov->s2[0], mp->tSCov->pOn, mp->tSCov->pOff);
    r = mp->rMatrices[0]; // only one is possible, this week
    for(i = 0; i < dp->dim; i++) {
      for(j = dp->dim; j < mp->dim; j++) {
	if((i + dp->dim) == j) {
	  r->bigR[i][j] = mp->tSCov->s1[0];
	} else {
	  r->bigR[i][j] = 0.0;
	}
      }
    }
    for(i = dp->dim; i < mp->dim; i++) {
      for(j = 0; j < dp->dim; j++) {
	if(i == (j + dp->dim)) {
	  r->bigR[i][j] = mp->tSCov->s2[0];
	} else {
	  r->bigR[i][j] = 0.0;
	}
      }
    }
    for(i = dp->dim; i < mp->dim; i++) {
      for(j = dp->dim; j < mp->dim; j++) {
	r->bigR[i][j] = 0.0;
      }
    }

    // The long comp, in mp->comps[0]->val, is calculated from the short comp, in mp->tSCov->halfComp
    for(i = 0; i< dp->dim; i++) {
      mp->comps[0]->val[i] = mp->tSCov->halfComp[i] * mp->tSCov->pOn;
      mp->comps[0]->val[i + dp->dim] = mp->tSCov->halfComp[i] * mp->tSCov->pOff;
    }
  }
	
#endif			


	
#if 0
  // Write out some model usage info. ie what nodes use what comps etc.
  printf("part %i\n", pNum);
  for(nNum = 0; nNum < aTree->nNodes; nNum++) {
    printf("    node %2i ", nNum);
    n = aTree->nodes[nNum];
    printf(" c%i", n->compNums[pNum]);
    if(n != aTree->root) {
      printf(" r%i", n->rMatrixNums[pNum]);
      printf(" g%i", n->gdasrvNums[pNum]);
    }
    printf("\n");
  }
	
#endif	

#if 0
  printf("bigR\n");
  dump_psdmatrix(aTree->model->parts[pNum]->rMatrices[0]->bigR, aTree->model->parts[pNum]->dim);
  printf("comp: ");
  dump_pdvector(aTree->model->parts[pNum]->comps[0]->val, aTree->model->parts[pNum]->dim);
#endif

#if 0
  {
    double sum=0.0;
    for(i = 0; i < aTree->model->parts[pNum]->dim; i++) {
      sum += aTree->model->parts[pNum]->comps[0]->val[i];
      printf(" %7.4f", sum);
    }
    printf("\n");
  }
			
#endif


  // bigQ
	
  /*
    Each part has a bigQAndEigThing.  It is nComps * nRMatrices.
    Each slot has a bigQAndEig struct, but at the start of play they
    are all empty.  They are only filled and the qEig calculated if
    needed.
  */

  //printf("about to set bigQ.\n");

  // Nasty surprise.  If a charFreq is zero, eg if the charFreq is
  // (0.4, 0.4, 0.2, 0.0), then the bigP's are *wrong* --- way
  // wrong.  However, its ok if eg the charFreq is (0.4, 0.4,
  // 0.199999, 0.000001).  In previous versions I used to check and
  // adjust, but here it is just prevented, elsewhere.  In
  // optimizations, the PIVEC_MIN is enforced.


#if 0
  // Check comps anyway.  This interferes with mcmc when I do verifyModelPrams().
  {
    double sum=0.0;

    mp = aTree->model->parts[pNum];
    //printf("    part %i, nComps=%i, nRMatrices=%i, nCat=%i\n", pNum, mp->nComps, mp->nRMatrices, mp->nCat);
    for(i = 0; i < mp->nComps; i++) {
      for(j = 0; j < mp->dim; j++) {
	if(mp->comps[i]->val[j] < (0.5 * PIVEC_MIN)) {
	  printf("p4_setPramsPart()  part %i, comp %i, value %i is %g   Bad.\n", pNum, i, j, mp->comps[i]->val[j]);
	  exit(1);
	}
      }
      // check that the sum is 1.0
      sum = 0.0;
      for(j = 0; j < mp->dim; j++) {
	sum += mp->comps[i]->val[j];
      }
      // normalize if needed
      if(fabs(sum - 1.0) > 1e-15) {
	printf("   ** p4_setPramsPart().  normalizing comp. diff=%f (%g)\n", fabs(sum - 1.0), fabs(sum - 1.0));
	if(fabs(sum - 1.0) > 0.001) {
	  printf("p4_setPramsPart()  part %i, comp %i, sum of vals = %f.  Bad.\n", pNum, i, sum);
	}
	for(j = 0; j < mp->dim; j++) {
	  mp->comps[i]->val[j] /= sum;
	}
      }
      // check again
      sum = 0.0;
      for(j = 0; j < mp->dim; j++) {
	sum += mp->comps[i]->val[j];
      }
      // if its wrong this time, give up  (diff of 1e-16 was too small, raised to 1e-14)
      if(fabs(sum - 1.0) > 1e-14) {
	printf("p4_setPramsPart()  part %i, comp %i, values do not sum to 1.0.  sum=%g\n", pNum, i, sum);
	printf("  sum - 1.0 = %g\n", sum - 1.0);
	exit(1);
      }
    }
  }
#endif

#if 1
  // Check comps anyway.  Same as above, but no normalizing.
  {
    double sum=0.0;

    mp = aTree->model->parts[pNum];
    //printf("    part %i, nComps=%i, nRMatrices=%i, nCat=%i\n", pNum, mp->nComps, mp->nRMatrices, mp->nCat);
    for(i = 0; i < mp->nComps; i++) {
      for(j = 0; j < mp->dim; j++) {
	if(mp->comps[i]->val[j] < (0.5 * PIVEC_MIN)) {
	  printf("p4_setPramsPart()  part %i, comp %i, value %i is %g   Bad.\n", pNum, i, j, mp->comps[i]->val[j]);
	  exit(1);
	}
      }
      // check
      sum = 0.0;
      for(j = 0; j < mp->dim; j++) {
	sum += mp->comps[i]->val[j];
      }
      // if its wrong this time, give up  (diff of 1e-16 was too small, raised to 1e-14)
      if(fabs(sum - 1.0) > 1e-14) {
	printf("**p4_setPramsPart()  part %i, comp %i, values do not sum to 1.0.  sum=%g\n", pNum, i, sum);
	printf("**  sum - 1.0 = %g\n", sum - 1.0);
	exit(1);
      }
    }
  }
#endif

  // turn on all needsReset
  mp = aTree->model->parts[pNum];
  //printf("    part %i, nComps=%i, nRMatrices=%i, nCat=%i\n", pNum, mp->nComps, mp->nRMatrices, mp->nCat);
  for(i = 0; i < mp->nComps; i++) {
    for(j = 0; j < mp->nRMatrices; j++) {
      mp->bQETneedsReset[(i * mp->nRMatrices) + j] = 1;
    }
  }
	

  mp = aTree->model->parts[pNum];
  dim = mp->dim;
  for(nNum = 0; nNum < aTree->nNodes; nNum++) {
    n = aTree->nodes[nNum];
    if(n != aTree->root) {
      cNum = n->compNums[pNum];
      rNum = n->rMatrixNums[pNum];
      //printf("cNum = %i, rNum = %i\n", cNum, rNum);
      c = mp->comps[cNum];
      r = mp->rMatrices[rNum];
      aQE = aTree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];
      if(mp->bQETneedsReset[(cNum * mp->nRMatrices) + rNum]) {
	//printf("Doing comp %i, rMatrix %i\n", cNum, rNum);
	if(!aQE->bigQ) {
	  aQE->bigQ = psdmatrix(dim);
	  setBigQFromRMatrixDotCharFreq(aQE->bigQ, r->bigR, c->val, dim);
	  normalizeBigQ(aQE->bigQ, c->val, dim);
	  aQE->qEig = allocEig(dim, aQE->bigQ);
	}
	else {
	  //dump_psdmatrix(r->bigR, dim);
	  //setBigQFromRMatrixDotCharFreq(m->bigQ, theRMatrix, theCharFreq, m->dim);
	  setBigQFromRMatrixDotCharFreq(aQE->bigQ, r->bigR, c->val, dim);
	  //dump_psdmatrix(aQE->bigQ, dim);
	  //normalizeBigQ(m->bigQ, theCharFreq, m->dim);
	  normalizeBigQ(aQE->bigQ, c->val, dim);
	  //dump_psdmatrix(aQE->bigQ, dim);
	}
	ret = eigensystem(aQE->qEig);
	if(ret) {
	  printf("There is a problem with the eigensystem.\n");
	  exit(1);
	}
	mp->bQETneedsReset[(cNum * mp->nRMatrices) + rNum] = 0;  // No need to do it again for this 
	// combination of cNum and rNum on another node.
      }
    }
  }
  

#if 0
  for(cNum = 0; cNum < aTree->model->parts[pNum]->nComps; cNum++) {
    for(rNum = 0; rNum < aTree->model->parts[pNum]->nRMatrices; rNum++) {
      aQE = aTree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];
      printf("cNum=%i, rNum=%i, mp->needsReset=%i, ", cNum, rNum, 
	     mp->bQETneedsReset[(cNum * mp->nRMatrices) + rNum]);
      printf("aQE->bigQ=%li\n", (long int)(aQE->bigQ));
    }
  }
#endif 
	

#if 0
  printf("bigQ\n");
  dump_psdmatrix(aTree->model->parts[pNum]->bigQAndEigThing[0][0]->bigQ, aTree->model->parts[pNum]->dim);
#endif
				
  //printf("about to set bigP\n");
  // for each node except the root, "calculateBigPDecks()"
  for(nNum = 0; nNum < aTree->nNodes; nNum++) {
    n = aTree->nodes[nNum];
    if(n != aTree->root) {
      p4_calculateBigPDecksPart(n, pNum);
    }
  }

#if 0
  printf("bigP for node 1\n");
  dump_psdmatrix(aTree->nodes[1]->bigPDecks[0][0], aTree->model->parts[pNum]->dim);	
  //dump_psdmatrix(aTree->nodes[1]->bigPDecks[0][1], aTree->model->parts[pNum]->dim);
  //printf("bigP for node 2\n");
  //dump_psdmatrix(aTree->nodes[2]->bigPDecks[0][0], aTree->model->parts[pNum]->dim);	
#endif


	

}


void p4_setPramsTest(p4_tree *aTree)
{
	int pNum;

	for(pNum = 0; pNum < aTree->model->nParts; pNum++) {
		p4_setPramsPartTest(aTree, pNum);
	}
}


void p4_setPramsPartTest(p4_tree *aTree, int pNum)
{

  int mNum;
  int nNum, cNum, rNum;
  int dim;
  int i, j;
  int ret;
  double alpha, beta;
  p4_modelPart   *mp;
  //part           *dp;
  p4_gdasrv      *g;
  p4_rMatrix     *r;
  p4_node        *n;
  p4_comp        *c;
  p4_bigQAndEig  *aQE;

  printf("a p4_setPramsPartTest here.  pNum=%i\n", pNum);

#if 1
  mp = aTree->model->parts[pNum];
  if(mp->nCat > 1) {
    // Set gdasrv->rates using Yang's DiscreteGamma function.
    for(mNum = 0; mNum < mp->nGdasrvs; mNum++) {
      g = mp->gdasrvs[mNum];
      // All rates were set when the gdasrv was given its
      // val, in pf_p4_setGdasrvVal() in pfmodule.c.  We
      // need only re-do the ones with a free shape.
      if(g->free) {
	printf("part %i, gdasrv %i, is free, with val %f\n", pNum, mNum, g->val[0]);
	// int DiscreteGamma (double freqK[], double rK[], double alfa, double beta, int K, int median)
	// Here alfa = beta.  The freqK array is used as a
	// working space for internal calcs by
	// DiscreteGamma, but boringly ends up as all 1.0/K
	DiscreteGamma(g->freqs, g->rates, g->val[0], g->val[0], mp->nCat, 0);
	if(1) {
	  printf("gdasrv[%i] = %f, free=%i\n", mNum, g->val[0], g->free);
	  for(i = 0; i < mp->nCat; i++) {
	    printf("  %i  %f   %f\n", i, g->freqs[i], g->rates[i]);
	  }
	}
      }
    }
  }
#endif

#if 1
  // Set rMatrix for 2p (ie nst=2) models (k2p and hky)
  mp = aTree->model->parts[pNum];
  dim = mp->dim;
  for(mNum = 0; mNum < mp->nRMatrices; mNum++) {
    r = mp->rMatrices[mNum];
    if(r->free && (r->spec == RMATRIX_2P)) {
      //printf("part %i, rMatrix %i, is free.  spec = %i\n", pNum, mNum, r->spec);
      alpha = 1.0 / 3.0;
      beta = alpha * r->kappa[0];
      r->bigR[0][0] = 0.0;
      r->bigR[0][1] = alpha;
      r->bigR[0][2] = beta;
      r->bigR[0][3] = alpha;
      r->bigR[1][0] = alpha;
      r->bigR[1][1] = 0.0;
      r->bigR[1][2] = alpha;
      r->bigR[1][3] = beta;
      r->bigR[2][0] = beta;
      r->bigR[2][1] = alpha;
      r->bigR[2][2] = 0.0;
      r->bigR[2][3] = alpha;
      r->bigR[3][0] = alpha;
      r->bigR[3][1] = beta;
      r->bigR[3][2] = alpha;
      r->bigR[3][3] = 0.0;
    }
  }
#endif
	
  // bigQ
	
  /*
    Each part has a bigQAndEigThing.  It is nComps * nRMatrices.
    Each slot has a bigQAndEig struct, but at the start of play they
    are all empty.  They are only filled and the qEig calculated if
    needed.
  */

  //printf("about to set bigQ.\n");

  // Nasty surprise.  If a charFreq is zero, eg if the charFreq is
  // (0.4, 0.4, 0.2, 0.0), then the bigP's are *wrong* --- way
  // wrong.  However, its ok if eg the charFreq is (0.4, 0.4,
  // 0.199999, 0.000001).  In previous versions I used to check and
  // adjust, but here it is just prevented, elsewhere.  In
  // optimizations, the PIVEC_MIN is enforced.


#if 0
  // Check comps anyway.  This interferes with mcmc when I do verifyModelPrams().
  {
    double sum=0.0;

    mp = aTree->model->parts[pNum];
    //printf("    part %i, nComps=%i, nRMatrices=%i, nCat=%i\n", pNum, mp->nComps, mp->nRMatrices, mp->nCat);
    for(i = 0; i < mp->nComps; i++) {
      for(j = 0; j < mp->dim; j++) {
	if(mp->comps[i]->val[j] < (0.5 * PIVEC_MIN)) {
	  printf("p4_setPramsPart()  part %i, comp %i, value %i is %g   Bad.\n", pNum, i, j, mp->comps[i]->val[j]);
	  exit(1);
	}
      }
      // check that the sum is 1.0
      sum = 0.0;
      for(j = 0; j < mp->dim; j++) {
	sum += mp->comps[i]->val[j];
      }
      // normalize if needed
      if(fabs(sum - 1.0) > 1e-15) {
	printf("   ** p4_setPramsPart().  normalizing comp. diff=%f (%g)\n", fabs(sum - 1.0), fabs(sum - 1.0));
	if(fabs(sum - 1.0) > 0.001) {
	  printf("p4_setPramsPart()  part %i, comp %i, sum of vals = %f.  Bad.\n", pNum, i, sum);
	}
	for(j = 0; j < mp->dim; j++) {
	  mp->comps[i]->val[j] /= sum;
	}
      }
      // check again
      sum = 0.0;
      for(j = 0; j < mp->dim; j++) {
	sum += mp->comps[i]->val[j];
      }
      // if its wrong this time, give up  (diff of 1e-16 was too small, raised to 1e-14)
      if(fabs(sum - 1.0) > 1e-14) {
	printf("p4_setPramsPart()  part %i, comp %i, values do not sum to 1.0.  sum=%g\n", pNum, i, sum);
	printf("  sum - 1.0 = %g\n", sum - 1.0);
	exit(1);
      }
    }
  }
#endif

#if 1
  // Check comps anyway.  Same as above, but no normalizing.  Its only
  // a check -- it doesn't actually do anything.
  {
    double sum=0.0;

    mp = aTree->model->parts[pNum];
    //printf("    part %i, nComps=%i, nRMatrices=%i, nCat=%i\n", pNum, mp->nComps, mp->nRMatrices, mp->nCat);
    for(i = 0; i < mp->nComps; i++) {
      for(j = 0; j < mp->dim; j++) {
	if(mp->comps[i]->val[j] < (0.5 * PIVEC_MIN)) {
	  printf("p4_setPramsPart()  part %i, comp %i, value %i is %g   Bad.\n", pNum, i, j, mp->comps[i]->val[j]);
	  exit(1);
	}
      }
      // check
      sum = 0.0;
      for(j = 0; j < mp->dim; j++) {
	sum += mp->comps[i]->val[j];
      }
      // if its wrong this time, give up  (diff of 1e-16 was too small, raised to 1e-14)
      if(fabs(sum - 1.0) > 1e-14) {
	printf("**p4_setPramsPart()  part %i, comp %i, values do not sum to 1.0.  sum=%g\n", pNum, i, sum);
	printf("**  sum - 1.0 = %g\n", sum - 1.0);
	exit(1);
      }
    }
  }
#endif

#if 1
  // turn on all needsReset
  mp = aTree->model->parts[pNum];
  //printf("    part %i, nComps=%i, nRMatrices=%i, nCat=%i\n", pNum, mp->nComps, mp->nRMatrices, mp->nCat);
  for(i = 0; i < mp->nComps; i++) {
    for(j = 0; j < mp->nRMatrices; j++) {
      mp->bQETneedsReset[(i * mp->nRMatrices) + j] = 1;
    }
  }
	

  mp = aTree->model->parts[pNum];
  dim = mp->dim;
  for(nNum = 0; nNum < aTree->nNodes; nNum++) {
    n = aTree->nodes[nNum];
    if(n != aTree->root) {
      cNum = n->compNums[pNum];
      rNum = n->rMatrixNums[pNum];
      //printf("cNum = %i, rNum = %i\n", cNum, rNum);
      c = mp->comps[cNum];
      r = mp->rMatrices[rNum];
      aQE = aTree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];
      if(mp->bQETneedsReset[(cNum * mp->nRMatrices) + rNum]) {
	//printf("Doing comp %i, rMatrix %i\n", cNum, rNum);
	if(!aQE->bigQ) {
	  aQE->bigQ = psdmatrix(dim);
	  setBigQFromRMatrixDotCharFreq(aQE->bigQ, r->bigR, c->val, dim);
	  normalizeBigQ(aQE->bigQ, c->val, dim);
	  aQE->qEig = allocEig(dim, aQE->bigQ);
	}
	else {
	  //dump_psdmatrix(r->bigR, dim);
	  //setBigQFromRMatrixDotCharFreq(m->bigQ, theRMatrix, theCharFreq, m->dim);
	  setBigQFromRMatrixDotCharFreq(aQE->bigQ, r->bigR, c->val, dim);
	  //dump_psdmatrix(aQE->bigQ, dim);
	  //normalizeBigQ(m->bigQ, theCharFreq, m->dim);
	  normalizeBigQ(aQE->bigQ, c->val, dim);
	  //dump_psdmatrix(aQE->bigQ, dim);
	}
	ret = eigensystem(aQE->qEig);
	if(ret) {
	  printf("There is a problem with the eigensystem.\n");
	  exit(1);
	}
	mp->bQETneedsReset[(cNum * mp->nRMatrices) + rNum] = 0;  // No need to do it again for this 
	// combination of cNum and rNum on another node.
      }
    }
  }
#endif 

#if 0
  for(cNum = 0; cNum < aTree->model->parts[pNum]->nComps; cNum++) {
    for(rNum = 0; rNum < aTree->model->parts[pNum]->nRMatrices; rNum++) {
      aQE = aTree->model->parts[pNum]->bigQAndEigThing[cNum][rNum];
      printf("cNum=%i, rNum=%i, mp->needsReset=%i, ", cNum, rNum, 
	     mp->bQETneedsReset[(cNum * mp->nRMatrices) + rNum]);
      printf("aQE->bigQ=%li\n", (long int)(aQE->bigQ));
    }
  }
#endif 
	

#if 0
  printf("bigQ\n");
  dump_psdmatrix(aTree->model->parts[pNum]->bigQAndEigThing[0][0]->bigQ, aTree->model->parts[pNum]->dim);
#endif

  
#if 1
  //printf("about to set bigP\n");
  // for each node except the root, "calculateBigPDecks()"
  for(nNum = 0; nNum < aTree->nNodes; nNum++) {
    if(nNum != NO_ORDER) {
      n = aTree->nodes[nNum];
      if(n != aTree->root) {
	p4_calculateBigPDecksPart(n, pNum);
      }
    }
  }
#endif

#if 0
  printf("bigP for node 1\n");
  dump_psdmatrix(aTree->nodes[1]->bigPDecks[0][0], aTree->model->parts[pNum]->dim);	
  //dump_psdmatrix(aTree->nodes[1]->bigPDecks[0][1], aTree->model->parts[pNum]->dim);
  //printf("bigP for node 2\n");
  //dump_psdmatrix(aTree->nodes[2]->bigPDecks[0][0], aTree->model->parts[pNum]->dim);	
#endif

}

void p4_calculateAllBigPDecksAllParts(p4_tree *aTree)
{
  int pNum, nNum;
  p4_node *n;

  for(pNum = 0; pNum < aTree->model->nParts; pNum++) {
    for(nNum = 0; nNum < aTree->nNodes; nNum++) {
      if(nNum != NO_ORDER) {
	n = aTree->nodes[nNum];
	if(n != aTree->root) {
	  p4_calculateBigPDecksPart(n, pNum);
	}
      }
    }
  }
}



/* 

Flags:
======

  - p4_tree->model->parts[pNum]->clNeedsUpdating     Not used anymore!
                              The prams have changed, so all cond
                              likes need updating in that part.

  - p4_node->clNeedsUpdating
                              If a brLen changed, then the cond likes
                              need updating from there down to the
                              root, but not above.  This flag says
                              which cl's to update.
                              It does not seem to be used except by newt.

 setPrams()
  
 setPramsPart() 
  - recalculates all Q-matrices and resets eigensystems
  - recalculates P-matrices for all the nodes,

 treeLogLike()
  - sets condLikes of interior nodes
  - calculates the logLike


*/






double p4_treeLogLike(p4_tree *aTree, int getSiteLikes)
{
  double lnL = 0.0;
  int    i, j, pNum;
  p4_node *aNode;
  part   *dp;
	
  for(j = 0; j < aTree->nNodes; j++) {
    i = aTree->postOrder[j];
    if(i != NO_ORDER) {
      aNode = aTree->nodes[i];
      if(!aNode->isLeaf) {
	p4_setConditionalLikelihoodsOfInteriorNode(aNode);
      }
    }
  }
		

#if 0
  //printf("cl for node 1\n");
  //p4_dumpCL(aTree->nodes[1]->cl);
  printf("cl for node 0\n");
  p4_dumpCL(aTree->nodes[0]->cl);
  lnL = 0.0;
  for(i = 0; i < 4; i++) {
    lnL += (aTree->nodes[0]->cl[0][0][i][0] * aTree->model->parts[0]->comps[0]->val[i]);
  }
  printf("got likelihood=%f, %g\n", lnL, lnL);
  lnL = log(lnL);
  printf("got logLike=%f\n", lnL);
  lnL = 0.0;
#endif

#if 0
  printf("p4_treeLogLike()  check for bigQ's\n");
  for(i = 0; i < aTree->model->parts[0]->nComps; i++) {
    printf("    compNum %2i  bigQ %i\n", i, aTree->model->parts[0]->bigQAndEigThing[i][0]->bigQ);
  }
#endif

	

  for(pNum = 0; pNum < aTree->nParts; pNum++) {
    dp = aTree->data->parts[pNum];
    lnL += p4_partLogLike(aTree, dp, pNum, getSiteLikes);

  }
  //printf("p4_treeLogLike. got lnL = %f\n", lnL);
  aTree->logLike = lnL;
  return lnL;
}

double p4_partLogLike(p4_tree *aTree, part *dp, int pNum, int getSiteLikes)
{
	double lnL = 0.0;
	double like = 0.0;
	//double rateLike = 0.0;
	double theSum = 0.0;
	int    i, seqPos, rate, cNum, rNum;
	double *patternLikes = NULL;
	double  oneMinusPInvar = 1.0;
	p4_modelPart   *mp = aTree->model->parts[pNum];


	//printf("p4_partLogLike(). rootComp is\n");
	//printf("    %f  %f  %f  %f\n", mp->comps[0]->val[0], mp->comps[0]->val[1], mp->comps[0]->val[2], mp->comps[0]->val[3]);

	if(getSiteLikes){
		if(dp->siteLikes == NULL) {
			dp->siteLikes = malloc(dp->nChar * sizeof(double));
			if(!dp->siteLikes) {
				printf("Failed to malloc siteLikes.\n");
				exit(0);
			}
		}
		patternLikes = malloc(dp->nPatterns * sizeof(double));
		if(!patternLikes) {
			printf("Failed to malloc patternLikes.\n");
			exit(0);
		}
		for(seqPos = 0; seqPos < dp->nChar; seqPos++) {
			dp->siteLikes[seqPos] = 0.0;
		}
	}

#if 0
	printf("Conditional Likes of the root\n");
	for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
		printf("    Pattern position %i\n", seqPos);
		for(rate = 0; rate < mp->nCat; rate++) {
			printf("            Rate %i: ", rate);
			for(i = 0; i < mp->dim; i++) {
				printf("  %g", aTree->root->cl[pNum][rate][i][seqPos]);
			}
			printf("\n");
		}
	}
#endif

	oneMinusPInvar = 1.0 - mp->pInvar->val[0];

	if(!mp->freqsTimesOneMinusPInvar) {
		mp->freqsTimesOneMinusPInvar = malloc(mp->nCat * sizeof(double));
		if(!mp->freqsTimesOneMinusPInvar) {
			printf("memory allocation error, modelpart->freqsTimesOneMinusPInvar\n");
			exit(0);
		}
	}

	if(mp->isMixture) {
		for(rate = 0; rate < mp->nCat; rate++) {
			mp->freqsTimesOneMinusPInvar[rate] = oneMinusPInvar * mp->mixture->freqs[rate];
		}
	}
	else {
		for(rate = 0; rate < mp->nCat; rate++) {
			//mp->freqsTimesOneMinusPInvar[rate] = 
			//	aTree->root->models[pNum]->gammaFreqs[rate] * oneMinusPInvar;
			mp->freqsTimesOneMinusPInvar[rate] = oneMinusPInvar / (double)(mp->nCat);
		}
	}

	for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
		like = 0.0;

		// Either do pInvar or not
		if(mp->pInvar->val[0]) { // do pInvar
			// first deal with the contribution due to the possibility that it is a variable site

			if(mp->isMixture) {
				rate = 0;
				for(cNum = 0; cNum < mp->nComps; cNum++) {
					for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
						theSum = 0.0;
						for(i = 0; i < mp->dim; i++) {
							theSum += mp->comps[cNum]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
						}
						theSum *= mp->freqsTimesOneMinusPInvar[rate];
						like += theSum;
						rate++;
					}
				}
			}
			else {
				for(rate = 0; rate < mp->nCat; rate++) {
					for(i = 0; i < mp->dim; i++) {
						like += mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
#if 0
						printf("a comp[%i]=%f, cl[%i]=%f;   comp*cl=%f * %f = %f; cumLike=%f\n", 
							   i,
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   i,
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos],
							   like
							   );
#endif
					}
				}
#if 0
				printf("scaling by pVariable=%f.  %f * %f = %f\n", 
					   mp->freqsTimesOneMinusPInvar[0], 
					   like, mp->freqsTimesOneMinusPInvar[0], 
					   like * mp->freqsTimesOneMinusPInvar[0]);
#endif
				like *=  mp->freqsTimesOneMinusPInvar[0];
			}
			// Now deal with the invariant site contribution
			if(dp->globalInvarSitesVec[seqPos] > 0) { // ie its an invar site
				//printf("treeLogLike: invarSitesVec[%i] = %i\n", seqPos, dp->globalInvarSitesVec[seqPos]);				
				//printf("treeLogLike: doing invarSite contribution\n");
				// See part.c for an explanation of the magic of the globalInvarSitesArray
				// The globalInvarSitesArray only goes to dp->dim, but mp->dim might be bigger if covarion.
				if(mp->isMixture) {
					for(i = 0; i < dp->dim; i++) {
						if(dp->globalInvarSitesArray[i][seqPos]) {
							// add up the comp->val[i] over all mixture components
							rate = 0;
							theSum = 0.0;
							for(cNum = 0; cNum < mp->nComps; cNum++) {
								for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
									theSum += mp->comps[cNum]->val[i] * mp->mixture->freqs[rate];
									rate++;
								}
							}
							like += theSum * mp->pInvar->val[0];
						}
					}
				}
				else {
					for(i = 0; i < dp->dim; i++) {
						if(dp->globalInvarSitesArray[i][seqPos]) {
							like += mp->comps[aTree->root->compNums[pNum]]->val[i] * mp->pInvar->val[0];
#if 0
							printf("b comp[%i]=%f, pInvar=%f;   comp*pInvar=%f * %f = %f; cumLike=%f\n", 
								   i,
								   mp->comps[aTree->root->compNums[pNum]]->val[i],
								   mp->pInvar->val[0],
								   mp->comps[aTree->root->compNums[pNum]]->val[i],
								   mp->pInvar->val[0],
								   mp->comps[aTree->root->compNums[pNum]]->val[i] * mp->pInvar->val[0],
								   like
								   );
#endif

						}
					}
				}
			}
		}


		else {  // pInvar is not part of the equation
			if(mp->isMixture) {
				rate = 0;
				for(cNum = 0; cNum < mp->nComps; cNum++) {
					for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
						theSum = 0.0;
						for(i = 0; i < mp->dim; i++) {
							theSum += mp->comps[cNum]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
						}
						theSum *= mp->mixture->freqs[rate];
						like += theSum;
						rate++;
					}
				}
				
			}
			else {
				for(rate = 0; rate < mp->nCat; rate++) {
					//rateLike = 0.0;
					for(i = 0; i < mp->dim; i++) {
						like += mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
						//rateLike += mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
#if 0
						printf("d comp[%i]=%f, cl[%i]=%f;   comp*cl=%f * %f = %f; cumLike=%f\n", 
							   i,
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   i,
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos],
							   like
							   );
#endif
					}
					//printf("seqPos %3i rate %i, rateLike=%.12f\n", seqPos, rate, rateLike);
				}
				// If gamma freqs are not the same, it should be something like this:
				//like = like * aTree->model->parts[pNum]->gammaFreqs[0];
				// But the following will do if we have equal gamma freqs.
				if(mp->nCat > 1) {
					like = like / (double)(mp->nCat);
					//printf("seqPos %3i   like=%.12f\n", seqPos, like);
				}
			}
		}

		if(like <= 0.0) {
			//printf("p4_tree.c: treeLogLike: (zero-based) seqPos %i, site like %g\n", seqPos, like);
			//printf("    Its <= 0.0, so returning -1.0e99\n");
			//printf("    dp->globalInvarSitesVec[%i] = %i\n", seqPos, dp->globalInvarSitesVec[seqPos]);
			if(getSiteLikes) {
				if(patternLikes) {
					free(patternLikes);
					patternLikes = NULL;
					printf("Site likelihoods requested, but one likelihood is zero or less,\n");
					printf("    so getSiteLikes is not on.\n");
				}
			}

			return -1.0e99;
		}
		//printf("finished rate cats: seqPos = %i, lnL = %7.4f, like = %7.4f\n", seqPos, lnL, like);
		//printf("site=%i, like=%f, logLike=%f\n", seqPos, like, log(like));
		lnL = lnL + (dp->patternCounts[seqPos] * log(like));
		//printf("  seqPos %i   patCount %i  like %f  logLike %f  timesPatCount %f     total lnL %f\n",
		//		seqPos, dp->patternCounts[seqPos], like, log(like), 
		//		dp->patternCounts[seqPos] * log(like), lnL);
		
		if(getSiteLikes) patternLikes[seqPos] = like;
	}

	if(0) {
		if(getSiteLikes) {
			printf("tree.c: treeLogLike: patternLikes\n");
			for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
				printf("%f\n", patternLikes[seqPos]);
			}
		}
	}

	if(getSiteLikes) {
		// I have the pattern likes, so now figure out the siteLikes, using sequencePositionPatternIndex
		for(seqPos = 0; seqPos < dp->nChar; seqPos++) {
			dp->siteLikes[seqPos] = patternLikes[dp->sequencePositionPatternIndex[seqPos]];
		}
		free(patternLikes);
	}

	//printf("p4_tree.p4_partLogLike().  %.6f\n", lnL);
	aTree->partLikes[pNum] = lnL;
	return lnL;

}

//double p4_partLogLikeWinningGammaCats(p4_tree *aTree, part *dp, int pNum, int getSiteLikes, int *winningGammaCats, double *work)

double p4_partLogLikeSiteRates(p4_tree *aTree, part *dp, int pNum, int getSiteLikes, double *siteRates, int *gammaCats, double *work)
{
	double lnL = 0.0;
	double like = 0.0;
	double theSum = 0.0;
	int    i, seqPos, rate, cNum, rNum;
	double *patternLikes = NULL;
	double  oneMinusPInvar = 1.0;
	double  temp;
	double *tempRates = NULL;
	double  biggest;             // for gammaCats
	int    *tempCategories = NULL;  // for winning patternPositions, not actual sites.
	p4_modelPart   *mp = aTree->model->parts[pNum];


	if(getSiteLikes){
		if(dp->siteLikes == NULL) {
			dp->siteLikes = malloc(dp->nChar * sizeof(double));
			if(!dp->siteLikes) {
				printf("Failed to malloc siteLikes.\n");
				exit(0);
			}
		}
		patternLikes = malloc(dp->nPatterns * sizeof(double));
		if(!patternLikes) {
			printf("Failed to malloc patternLikes.\n");
			exit(0);
		}
		for(seqPos = 0; seqPos < dp->nChar; seqPos++) {
			dp->siteLikes[seqPos] = 0.0;
		}
	}

	if(mp->nGdasrvs != 1) {
		printf("nGdasrvs is %i, should be 1.\n", mp->nGdasrvs);
		exit(0);
	}

	tempCategories = malloc(dp->nPatterns * sizeof(int));
	if(!tempCategories) {
		printf("Failed to malloc tempCategories.\n");
		exit(0);
	}
	tempRates = malloc(dp->nPatterns * sizeof(double));
	if(!tempRates) {
		printf("Failed to malloc tempRates.\n");
		exit(0);
	}

#if 0
	printf("Conditional Likes of the root\n");
	for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
		printf("    Pattern position %i\n", seqPos);
		for(rate = 0; rate < mp->nCat; rate++) {
			printf("            Rate %i: ", rate);
			for(i = 0; i < mp->dim; i++) {
				printf("  %g", aTree->root->cl[pNum][rate][i][seqPos]);
			}
			printf("\n");
		}
	}
#endif

	oneMinusPInvar = 1.0 - mp->pInvar->val[0];

	if(!mp->freqsTimesOneMinusPInvar) {
		mp->freqsTimesOneMinusPInvar = malloc(mp->nCat * sizeof(double));
		if(!mp->freqsTimesOneMinusPInvar) {
			printf("memory allocation error, modelpart->freqsTimesOneMinusPInvar\n");
			exit(0);
		}
	}

	if(mp->isMixture) {
		for(rate = 0; rate < mp->nCat; rate++) {
			mp->freqsTimesOneMinusPInvar[rate] = oneMinusPInvar * mp->mixture->freqs[rate];
		}
	}
	else {
		for(rate = 0; rate < mp->nCat; rate++) {
			//mp->freqsTimesOneMinusPInvar[rate] = 
			//	aTree->root->models[pNum]->gammaFreqs[rate] * oneMinusPInvar;
			mp->freqsTimesOneMinusPInvar[rate] = oneMinusPInvar / (double)(mp->nCat);
		}
	}

	for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
		like = 0.0;

		// Either do pInvar or not
		if(mp->pInvar->val[0]) { // do pInvar
			// first deal with the contribution due to the possibility that it is a variable site

			if(mp->isMixture) {
				rate = 0;
				for(cNum = 0; cNum < mp->nComps; cNum++) {
					for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
						theSum = 0.0;
						for(i = 0; i < mp->dim; i++) {
							theSum += mp->comps[cNum]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
						}
						theSum *= mp->freqsTimesOneMinusPInvar[rate];
						like += theSum;
						rate++;
					}
				}
			}
			else {
				for(rate = 0; rate < mp->nCat; rate++) {
					work[rate] = 0.0;
					for(i = 0; i < mp->dim; i++) {
						like += mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
						work[rate] += mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
#if 0
						printf("a comp[%i]=%f, cl[%i]=%f;   comp*cl=%f * %f = %f; cumLike=%f\n", 
							   i,
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   i,
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos],
							   like
							   );
#endif
					}
				}
#if 0
				printf("scaling by pVariable=%f.  %f * %f = %f\n", 
					   mp->freqsTimesOneMinusPInvar[0], 
					   like, mp->freqsTimesOneMinusPInvar[0], 
					   like * mp->freqsTimesOneMinusPInvar[0]);
#endif
				like *=  mp->freqsTimesOneMinusPInvar[0];
			}
			// Now deal with the invariant site contribution
			if(dp->globalInvarSitesVec[seqPos] > 0) { // ie its an invar site
				//printf("treeLogLike: invarSitesVec[%i] = %i\n", seqPos, dp->globalInvarSitesVec[seqPos]);				
				//printf("treeLogLike: doing invarSite contribution\n");
				// See part.c for an explanation of the magic of the globalInvarSitesArray
				// The globalInvarSitesArray only goes to dp->dim, but mp->dim might be bigger if covarion.
				if(mp->isMixture) {
					for(i = 0; i < dp->dim; i++) {
						if(dp->globalInvarSitesArray[i][seqPos]) {
							// add up the comp->val[i] over all mixture components
							rate = 0;
							theSum = 0.0;
							for(cNum = 0; cNum < mp->nComps; cNum++) {
								for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
									theSum += mp->comps[cNum]->val[i] * mp->mixture->freqs[rate];
									rate++;
								}
							}
							like += theSum * mp->pInvar->val[0];
						}
					}
				}
				else {
					for(i = 0; i < dp->dim; i++) {
						if(dp->globalInvarSitesArray[i][seqPos]) {
							like += mp->comps[aTree->root->compNums[pNum]]->val[i] * mp->pInvar->val[0];
#if 0
							printf("b comp[%i]=%f, pInvar=%f;   comp*pInvar=%f * %f = %f; cumLike=%f\n", 
								   i,
								   mp->comps[aTree->root->compNums[pNum]]->val[i],
								   mp->pInvar->val[0],
								   mp->comps[aTree->root->compNums[pNum]]->val[i],
								   mp->pInvar->val[0],
								   mp->comps[aTree->root->compNums[pNum]]->val[i] * mp->pInvar->val[0],
								   like
								   );
#endif

						}
					}
				}
			}
		}


		else {  // pInvar is not part of the equation
			if(mp->isMixture) {
				rate = 0;
				for(cNum = 0; cNum < mp->nComps; cNum++) {
					for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
						theSum = 0.0;
						for(i = 0; i < mp->dim; i++) {
							theSum += mp->comps[cNum]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
						}
						theSum *= mp->mixture->freqs[rate];
						like += theSum;
						rate++;
					}
				}
				
			}
			else {
				for(rate = 0; rate < mp->nCat; rate++) {
					work[rate] = 0.0;
					for(i = 0; i < mp->dim; i++) {
						like += mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
						work[rate] += mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos];
#if 0
						printf("d comp[%i]=%f, cl[%i]=%f;   comp*cl=%f * %f = %f; cumLike=%f\n", 
							   i,
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   i,
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i],
							   aTree->root->cl[pNum][rate][i][seqPos],
							   mp->comps[aTree->root->compNums[pNum]]->val[i] * aTree->root->cl[pNum][rate][i][seqPos],
							   like
							   );
#endif
					}
				}
				// If gamma freqs are not the same, it should be something like this:
				//like = like * aTree->model->parts[pNum]->gammaFreqs[0];
				// But the following will do if we have equal gamma freqs.
				if(mp->nCat > 1) {
					like = like / (double)(mp->nCat);
				}
			}
		}

		if(like <= 0.0) {
			//printf("p4_tree.c: treeLogLike: (zero-based) seqPos %i, site like %g\n", seqPos, like);
			//printf("    Its <= 0.0, so returning -1.0e99\n");
			//printf("    dp->globalInvarSitesVec[%i] = %i\n", seqPos, dp->globalInvarSitesVec[seqPos]);
			if(getSiteLikes) {
				if(patternLikes) {
					free(patternLikes);
					patternLikes = NULL;
					printf("Site likelihoods requested, but one likelihood is zero or less,\n");
					printf("    so getSiteLikes is not on.\n");
				}
			}

			return -1.0e99;
		}
		//printf("finished rate cats: seqPos = %i, lnL = %7.4f, like = %7.4f\n", seqPos, lnL, like);
		//printf("site=%i, like=%f, logLike=%f\n", seqPos, like, log(like));
		lnL = lnL + (dp->patternCounts[seqPos] * log(like));
		//printf("  seqPos %i   patCount %i  like %f  logLike %f  timesPatCount %f     total lnL %f\n",
		//		seqPos, dp->patternCounts[seqPos], like, log(like), 
		//		dp->patternCounts[seqPos] * log(like), lnL);
		
		if(getSiteLikes) patternLikes[seqPos] = like;

		biggest = work[0];
		tempCategories[seqPos] = 0;
		for(rate = 0; rate < mp->nCat; rate++) {
			if(work[rate] > biggest) {
				biggest = work[rate];
				tempCategories[seqPos] = rate;
			}
		}


		temp = 0.0;
		for(rate = 0; rate < mp->nCat; rate++) {
			temp += work[rate] * mp->gdasrvs[0]->freqs[rate] * mp->gdasrvs[0]->rates[rate];
		}
		tempRates[seqPos] = temp/like;

	}

	if(0) {
		if(getSiteLikes) {
			printf("tree.c: treeLogLike: patternLikes\n");
			for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
				printf("%f\n", patternLikes[seqPos]);
			}
		}
	}

	if(getSiteLikes) {
		// I have the pattern likes, so now figure out the siteLikes, using sequencePositionPatternIndex
		for(seqPos = 0; seqPos < dp->nChar; seqPos++) {
			dp->siteLikes[seqPos] = patternLikes[dp->sequencePositionPatternIndex[seqPos]];
		}
		free(patternLikes);
	}

	//for(seqPos = 0; seqPos < 10; seqPos++) {
	//	printf("%2i  seqPositionPatternIndex %3i\n", seqPos, dp->sequencePositionPatternIndex[seqPos]);
	//}


	for(seqPos = 0; seqPos < dp->nChar; seqPos++) {
		siteRates[seqPos] = tempRates[dp->sequencePositionPatternIndex[seqPos]];
		gammaCats[seqPos] = tempCategories[dp->sequencePositionPatternIndex[seqPos]];
	}
	free(tempRates);
	free(tempCategories);
	

	//printf("p4_tree.p4_partLogLikeWinningGammaCats().  %f\n", lnL);
	aTree->partLikes[pNum] = lnL;
	return lnL;


}


void p4_getPreOrderNodeNumsAbove(p4_tree *aTree, p4_node *aNode)
{
	int stackIndx, intsIndx;
	p4_node *p, *q;

	//printf("    p4_getPreOrderNodeNumsAbove() here.\n");
	
	// If there are no nodes above, return early
	if(!aNode->leftChild) {
		aTree->ints[0] = -1;
		return;
	}

	aTree->stack[0] = aNode;
	stackIndx = 1; // next
	intsIndx = 0;  // next
	while(stackIndx) {
		p = aTree->stack[stackIndx - 1];
		if(p->leftChild) {
			aTree->stack[stackIndx] =  p->leftChild;
			stackIndx++;
			//printf("l preOrder appending %i\n", p->leftChild->nodeNum);
			aTree->ints[intsIndx] = p->leftChild->nodeNum;
			intsIndx++;
		}
		else if(p->sibling) {
			aTree->stack[stackIndx - 1] = p->sibling;
			//printf("s preOrder appending %i\n", p->sibling->nodeNum);
			aTree->ints[intsIndx] = p->sibling->nodeNum;
			intsIndx++;
		}
		else {
			stackIndx--;
			//printf("a stackIndx=%i\n", stackIndx);
			if(stackIndx == 0) {
				break;
			}
			q = aTree->stack[stackIndx - 1];
			//printf("q is node %i\n", q->nodeNum);
			stackIndx--;
			//printf("b stackIndx=%i\n", stackIndx);
			while(!(q->sibling)) {
				if(stackIndx == 0) {
					break;
				}
				q = aTree->stack[stackIndx - 1];
				//printf("q is node %i\n", q->nodeNum);
				stackIndx--;
				//printf("c stackIndx=%i\n", stackIndx);
			}
			if(stackIndx == 0) {
				break;
			}
			if(q->sibling) {
				//printf("d  q is node %i, its sib is node %i, stackIndx=%i, intsIndx=%i\n", 
				//	   q->nodeNum, q->sibling->nodeNum, stackIndx, intsIndx);
				aTree->stack[stackIndx] = q->sibling;
				stackIndx++;
				//printf("z preOrder appending %i\n", q->sibling->nodeNum);
				aTree->ints[intsIndx] = q->sibling->nodeNum;
				intsIndx++;
			}
		} 
	}
	aTree->ints[intsIndx] = -1;
	
}


