#include "pftypes.h"
#include "p4_tree.h"
#include "p4_node.h"
#include "pmatrices.h"
#include "defines.h"
#include "util.h"

//static int getFirstLike;
//static int newtCount;

void p4_newtSetup(p4_tree *aTree)
{
  int nNum, i, j;
  p4_node  *aNode;

  // malloc cl2, if needed.
  if(!aTree->root->cl2) {
    for(nNum = 0; nNum < aTree->nNodes; nNum++) {
      aNode = aTree->nodes[nNum];
      aNode->cl2 = (double ****)malloc(aNode->nParts * sizeof(double ***));
      if(!aNode->cl2) {
	printf("Failed to allocate memory for cl2.\n");
	exit(1);
      }
      for(i = 0; i < aNode->nParts; i++) {
	aNode->cl2[i] = (double ***)malloc(aNode->tree->model->parts[i]->nCat * sizeof(double **));
	if(!aNode->cl2[i]) {
	  printf("Failed to allocate memory for cl2[i].\n");
	  exit(1);
	}
	for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
	  aNode->cl2[i][j] = pdmatrix(aNode->tree->model->parts[i]->dim, aNode->tree->data->parts[i]->nChar);
	}
      }
    }
  }

  // malloc bigPDecks_1stD and bigPDecks_2ndD
  if(!aTree->root->bigPDecks_1stD) {
    for(nNum = 0; nNum < aTree->nNodes; nNum++) {
      aNode = aTree->nodes[nNum];

      aNode->bigPDecks_1stD = (double ****)malloc(aNode->nParts * sizeof(double ***));
      if(!aNode->bigPDecks_1stD) {
	printf("Failed to allocate memory for bigPDecks_1stD.\n");
	exit(1);
      }
      for(i = 0; i < aNode->nParts; i++) {
	aNode->bigPDecks_1stD[i] = (double ***)malloc(aNode->tree->model->parts[i]->nCat * sizeof(double **));
	if(!aNode->bigPDecks_1stD[i]) {
	  printf("Failed to allocate memory for bigPDecks_1stD[i].\n");
	  exit(1);
	}
	for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
	  aNode->bigPDecks_1stD[i][j] = psdmatrix(aNode->tree->model->parts[i]->dim);
	}
      }

      aNode->bigPDecks_2ndD = (double ****)malloc(aNode->nParts * sizeof(double ***));
      if(!aNode->bigPDecks_2ndD) {
	printf("Failed to allocate memory for bigPDecks_2ndD.\n");
	exit(1);
      }
      for(i = 0; i < aNode->nParts; i++) {
	aNode->bigPDecks_2ndD[i] = (double ***)malloc(aNode->tree->model->parts[i]->nCat * sizeof(double **));
	if(!aNode->bigPDecks_2ndD[i]) {
	  printf("Failed to allocate memory for bigPDecks_2ndD[i].\n");
	  exit(1);
	}
	for(j = 0; j < aNode->tree->model->parts[i]->nCat; j++) {
	  aNode->bigPDecks_2ndD[i][j] = psdmatrix(aNode->tree->model->parts[i]->dim);
	}
      }
    }
  }
}


void p4_newtAround(p4_tree *aTree, double epsilon, double likeDelta)
{
  int nNum, i, j;
  double thisLogLike, previousLogLike, diff;
  p4_node *aNode, *p;

  previousLogLike = p4_treeLogLike(aTree, 0);
  //printf("Starting p4_newtAround with logLike=%f, epsilon=%f, likeDelta=%g\n", previousLogLike, epsilon, likeDelta);
  //getFirstLike = 1;
  //newtCount = 0;
  for(i = 0; i < 20; i++) {                   // default 20, raise for difficult opts.
    // At first, set all cl2NeedsUpdating
    for(nNum = 0; nNum < aTree->nNodes; nNum++) {
      aNode = aTree->nodes[aTree->postOrder[nNum]];
      if(aNode != aTree->root) {
	aNode->cl2NeedsUpdating = 1;
      }
    }

    for(nNum = 0; nNum < aTree->nNodes; nNum++) {
      aNode = aTree->nodes[aTree->postOrder[nNum]];
      if(aNode != aTree->root) {
	//printf("\n=== do node %i\n", aNode->nodeNum);

	// Update cl, if needed.  It will only be one node
	// that needs it.  Thats cuz we are going in
	// post-order.
	if(aNode->clNeedsUpdating) {
	  //printf("    update cl on node %i\n", aNode->nodeNum);
	  p4_setConditionalLikelihoodsOfInteriorNode(aNode);
	}
				
	// Update cl2, if needed.  Need to do any nodes with
	// cl2NeedsUpdating set, going from aNode down to the
	// root.  The path from aNode down to the root can be
	// put in aTree->ints, a utility vector.
	if(aNode->cl2NeedsUpdating) {
	  p = aNode;
	  j = 0; // index for aTree->ints
	  while(p->parent && p->cl2NeedsUpdating) {
	    aTree->ints[j] = p->nodeNum;
	    p = p->parent;
	    j++;
	  }
	  aTree->ints[j] = -1; // end of the path
	  //printf("    got path = ");
	  for(j = 0; j < aTree->nNodes; j++) {
	    if(aTree->ints[j] >= 0) {
	      //printf(" %i", aTree->ints[j]);
	    }
	    else {
	      break;
	    }
	  }
	  //printf("\n");
	  // we need j on the -1
	  if(aTree->ints[j] != -1) {
	    printf("Bad.\n");
	    exit(1);
	  }
	  while(j) {
	    j--;
	    p = aTree->nodes[aTree->ints[j]];
	    //printf("    update cl2 for node %i\n", p->nodeNum);
	    p4_setNodeCL2(aTree, p);
	  }
	}

	// Newt it
	//printf("    Now newtNode %i\n", aNode->nodeNum);
	p4_newtNode(aNode, epsilon);

	// Flag nodes above aNode->parent that cl2NeedsUpdating
	p4_getPreOrderNodeNumsAbove(aTree, aNode); // gets put into aTree->ints
	//printf("    got node nums above = ");
	for(j = 0; j < aTree->nNodes; j++) {
	  if(aTree->ints[j] >= 0) {
	    //printf(" %i", aTree->ints[j]);
	  }
	  else {
	    break;
	  }
	}
	//printf("\n");
				

	// Flag nodes below aNode that cl clNeedsUpdating
	p = aNode->parent;
	while(p) {
	  p->clNeedsUpdating = 1;
	  p = p->parent;
	}
      }
    }

    thisLogLike = p4_treeLogLike(aTree, 0);
    diff = thisLogLike - previousLogLike;
#if 0
    printf("%i end p4_newtAround() logLike = %f, diff=%f\n", i, thisLogLike, diff);
    //exit(0);
#endif
#if 0
    printf("%i p4_newtAround() logLike = %f, diff=%f, newtCount=%i\n", i, thisLogLike, diff, newtCount);
    //exit(0);
#endif
    if(diff < -1.0e-5) {  // changed from -1.0e-6, which seemed to be too sensitive.
      printf("Bad newt.  Likelihood decreased.  Diff=%f (%g) ...continuing anyway\n", diff, diff);
      //for(nNum = 0; nNum < aTree->nNodes; nNum++) {
      //	aNode = aTree->nodes[nNum];
      //	if(aNode != aTree->root) {
      //		printf("        %2i    %.18f (%g)\n", aNode->nodeNum, aNode->brLen[0], aNode->brLen[0]);
      //	}
      //}
			
      //exit(0);
    }
    if(fabs(diff) < likeDelta) {
      //if(1) {
      //printf("finished p4_newtAround() did %i loops, logLike = %f, diff=%f\n", i + 1, thisLogLike, diff);
      return;
    }
    previousLogLike = thisLogLike;
  }
  //printf("p4_newtAround.  max iterations reached. i=%i, logLike=%f, diff=%f\n", i, thisLogLike, diff);
  //return newtCount;
}


void p4_newtNode(p4_node *aNode, double epsilon)
{
  double oldBrLen, currentGuess, previousGuess, nextGuess;
  double lnL, firstD, secondD;
  double like, first, second;    // for each rate
  double likeS, firstS, secondS; // S for site
  //double oneMinusPInvar;
  double temp2;
  int from, to, charCode;
  int  iter, pNum, rate, seqPos;
  p4_modelPart   *mp;
  part           *dp;
  int isN;
  int cNum, rNum;
  //double startingLike;

  oldBrLen = aNode->brLen[0];
  currentGuess = oldBrLen;
  previousGuess = oldBrLen;
  nextGuess = 0.0;
  lnL = 0.0;
  iter = 0;
	
  like = first = second = 0.0;

#if 0
  // only for debugging
  startingLike = p4_treeLogLike(aNode->tree, 0);
#endif

  firstD = 1.0;
  while(1) {
    //newtCount++;
    p4_calculateBigPDecks(aNode);
    p4_calculateBigPDecks_1stD(aNode);
    p4_calculateBigPDecks_2ndD(aNode);

#if 0
    printf("bigP, node %i\n", aNode->nodeNum);
    dump_psdmatrix(aNode->bigPDecks[0][0], 4);
    printf("bigP_1stD, node %i\n", aNode->nodeNum);
    dump_psdmatrix(aNode->bigPDecks_1stD[0][0], 4);    
    printf("bigP_2ndD, node %i\n", aNode->nodeNum);
    dump_psdmatrix(aNode->bigPDecks_2ndD[0][0], 4);
#endif

    lnL = 0.0;
    firstD = 0.0;
    secondD = 0.0;

    for(pNum = 0; pNum < aNode->nParts; pNum++) {
      //printf("    pNum %i\n", pNum);
      mp = aNode->tree->model->parts[pNum];
      dp = aNode->tree->data->parts[pNum];

      //oneMinusPInvar = 1.0 - mp->pInvar->val[0];

      if(!mp->freqsTimesOneMinusPInvar) {
	mp->freqsTimesOneMinusPInvar = malloc(mp->nCat * sizeof(double));
	if(!mp->freqsTimesOneMinusPInvar) {
	  printf("memory allocation error, part->freqsTimesOneMinusPInvar\n");
	  exit(0);
	}
      }
      if(mp->isMixture) {
	for(rate = 0; rate < mp->nCat; rate++) {
	  mp->freqsTimesOneMinusPInvar[rate] = (1.0 - mp->pInvar->val[0]) *  mp->mixture->freqs[rate];
	}
      }
      else {
	for(rate = 0; rate < mp->nCat; rate++) {
	  mp->freqsTimesOneMinusPInvar[rate] = (1.0 - mp->pInvar->val[0]) / (double)(mp->nCat);
	}
      }
			
      for(seqPos = 0; seqPos < dp->nPatterns; seqPos++) {
	//printf("        seqPos %i\n", seqPos);
	likeS = firstS = secondS = 0.0;

	// Either do pInvar or not
	if(mp->pInvar->val[0]) { // do pInvar
	  // First deal with the contribution due to the possibility that it is a variable site
	  for(rate = 0; rate < mp->nCat; rate++) {
	    like = first = second = 0.0;
	    if(aNode->isLeaf) {
	      charCode = dp->patterns[aNode->seqNum][seqPos];
	      if(charCode >= 0) {
		if(mp->doTSCovarion) {
		  // not done
		}
		else {
		  for(from = 0; from < mp->dim; from++) {
		    // Using temp2 speeds things up a tiny, tiny bit.
		    temp2 = aNode->cl2[pNum][rate][from][seqPos];
		    like += temp2 * aNode->bigPDecks[pNum][rate][from][charCode];
		    first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][charCode];
		    second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][charCode];
		  }
		}
	      }

	      else if (charCode == GAP_CODE || charCode == QMARK_CODE) {
		for(from = 0; from < mp->dim; from++) {
		  temp2 = aNode->cl2[pNum][rate][from][seqPos];
		  like += temp2;
		  // If it was just likes, I would
		  // not have to do this following
		  // loop, cuz sum of a bigP row is
		  // 1.0.  However, that is not so
		  // for the 1stD.  Or 2ndD?
		  for(to = 0; to < mp->dim; to++) {
		    //like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
		    first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
		    second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		  }
		}
	      }
	      else if(charCode >= EQUATES_BASE && 
		      charCode < EQUATES_BASE + dp->nEquates) {
		// If we get this to this point, most
		// usually, it will be an "n" or n-like
		// character (eg x).  Assume it is so unless
		// datapart->equates says otherwise.
		isN = 1;
		for(from = 0; from < dp->dim; from++) {
		  if(!dp->equates[charCode - EQUATES_BASE][from]) {
		    isN = 0;
		    break;
		  }
		}
		//printf("isN = %i\n", isN);
		if(isN) {
		  for(from = 0; from < mp->dim; from++) {
		    temp2 = aNode->cl2[pNum][rate][from][seqPos];
		    for(to = 0; to < mp->dim; to++) {
		      like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
		      first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
		      second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		    }
		  }
		}
		else {
		  for(from = 0; from < mp->dim; from++) {
		    temp2 = aNode->cl2[pNum][rate][from][seqPos];
		    for(to = 0; to < dp->dim; to++) {
		      if(dp->equates[charCode - EQUATES_BASE][to]) {
			like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
			first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
			second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		      }
		    }
		  }
		}
	      }



	      else {
		printf("p4_newtNode(%i):\n", aNode->nodeNum);
		printf("   Programming error.  This shouldn't happen\n");
		exit(0);
	      }

	    }
	    else {  // its not a leaf, use cond likes
	      for(from = 0; from < mp->dim; from++) {
		for(to = 0; to < mp->dim; to++) {
		  temp2 = aNode->cl2[pNum][rate][from][seqPos] * aNode->cl[pNum][rate][to][seqPos];
		  like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
		  first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
		  second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		}
	      }
							
	    }
						
	    // If it is a mixture, now we scale like by freq[rate]
	    if(mp->isMixture) {
	      //if(0) {
	      like   *= mp->freqsTimesOneMinusPInvar[rate];
	      first  *= mp->freqsTimesOneMinusPInvar[rate];
	      second *= mp->freqsTimesOneMinusPInvar[rate];
	    }
						
	    likeS += like;
	    firstS += first;
	    secondS += second;
						
	  }  // for(rate)
	  if(!mp->isMixture) {
	    likeS   *= mp->freqsTimesOneMinusPInvar[0];
	    firstS  *= mp->freqsTimesOneMinusPInvar[0];
	    secondS *= mp->freqsTimesOneMinusPInvar[0];
	  }

	  // Now deal with the invariant site contribution
	  if(dp->globalInvarSitesVec[seqPos] > 0) {  // ie its an invar site
	    //printf("a non-varied sited.\n");
	    if(mp->isMixture) {
	      for(from = 0; from < mp->dim; from++) {
		if(dp->globalInvarSitesArray[from][seqPos]) {
		  rate = 0;
		  temp2 = 0.0;
		  for(cNum = 0; cNum < mp->nComps; cNum++) {
		    for(rNum = 0; rNum < mp->nRMatrices; rNum++) {
		      temp2 += mp->comps[cNum]->val[from] * mp->mixture->freqs[rate];
		      rate++;
		    }
		  }
		  likeS += temp2 * mp->pInvar->val[0];
		}
	      }
	    }
	    else {
	      for(from = 0; from < mp->dim; from++) {
		if(dp->globalInvarSitesArray[from][seqPos]) {
		  likeS += mp->comps[aNode->tree->root->compNums[pNum]]->val[from] * mp->pInvar->val[0];
		  // first and second are not affected
		}
	      }
	    }
	  }

					
	}
	else { // pInvar does not apply
	  for(rate = 0; rate < mp->nCat; rate++) {
	    like = first = second = 0.0;
	    if(aNode->isLeaf) {
	      charCode = dp->patterns[aNode->seqNum][seqPos];
	      if(charCode >= 0) {
		if(mp->doTSCovarion) {
		  printf("check me.\n");
		  for(from = 0; from < mp->dim; from++) {
		    temp2 = aNode->cl2[pNum][rate][from][seqPos];
		    like += temp2 * (aNode->bigPDecks[pNum][rate][from][charCode] +
				     aNode->bigPDecks[pNum][rate][from][charCode + dp->dim]);
		    first += temp2 * (aNode->bigPDecks_1stD[pNum][rate][from][charCode] +
				      aNode->bigPDecks_1stD[pNum][rate][from][charCode + dp->dim]);
		    second += temp2 * (aNode->bigPDecks_2ndD[pNum][rate][from][charCode] +
				       aNode->bigPDecks_2ndD[pNum][rate][from][charCode + dp->dim]);
		  }

		}
		else {
		  for(from = 0; from < mp->dim; from++) {
		    // Using temp2 speeds things up a tiny, tiny bit.
		    temp2 = aNode->cl2[pNum][rate][from][seqPos];
		    //temp2 *= mp->mixture->freqs[rate];
		    like += temp2 * aNode->bigPDecks[pNum][rate][from][charCode];
		    first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][charCode];
		    second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][charCode];
		  }
		}
	      }

	      else if (charCode == GAP_CODE || charCode == QMARK_CODE) {
		for(from = 0; from < mp->dim; from++) {
		  temp2 = aNode->cl2[pNum][rate][from][seqPos];
		  like += temp2;
		  // If it was just likes, I would
		  // not have to do this following
		  // loop, cuz sum of a bigP row is
		  // 1.0.  However, that is not so
		  // for the 1stD.  Or 2ndD?
		  for(to = 0; to < mp->dim; to++) {
		    //like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
		    first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
		    second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		  }
		}
	      }
	      else if(charCode >= EQUATES_BASE && 
		      charCode < EQUATES_BASE + dp->nEquates) {
		// If we get this to this point, most
		// usually, it will be an "n" or n-like
		// character (eg x).  Assume it is so unless
		// datapart->equates says otherwise.
		isN = 1;
		for(from = 0; from < dp->dim; from++) {
		  if(!dp->equates[charCode - EQUATES_BASE][from]) {
		    isN = 0;
		    break;
		  }
		}
		//printf("isN = %i\n", isN);
		if(isN) {
		  for(from = 0; from < mp->dim; from++) {
		    temp2 = aNode->cl2[pNum][rate][from][seqPos];
		    for(to = 0; to < mp->dim; to++) {
		      like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
		      first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
		      second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		    }
		  }
		}
		else {
		  for(from = 0; from < mp->dim; from++) {
		    temp2 = aNode->cl2[pNum][rate][from][seqPos];
		    for(to = 0; to < dp->dim; to++) {
		      if(dp->equates[charCode - EQUATES_BASE][to]) {
			like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
			first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
			second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		      }
		    }
		  }
		}
	      }



	      else {
		printf("p4_newtNode(%i):\n", aNode->nodeNum);
		printf("   Programming error.  This shouldn't happen\n");
		exit(0);
	      }

	    }  // if(aNode->isLeaf)
	    else {  // aNode is not aLeaf, so use aNode->cl
	      for(from = 0; from < mp->dim; from++) {
		for(to = 0; to < mp->dim; to++) {
		  temp2 = aNode->cl2[pNum][rate][from][seqPos] * aNode->cl[pNum][rate][to][seqPos];
		  //temp2 *= mp->mixture->freqs[rate];
		  like += temp2 * aNode->bigPDecks[pNum][rate][from][to];
		  first += temp2 * aNode->bigPDecks_1stD[pNum][rate][from][to];
		  second += temp2 * aNode->bigPDecks_2ndD[pNum][rate][from][to];
		}
	      }
	    }
						
	    // If it is a mixture, now we scale like by freq[rate]
	    if(mp->isMixture) {
	      //if(0) {
	      like *= mp->mixture->freqs[rate];
	      first *= mp->mixture->freqs[rate];
	      second *= mp->mixture->freqs[rate];
	    }
						
	    likeS += like;
	    firstS += first;
	    secondS += second;
						
	  } // for(rate)
	  if(!mp->isMixture) {
	    // If it is not a mixture, we can do this outside the for(rate) loop.
	    // If gamma freqs are not the same, it should be something like this:
	    // like *= aTree->model->parts[pNum]->gammaFreqs[0];
	    // But the following will do since we have equal gamma freqs.
	    if(mp->nCat > 1) {
	      likeS   /=  (double)(mp->nCat);
	      firstS  /=  (double)(mp->nCat);
	      secondS /=  (double)(mp->nCat);
	    }
	  }
	} // else pInvar does not apply
	if(likeS < 1.0e-300) {
	  lnL += dp->patternCounts[seqPos] *      -100000;
	  firstD += dp->patternCounts[seqPos] *   1000000;
	  secondD += dp->patternCounts[seqPos] * 10000000;
	}
	else {
	  //lnL += dp->patternCounts[seqPos] * log(like);  // not really needed
	  //firstD += dp->patternCounts[seqPos] * (first / like);
	  //secondD += dp->patternCounts[seqPos] * ((second * like - first * first) / (like * like));
	  lnL += dp->patternCounts[seqPos] * log(likeS);  // not really needed
	  firstD += dp->patternCounts[seqPos] * (firstS / likeS);
	  secondD += dp->patternCounts[seqPos] * ((secondS * likeS - firstS * firstS) / (likeS * likeS));
	}
      }  // for(seqPos)
				
    } // for(pNum)
#if 0
    printf("x    end newtNode(nodeNum %i)  logLike=%f \n", aNode->nodeNum, lnL);
    //p4_dumpCL(aNode->cl);
    //p4_dumpCL(aNode->cl2);
    //break;
    exit(0);
#endif

#if 0
    if(getFirstLike) {
      printf("        p4_newtNode(%i) firstLike=%f\n", aNode->nodeNum, lnL);
      if(fabs(lnL - startingLike) > 1.0e-6) {
	printf("        p4_newtNode(%i) startingLike=%f, firstLike=%f\n", aNode->nodeNum, startingLike, lnL);
	exit(0);
      }
      getFirstLike = 0;
    }
#endif
			
    //printf("currentGuess = %f, logLike= %f\n", currentGuess, lnL);
    nextGuess = currentGuess - (firstD / secondD);
    if(secondD >= 0.0) {
      nextGuess = currentGuess / 5.0;
    }

    // break if the brLen is too short
    if(nextGuess <= BRLEN_MIN) {
      nextGuess = previousGuess / 4.0;
      if(nextGuess <= BRLEN_MIN) {
	aNode->brLen[0] = BRLEN_MIN - (BRLEN_MIN * 0.1) + ((BRLEN_MIN * 0.1) * ranDoubleUpToOne());
	break;
      }
    }

    // We break and return if the brLen is too long, or if we have
    // iterated more than the max, or if the likelihood curve is
    // flat enough.
    if(nextGuess >= 5.0 * oldBrLen) {
      aNode->brLen[0] = 5.0 * oldBrLen;
      break;
    }
    if(nextGuess > BRLEN_MAX) {
      aNode->brLen[0] = BRLEN_MAX;
      break;
    }
    if(iter > 20) {
      aNode->brLen[0] = currentGuess;
      //printf("newtNode %i, exceeded iter.\n", aNode->nodeNum);
      break;
    }
    iter++;

    if(fabs(firstD) < epsilon) {
      aNode->brLen[0] = currentGuess;
      break;
    }

    previousGuess = currentGuess;
    currentGuess = nextGuess;
    aNode->brLen[0] = nextGuess;
		
  }
  // After a break, we need to make sure that bigP is up to date, as
  // it is needed for cl calcs and for brent.  We do not need to
  // worry about bigP_1stD and bigP_2ndD.
  p4_calculateBigPDecks(aNode);
  //printf("      newtNode finishing. node %i, did %i loops, brLen=%f\n", aNode->nodeNum, iter + 1, aNode->brLen[0]);
  //printf(" %i", iter + 1);
}

void p4_setNodeCL2(p4_tree *aTree, p4_node *aNode)
{
  p4_node *p;

  if(aNode->parent == aTree->root) {
    //printf("      initializing node %i\n", aNode->nodeNum);
    p4_initializeCL2ToRootComp(aNode);
  }
  else {
    //printf("      do node %i going up, new on node %i\n", aNode->parent->nodeNum, aNode->nodeNum);
    p4_setCL2Up(aNode);
  }
	
  p = aNode->parent->leftChild;
  while(p) {
    if(p != aNode) {
      //printf("      do node %i, going down, combining with node %i\n", p->nodeNum, aNode->nodeNum);
      p4_setCL2Down(aNode, p);
    }
    p = p->sibling;
  }

  aNode->cl2NeedsUpdating = 0;
}




