#include "pftypes.h"
#include "p4_tree.h"
#include "pmatrices.h"
#include "defines.h"
// #include "util.h"
#include "p4_node.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>




void p4_simulate(p4_tree *t, p4_tree *refTree, const gsl_rng  *g)
{

    double  theRanDoubleUpToOne;
    int i, j, k, l;
    int pNum, nNum, seqPos;
    int simCat, parentState;
    data  *d;
    part  *dp;
    int   internSeqCounter = 0, tOrd = 0;
    p4_node *n;
    p4_modelPart  *mp;
	

    //printf("p4_simulate() here.\n");
    //return;

    /*


      Each part has space for sequences, which will be fine for the
      terminal nodes.  But they do not have enough room for the
      internal nodes, needed for the simulation.  So here is what I
      do: I make t->simSequences, with one t->nNodes by
      t->data->parts[i]->nChar integer matrix for each part, but
      not using pimatrix to make it, nor free_pimatrix to free it.  It
      is a ***int as usual.  For each part, it has nNodes int vectors,
      each of which are nChar long.  Here is the tricky part- the
      terminal nodes get the int vectors from the part->sequences, and
      the internal nodes get the int vectors from t->internalSequences.

      But the advantage of it all is that after a simulation the
      terminal node sequences are already in the part sequences array,
      and I can then do things like make patterns and calculate
      likelihoods.

      The picker is a vector of doubles.  It is the cumulative sum of
      the charFreq or the gammaFreqs eg for a charFreq = (0.4, 0.3, 0.2,
      0.1), the picker would be (0.4, 0.7, 0.9, 1.0).  So if I pick a
      random float between zero and 1.0, I can use the picker to pick
      a base that gives the composition of charFreq.

      The picker deck works the same way, the running sum of bigP

    */ 



    // Since we are over-writing t->data->parts[i]->sequences, the 
    // t->part->patterns is no longer representative.  My signal
    // that that is so is to zero the patternCount.

    d = t->data;
    for(pNum = 0; pNum < d->nParts; pNum++) {
        //printf("    part %i, previous nPatterns = %i\n", i, d->parts[pNum]->nPatterns);
        d->parts[pNum]->nPatterns = 0;
    }


    // If we have not allocated simSequences yet, do so now.
    if(t->simSequences == NULL) {
        t->simSequences = malloc(d->nParts * sizeof(int **));
        if(!t->simSequences) {
            printf("Problem allocating memory for simSequences (1)\n");
            exit(1);
        }
        for(pNum = 0; pNum < d->nParts; pNum++) {
            t->simSequences[pNum] = malloc(t->nNodes * sizeof(int *));
            if(!t->simSequences[pNum]) {
                printf("Problem allocating memory for simSequences (2)\n");
                exit(1);
            }
        }
    }

    // If we have not allocated internalSequences yet, do so now.
    if(t->internalSequences == NULL) {
        // First, simply find the number of internal nodes.
        internSeqCounter = 0;
        for(nNum = 0; nNum < t->nNodes; nNum++) {
            if(!t->nodes[nNum]->isLeaf) {
                internSeqCounter++;
            }
        }
        //printf("There are %i internal nodes\n", internSeqCounter);

        t->internalSequences = malloc(d->nParts * sizeof(int **));
        if(!t->internalSequences) {
            printf("Problem allocating memory for internalSequences (1)\n");
            exit(1);
        }
        for(pNum = 0; pNum < d->nParts; pNum++) {
            dp = d->parts[pNum];
            t->internalSequences[pNum] = malloc(internSeqCounter * sizeof(int *));
            if(!t->internalSequences[pNum]) {
                printf("Problem allocating memory for internalSequences (2)\n");
                exit(1);
            }
            for(j = 0; j < internSeqCounter; j++) {
                t->internalSequences[pNum][j] = pivector(dp->nChar);
            }
        }
    }
	
    // Go thru the nodes in node order, and assign sequences to
    // simSequences.  The sequences are either
    // part->sequences[node->sequenceNum] for the terminal nodes, or
    // t->internalSequences[i][internSeqCounter] for the
    // internal nodes.

    for(pNum = 0; pNum < d->nParts; pNum++) {
        dp = d->parts[pNum];
        internSeqCounter = 0;
        for(nNum = 0; nNum < t->nNodes; nNum++) {
            if(t->nodes[nNum]->isLeaf) {
                t->simSequences[pNum][nNum] = dp->sequences[t->nodes[nNum]->seqNum];
            }
            else {
                t->simSequences[pNum][nNum] = t->internalSequences[pNum][internSeqCounter];
                internSeqCounter++;
            }
        }
    }

#if 0
    printf("\n--------------------------------\n");
    printf("\npart sequences:\n");
    for(i = 0; i < d->nParts; i++) {
        dp = d->parts[i];
        printf("    part %i\n", i);
        for(j = 0; j < dp->nTax; j++) {
            printf("        taxon %i: %li\n", j, (long int)dp->sequences[j]);
        }
    }
	
    internSeqCounter = 0;
    for(j = 0; j < t->nNodes; j++) {
        if(!t->nodes[j]->isLeaf) {
            internSeqCounter++;
        }
    }
    printf("internalSequences:\n");
    for(i = 0; i < d->nParts; i++) {
        printf("    part %i\n", i);
        for(j = 0; j < internSeqCounter; j++) {
            printf("        internalSeq %i: %li\n", j, (long int)t->internalSequences[i][j]);
        }
    }

    printf("simSequences:\n");
    for(i = 0; i < d->nParts; i++) {
        printf("    part %i\n", i);
        for(j = 0; j < t->nNodes; j++) {
            printf("        simSeq %i: %li\n", j, (long int)t->simSequences[i][j]);
        }
    }
#endif		


    // There's at least a couple of ways to do the "change in gdasrv
    // over the tree" thing.  The way I'll be doing it is to assign
    // each site to a category, and that site stays in that category.
    // However, the gamma shape might change over the tree, and so the
    // rate in a given category might change also.  But the cat that a
    // given site is in does not change.  So it makes sense to define
    // which category in data->parts[pNum]->simCats, which is
    // what I do.

    // An alternative would be to allow the category to which a site
    // belongs to change over the tree.  Not done.

    // Example mean rates for 4 cats for gamma shape = 0.5 and 1.5.
    //                 cat      0.5        1.5
    //                  0    0.033388   0.225323
    //                  1    0.251916   0.588556
    //                  2    0.820268   1.050422
    //                  3    2.894428   2.135699

    // Example mean rates for 2 cats
    //                 cat      0.5
    //                  0    0.142652
    //                  1    1.857348
    //                       --------
    //                  sum  2.0

    // malloc mp->charStatePicker
    for(pNum = 0; pNum < d->nParts; pNum++) {
        mp = t->model->parts[pNum];
        if(!mp->charStatePicker) {
            mp->charStatePicker = malloc(mp->dim * sizeof(double));
        }
        if(!mp->charStatePicker) {
            printf("Problem allocating memory for modelPart->charStatePicker.\n");
            exit(1);
        }
    }

    // allocate memory for data->part[pNum]->simCats
    for(pNum = 0; pNum < d->nParts; pNum++) {
        dp = d->parts[pNum];
        if(!dp->simCats) {
            dp->simCats = malloc(dp->nChar * sizeof(int));
            if(!dp->simCats) {
                printf("p4_tree simulate().  Failed to allocate simCats.\n");
                exit(1);
            }
        }
    }

    if(refTree) {
        if(t == refTree) {
            printf("Sim tree is the same as the refTree --- can't do that.\n");
            exit(1);
        }
        
        // printf("simulate() with refTree %li\n", (long int)refTree);
        for(pNum = 0; pNum < refTree->data->nParts; pNum++) {
            dp = refTree->data->parts[pNum];

            // printf("There are %i patterns\n", dp->nPatterns);
            // for(cNum = 0; cNum < dp->nChar; cNum++) {
            //     printf(" %i", dp->sequencePositionPatternIndex[cNum]);
            // }
            // printf("\n");

            for(seqPos = 0; seqPos < dp->nChar; seqPos++) {
                p4_drawAncState(refTree, pNum, seqPos);
                if(dp->drawAncStResults->isInvar) {
                    // printf("isInvar %i, invarChNum %i\n",
                    //        dp->drawAncStResults->isInvar,
                    //        dp->drawAncStResults->invarChNum);
                    t->simSequences[pNum][t->root->nodeNum][seqPos] = dp->drawAncStResults->invarChNum;
                    d->parts[pNum]->globalInvarSitesVec[seqPos] = 1;
                }
                else {
                    // printf("catNum %i, chStNum %i\n",
                    //        dp->drawAncStResults->catNum,
                    //        dp->drawAncStResults->chStNum);
                    d->parts[pNum]->simCats[seqPos] = dp->drawAncStResults->catNum;
                    t->simSequences[pNum][t->root->nodeNum][seqPos] = dp->drawAncStResults->chStNum;
                    d->parts[pNum]->globalInvarSitesVec[seqPos] = 0;
                }
            }

            // print the simCats, root sequence, and pInvar.
            // Note we have dp from refTree->data, and we have d->parts[pNum] where the sim results are.
            if((0)) {
                if(d->parts[pNum]->simCats) {
                    //printf("part %i dp->simCats.  nCat=%i\n", pNum, refTree->model->parts[pNum]->nCat);
                    for(i = 0; i < dp->nChar; i++) {
                        printf("%1i", d->parts[pNum]->simCats[i]);
                    }
                    printf("\n");
                }
                //printf("Part %i root node sequence \n", pNum);
                for(i = 0; i < dp->nChar; i++) {
                    printf("%1i", t->simSequences[pNum][t->root->nodeNum][i]);
                }
                printf("\n");
                if(d->parts[pNum]->globalInvarSitesVec) {
                    //printf("globalInvarSitesVec\n");
                    for(i = 0; i < dp->nChar; i++) {
                        printf("%1i", d->parts[pNum]->globalInvarSitesVec[i]);
                    }
                    printf("\n");
                }
            }

        }

    }   // if(refTree)

    else {
        // Fill data->parts[pNum]->simCats, a vector of ints, nChar long.
        // We are going to assume here that the frequencies of each cat
        // are all equal.  If that is not the case, I would need to refer
        // to gdasrv->freqs.
        for(pNum = 0; pNum < d->nParts; pNum++) {
            dp = d->parts[pNum];
            mp = t->model->parts[pNum];
            if(mp->nCat > 1) {
                for(i = 0; i < dp->nChar; i++) {
                    dp->simCats[i] = (int)floor(((double)mp->nCat) * gsl_rng_uniform(g));
                }
            } else {
                for(i = 0; i < dp->nChar; i++) {
                    dp->simCats[i] = 0;
                }
            }
            
            if((0)) {
                printf("part %i dp->simCats.  nCat=%i\n", pNum, mp->nCat);
                for(i = 0; i < dp->nChar; i++) {
                    printf("%1i", dp->simCats[i]);
                }
                printf("\n");
            }
        }
	

        // fill picker mp->charStatePicker based on the composition of part[pNum], then fill the root node simSequences
        for(pNum = 0; pNum < d->nParts; pNum++) {
            dp = d->parts[pNum];
            mp = t->model->parts[pNum];
            mp->charStatePicker[0] = mp->comps[t->root->compNums[pNum]]->val[0];
            for(j = 1; j < mp->dim - 1; j++) {
                mp->charStatePicker[j] = mp->charStatePicker[j - 1] + mp->comps[t->root->compNums[pNum]]->val[j];
            }
            mp->charStatePicker[mp->dim -1] = 1.0;

            if((0)) {
                printf("Part %i mp->charStatePicker is\n", pNum);
                for(j = 0; j < mp->dim; j++) {
                    printf("%7.4f", mp->charStatePicker[j]);
                }
                printf("\n");
            }

            // simulate root sequence
            for(j = 0; j < d->parts[pNum]->nChar; j++) {
                theRanDoubleUpToOne = gsl_rng_uniform(g);
                for(k = 0; k < mp->dim; k++) {
                    if(theRanDoubleUpToOne < mp->charStatePicker[k]) {
                        t->simSequences[pNum][t->root->nodeNum][j] = k;
                        break;
                    }
                }
            }
        }


        // printf("Root node is nodeNum %i\n", t->root->nodeNum);
        // for(pNum = 0; pNum < d->nParts; pNum++) {
        //     printf("Part %i root node sequence is (%li)\n", 
        //            pNum, (long int)t->simSequences[pNum][t->root->nodeNum]);
        //     for(j = 0; j < d->parts[pNum]->nChar; j++) {
        //         printf("%1i", t->simSequences[pNum][t->root->nodeNum][j]);
        //     }
        //     printf("\n");
        // }
        // //exit(1);

        /*
          Fill the invarSites vectors.  1 for invar, 0 for potentially
          variable.

          Each data part has a globalInvarSitesVec.

        */
        for(pNum = 0; pNum < d->nParts; pNum++) {
            dp = d->parts[pNum];
            mp = t->model->parts[pNum];
            if(mp->pInvar->val[0] > 0.0) {
                for(i = 0; i < dp->nChar; i++) {
                    theRanDoubleUpToOne = gsl_rng_uniform(g);
                    if(theRanDoubleUpToOne < mp->pInvar->val[0]) {
                        dp->globalInvarSitesVec[i] = 1;
                    }
                    else {
                        dp->globalInvarSitesVec[i] = 0;
                    }
                }
            } else {
                for(i = 0; i < dp->nChar; i++) {
                    dp->globalInvarSitesVec[i] = 0;
                }
            }
            // for(i = 0; i < dp->nChar; i++) {
            // 	printf("%1i", dp->globalInvarSitesVec[i]);
            // }
            // printf("\n");
        }
    }
		
    /*
      What the root node sequence (which already exists, from above)
      mutates into along the tree is given by the bigP probability
      matrices.  There is a different bigP for each gamma distributed
      rate cat, and so we have a bigPDeck.  And there is a bigPDeck
      for each data partition, so we have bigPDecks.
	  
      In order to actually use it, it is converted to a pickerDecks,
      where the numbers are the running sums of each row of the
      bigPDecks.  That way I can use a random number from 0 to less than
      1 to choose what the sequences mutate to.

      See p4_node.c for these functions.
    */
	

    for(i = 0; i < t->nNodes; i++) {
        n = t->nodes[i];
        if(n != t->root) {
            p4_calculateBigPDecks(n);
            p4_calculatePickerDecks(n); // allocates if needed
        }
    }

#if 0
    n = t->nodes[1];
    //void dump_psdmatrix(double **m, int dim); // prints a square double matrix
    dump_psdmatrix(n->bigPDecks[0][0], 4);
#endif

    /*
      And now the moment we have all been waiting for.  We mutate the
      root sequence along the tree.  The nodes and simSequences are in
      t->nodes, but we traverse the tree in preOrder.

      If its an invariant site, leave it as the parent state.
      If its a variable site, mutate it using the pickerDecks.

    */
	
    //printf("About to simulate\n"); fflush(stdout);


    //printf("t->nNodes = %i\n", t->nNodes);
    //for(nNum = 0; nNum < t->nNodes; nNum++) {
    //	printf("preOrder[%2i] = %i\n", nNum, t->preOrder[nNum]);
    //}

    for(pNum = 0; pNum < d->nParts; pNum++) {
        mp = t->model->parts[pNum];
        for(nNum = 0; nNum < t->nNodes; nNum++) {
            tOrd = t->preOrder[nNum];
            //printf("got tOrd=%i\n", tOrd);
            if(tOrd != NO_ORDER) {
                n = t->nodes[tOrd];
                if(n != t->root) {
                    //printf("\nnNum = %i, Node %i, Taking parentState from t->simSequences %li\n", 
                    //	   nNum, n->nodeNum, (long int)t->simSequences[pNum][n->parent->nodeNum]);
                    for(k = 0; k < d->parts[pNum]->nChar; k++) {
                        parentState = t->simSequences[pNum][n->parent->nodeNum][k];
                        if(d->parts[pNum]->globalInvarSitesVec[k]) { // if its an invariable site
                            t->simSequences[pNum][tOrd][k] = parentState;
                        }
                        else {
                            simCat = d->parts[pNum]->simCats[k]; // if nCat is 1, this will be zero.
                            theRanDoubleUpToOne = gsl_rng_uniform(g);
                            for(l = 0; l < mp->dim; l++) {
                                if(theRanDoubleUpToOne < n->pickerDecks[pNum][simCat][parentState][l]) {
                                    t->simSequences[pNum][tOrd][k] = l;
                                    break;
                                }
                            }
                            //if(parentState != l) {
                            //	printf("[%i, %i], ", parentState, l);
                            //}
                        }
                    }
                }
            }
        }
    }
    //printf("...finished simulate\n"); fflush(stdout);

#if 0
    printf("\n");
    for(i = 0; i < d->nParts; i++) {
        printf("Part %i: simSequences for each node\n", i);
        for(j = 0; j < t->nNodes; j++) {
            tOrd = t->preOrder[j];
            printf("Node %i, tOrd %i\n", j, tOrd);
            for(k = 0; k < d->parts[i]->nChar; k++) {
                printf("%1i", t->simSequences[i][j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
    printf("\n");
    for(i = 0; i < d->nParts; i++) {
        printf("Part %i: sequences for each taxon\n", i);
        for(j = 0; j < d->nTax; j++) {
            printf("Taxon %i\n", j);
            for(k = 0; k < d->parts[i]->nChar; k++) {
                printf("%1i", d->parts[i]->sequences[j][k]);
            }
            printf("\n");
        }
        printf("\n");
    }
    printf("\n");
#endif

#if 0
    // Some diagnostics code.
    {
        int sumInvar = 0;
        int hits = 0;
        int *catSums;
        int *charSums;
        //int hits2 = 0;
        
        for(pNum = 0; pNum < d->nParts; pNum++) {
            dp = d->parts[pNum];
            mp = t->model->parts[pNum];
            catSums = pivector(mp->nCat);
            for(i = 0; i < mp->nCat; i++) {
                catSums[i] = 0;
            }
            charSums = pivector(dp->dim);
            for(i = 0; i < dp->dim; i++) {
                charSums[i] = 0;
            }
            
            for(seqPos = 0; seqPos < dp->nChar; seqPos++) {
                // Look for a certain pattern; only for three=node trees.
                if(t->simSequences[pNum][1][seqPos] == 0 && t->simSequences[pNum][2][seqPos] == 1) {
                    hits += 1;
                    //if(dp->simCats[seqPos] == 0 && t->simSequences[pNum][t->root->nodeNum][seqPos] == 0) {
                    //    hits2 += 1;
                    //}
                    if(dp->globalInvarSitesVec[seqPos]) {
                        sumInvar += 1;
                    } else {
                        catSums[dp->simCats[seqPos]] += 1;
                    }
                    charSums[t->simSequences[pNum][t->root->nodeNum][seqPos]] += 1;
                }
            }
            
            printf("hits %i, pInvar %i\n", hits, sumInvar);
            //printf("hits2 %i\n", hits2);
            printf("pInvar / hits : %.4f \n", (double)sumInvar / (double)hits);
            printf("catSums ");
            for(i = 0; i < mp->nCat; i++) {
                printf("%6i ", catSums[i]);
            }
            printf("\n");

            printf("catSums / hits : ");
            for(i = 0; i < mp->nCat; i++) {
                printf("%.4f ", (double)catSums[i] / (double)hits);
            }
            printf("\n");

            printf("charSums: ");
            for(i = 0; i < dp->dim; i++) {
                printf("%6i ", charSums[i]);
            }
            printf("\n");
            printf("charSums / hits : ");
            for(i = 0; i < dp->dim; i++) {
                printf("%.4f ", (double)charSums[i] / (double)hits);
            }
            printf("\n");
            
            free(catSums);
            catSums = NULL;
            free(charSums);
            charSums = NULL;


        }
    }
#endif

#if 0
    printf("\n");
    printf("part[0] is %li\n", (long int)d->parts[0]);
    printf("Part 0: sequence for first taxon\n");
    for(k = 0; k < d->parts[0]->nChar; k++) {
        printf("%1i", d->parts[0]->sequences[0][k]);
    }
    printf("\n");
#endif

    //printf("p4_treeSim.c: finished simulate\n");
}

//==================================================
//==================================================


void p4_drawAncState(p4_tree *t, int partNum, int seqPos)
{
    part  *dp;
    p4_modelPart  *mp;
    int catNum, chStNum, patNum;
    int gotIt, isInvar;
    //double diff;
    double sLike, sLikeC;
    double mySum, mySum2, myPInvarSum;
    int i,k;
    //int j;
    double theRanDouble;


    dp = t->data->parts[partNum];
    mp = t->model->parts[partNum];

    patNum = dp->sequencePositionPatternIndex[seqPos];

    // printf("There are %i patterns\n", dp->nPatterns);
    // int cNum;
    // for(cNum = 0; cNum < dp->nChar; cNum++) {
    //     printf(" %i", dp->sequencePositionPatternIndex[cNum]);
    // }
    // printf("\n");
            
    // Make the ancStPicker.  It should be at least dim * nCat, and
    // optionally more if pInvar.  Usually only one more space will be
    // needed if pInvar, but if ambigs are involved then more will be
    // needed.  The maximum that pInvar needs is dim --- , for a site
    // composed of constant Ns or Xs, then dim extra will be needed.
    if(!mp->ancStPicker) {
        i = mp->dim * (mp->nCat + 1);
        mp->ancStPicker = malloc(i * sizeof(double));
        if(!mp->ancStPicker) {
            printf("Problem allocating memory for modelPart->ancStPicker.\n");
            exit(1);
        }
    }
    
    for(i = 0; i <  (mp->dim * (mp->nCat + 1)); i++) {
        mp->ancStPicker[i] = 0.0;
    }
    
    if(!dp->drawAncStResults) {
        dp->drawAncStResults = malloc(sizeof(p4_drawAncStResults));
        if(!dp->drawAncStResults) {
            printf("Problem allocating memory for dataPart->drawAncStResults.\n");
            exit(1);
        }
    }
    
    // printf("pat %3i\n ", patNum);
    // for(catNum = 0; catNum < mp->nCat; catNum++) {
    //     printf("cat %i:", catNum);
    //     for(chStNum = 0; chStNum < dp->dim; chStNum++) {
    //         printf(" %8.5f", t->root->cl[partNum][catNum][chStNum][patNum]);
    //     }
    //     printf("\n");
    // }
            
    // printf("mp->freqsTimesOneMinusPInvar[catNum]: ");
    // for(catNum = 0; catNum < mp->nCat; catNum++) {
    //     printf("[catNum %i] %8.6f ", catNum, mp->freqsTimesOneMinusPInvar[catNum]);
    // }
    // printf("\n");

    // printf("dp->globalInvarSitesArray[chStNum][patNum]\n");
    // for(chStNum = 0; chStNum < dp->dim; chStNum++) {
    //     for(k = 0; k < dp->nPatterns; k++) {
    //         printf("%i", dp->globalInvarSitesArray[chStNum][k]);
    //     }
    //     printf("\n");
    // }

    // double sLike, sLikeC
    sLike = 0.0;
    i = 0;  // counter for mp-ancStPicker
    for(catNum = 0; catNum < mp->nCat; catNum++) {
        for(chStNum = 0; chStNum < dp->dim; chStNum++) {
            sLikeC = mp->comps[t->root->compNums[partNum]]->val[chStNum] * 
                t->root->cl[partNum][catNum][chStNum][patNum];
            sLikeC *= mp->freqsTimesOneMinusPInvar[catNum];
            sLike += sLikeC;
            mp->ancStPicker[i] = sLike;
            i++;
            // printf("    patNum %i  catNum %i chStNum %i  sLikeC %10.8f  sLike %10.8f\n",
            //         patNum, catNum, chStNum, sLikeC, sLike);
        }
    }
    if(mp->pInvar->val[0]) { // do pInvar
        for(chStNum = 0; chStNum < dp->dim; chStNum++) {
            if(dp->globalInvarSitesArray[chStNum][patNum]) {
                sLikeC = mp->comps[t->root->compNums[partNum]]->val[chStNum] * mp->pInvar->val[0];
                sLike += sLikeC;
                mp->ancStPicker[i] = sLike;
            }
            i++;
            //printf("    patNum %i  pInvar   chStNum %i  sLikeC %10.8f  sLike %10.8f\n",
            //       patNum, chStNum, sLikeC, sLike);
            
        }
        //printf("\n");
    }

    if((0)) {
        // Some diagnostics ...
        for(i = 0; i <  (mp->dim * (mp->nCat + 1)); i++) {
            printf("    %3i  %10.8f \n", i, mp->ancStPicker[i]);
        }

        // Calculate the prob of each char state, summed over rateCats.
        i = 0;
        mySum2 = 0.0;
        myPInvarSum = 0.0;
        for(chStNum = 0; chStNum < dp->dim; chStNum++) {
            mySum = 0.0;
            for(catNum = 0; catNum < mp->nCat; catNum++) {
                sLikeC = mp->comps[t->root->compNums[partNum]]->val[chStNum] * 
                    t->root->cl[partNum][catNum][chStNum][patNum];
                sLikeC *= mp->freqsTimesOneMinusPInvar[catNum];
                //printf("chStNum %2i, catNum %i, sLikeC = %10.8f \n", chStNum, catNum, sLikeC);
                mySum += sLikeC;
            }

            // if it is invar
            // Using the globalInvarSitesArray accommodates ambiguities.
            if(dp->globalInvarSitesArray[chStNum][patNum]) {
                sLikeC = mp->comps[t->root->compNums[partNum]]->val[chStNum] * mp->pInvar->val[0];
                //printf("chStNum %2i, pInvar, sLikeC = %10.8f \n", chStNum, sLikeC);
                mySum += sLikeC;
                myPInvarSum += sLikeC;
            }

            mySum /= sLike;
            mySum2 += mySum;
            printf("Char state %2i  prob %10.8f \n", chStNum, mySum);
        }
        myPInvarSum /= sLike;
        printf("pInvar prob %10.8f \n", myPInvarSum);
        printf("Cumulative prob of char states %10.8f \n", mySum2);

        // Calculate the prob of each rate category
        mySum2 = myPInvarSum;   // cumulative over site rates
        for(catNum = 0; catNum < mp->nCat; catNum++) {
            mySum = 0.0;
            for(chStNum = 0; chStNum < dp->dim; chStNum++) {
                sLikeC = mp->comps[t->root->compNums[partNum]]->val[chStNum] * 
                    t->root->cl[partNum][catNum][chStNum][patNum];
                sLikeC *= mp->freqsTimesOneMinusPInvar[catNum];
                //printf("chStNum %2i, catNum %i, sLikeC = %10.8f \n", chStNum, catNum, sLikeC);
                mySum += sLikeC;
            }
            mySum /= sLike;
            printf("catNum %i prob %10.8f \n", catNum, mySum);
            mySum2 += mySum;
        }
        printf("Cumulative prob of rate cats, including pInvar:  %10.8f \n", mySum2);
    }


    // check that dp->siteLikes has the same numbers as calculated
    // here.  This would only work if siteLikes have been calculated
    // in the t, and we cannot be sure of that, so this bit is
    // only for debugging. 

    // One little problem is that dp->siteLikes are for sites, not
    // for patterns.  So I need to find the first i where
    // dp->sequencePositionPatternIndex[i] is equal to patNum.

    // for(j = 0; j < dp->nChar; j++) {
    //     if(dp->sequencePositionPatternIndex[j] == patNum) {
    //         break;
    //     }
    // }
    // diff = dp->siteLikes[j] - sLike;
    // if(fabs(diff) > 1.0e-12) {
    //     printf("p4_simulate() patNum %i, site like here (position %i %g) and t site like (%g) differ by %10.8f, %g\n", 
    //            patNum, j, sLike, dp->siteLikes[patNum], diff, diff);
    //     printf("p4_simulate() patNum %i, site like here and t site like differ by %10.8f, %g\n", 
    //            patNum, diff, diff);
    //     exit(1);
    // }

    // The random number is not from zero to one, it is from zero to sLike.
    // That way we do not need to change the ancStPicker
    //theRanDouble = ((double)random()) / ((double)((long)RAND_MAX));
    theRanDouble = (double)random() / ((double)(RAND_MAX) + 1.0) ;
    theRanDouble *= sLike;
    // printf("Got theRanDouble %.10f\n", theRanDouble);
                
    // Now pick an ancestral state --- cat, chSt, invar, and invarChar.
    gotIt = 0;
    k = 0;  // another counter for mp-ancStPicker
    for(catNum = 0; catNum < mp->nCat; catNum++) {
        for(chStNum = 0; chStNum < dp->dim; chStNum++) {
            if(theRanDouble < mp->ancStPicker[k]) {
                gotIt = 1;
                //printf("xGot anc state: patNum %i, catNum %i, chStNum %i\n", patNum, catNum, chStNum);
            }
            k++;
            if(gotIt) break;
        }
        if(gotIt) break;
    }
    isInvar = 0;
    if(!gotIt) {
        if(mp->pInvar->free) { // do pInvar
            for(chStNum = 0; chStNum < dp->dim; chStNum++) {
                if(dp->globalInvarSitesArray[chStNum][patNum]) {
                    if(theRanDouble < mp->ancStPicker[k]) {
                        isInvar = 1;
                        //printf("yGot anc state: patNum %i, pInvar, chStNum %i\n", patNum, chStNum);
                        gotIt = 1;
                    }
                }
                k++;
                if(gotIt) break;
            }
        }
    }
    
    if(!gotIt) {
        printf("Something is wrong with the ancestral state picker. gotIt is zero.\n");
        exit(1);
    }
    // if(!isInvar) {
    //     printf(" Got anc state: patNum %i, catNum %i, chStNum %i\n", patNum, catNum, chStNum);
    // }
    // else {
    //     printf(" Got anc state: patNum %i, pInvar chStNum %i\n", patNum, chStNum);
    // }

    if(!isInvar) {
        dp->drawAncStResults->catNum = catNum;
        dp->drawAncStResults->chStNum = chStNum;
        dp->drawAncStResults->isInvar = 0;
        dp->drawAncStResults->invarChNum = -1;
    }
    else {
        dp->drawAncStResults->isInvar = isInvar;
        dp->drawAncStResults->invarChNum = chStNum;
        dp->drawAncStResults->catNum = -1;
        dp->drawAncStResults->chStNum = -1;
    }
    if((0)) {
        printf("chStNum %i, catNum %i, isInvar %i, invarChNum %i\n", 
               dp->drawAncStResults->chStNum,
               dp->drawAncStResults->catNum,
               dp->drawAncStResults->isInvar,
               dp->drawAncStResults->invarChNum);
    }

}

// This is called by a pf method, in pfmodule.c
void p4_drawAncStateP(p4_tree *t, int partNum, int seqPos, int *draw)
{
    part *dp;

    dp = t->data->parts[partNum];
    p4_drawAncState(t, partNum, seqPos);
    draw[0] = dp->drawAncStResults->chStNum;
    draw[1] = dp->drawAncStResults->catNum;
    draw[2] = dp->drawAncStResults->isInvar;
    draw[3] = dp->drawAncStResults->invarChNum;
}

PyObject *p4_expectedCompositionCounts(p4_tree *t, int partNum)
{
    int i, j, k;
    int nNum;
    //PyObject	*bigTup;
    PyObject	*mediumTup;
    PyObject	*smallTup;
    int terminalNodeCount = 0;
    double    **expPi;
    int nGaps;
    p4_modelPart  *mp = t->model->parts[partNum];
    part          *dp = t->data->parts[partNum];
	
    // We want the expected counts of the terminal nodes only.  They
    // should be in the order of the sequences.  Ie the first set
    // should be from node->seqNum=0, the second set should be
    // from node->seqNum=1, and so on.

    //printf("p4_treeSim.c: p4_expectedCompositionCounts here\n");

    if((0)) {
        printf("    %15s ", "preOrder");
        for(i = 0; i < t->nNodes; i++) {
            printf("%i ", t->preOrder[i]);
        }
        printf("\n");
    }

    // Now tell each node, in tree order, to p4_calculateExpectedCharFreq.
    // This is inefficient, because each node does the calculation 
    // for each part.  But this function (ie p4_expectedCharFreq()) asks 
    // for the expected composition for a single part.

    for(nNum = 0; nNum < t->nNodes; nNum++) {
        if(t->preOrder[nNum] != NO_ORDER) {
            p4_calculateExpectedComp(t->nodes[t->preOrder[nNum]]);
        }
    }
	
    expPi = pdmatrix(t->nNodes, mp->dim);
    for(nNum = 0; nNum < t->nNodes; nNum++) {
        for(j = 0; j < mp->dim; j++) {
            expPi[nNum][j] = 0.0;
        }
    }

    // Add up the charFreqs in the rate categories ...
    for(nNum = 0; nNum < t->nNodes; nNum++) {
        for(j = 0; j < mp->nCat; j++) {
            for(k = 0; k < mp->dim; k++) {
                expPi[nNum][k] += t->nodes[nNum]->expectedComp[partNum][j][k];
            }
        }
    }

    // ... and divide by the number of rate categories to get mean.
    for(nNum = 0; nNum < t->nNodes; nNum++) {
        for(k = 0; k < mp->dim; k++) {
            expPi[nNum][k] /= (double)mp->nCat;
        }
    }

    for(nNum = 0; nNum < t->nNodes; nNum++) {
        if(t->nodes[nNum]->isLeaf) {
            terminalNodeCount++;
        }
    }
    mediumTup = PyTuple_New(terminalNodeCount);
    for(nNum = 0; nNum < t->nNodes; nNum++) {
        if(t->nodes[nNum]->isLeaf) {
            smallTup = PyTuple_New(mp->dim);
			
            nGaps = 0;
            for(k = 0; k < dp->nChar; k++) {
                if(dp->sequences[t->nodes[nNum]->seqNum][k] == GAP_CODE ||
                   dp->sequences[t->nodes[nNum]->seqNum][k] == QMARK_CODE)
                    nGaps++;
            }
				   
            for(j = 0; j < mp->dim; j++) {
                PyTuple_SetItem(smallTup, j, 
                                PyFloat_FromDouble((double)(dp->nChar - nGaps) * expPi[nNum][j]));
            }
            PyTuple_SetItem(mediumTup, t->nodes[nNum]->seqNum, smallTup);
        }
    }
    free_pdmatrix(expPi);
    expPi = NULL;
    return mediumTup;
}


PyObject *p4_expectedComposition(p4_tree *t)
{
    int i, j, k, partNum;
    int nNum;
    PyObject	*bigTup;
    PyObject	*mediumTup;
    PyObject	*smallTup;
    int terminalNodeCount = 0;
    double    **expPi = NULL;
    //int nGaps;
    p4_modelPart  *mp = NULL;
    //part          *dp = NULL;
	
    // We want the expected comp of the terminal nodes only.  They
    // should be in the order of the sequences.  Ie the first set
    // should be from node->seqNum=0, the second set should be
    // from node->seqNum=1, and so on.

    //printf("p4_treeSim.c: p4_expectedCompositionC here\n");

    if((0)) {
        printf("    %15s ", "preOrder");
        for(i = 0; i < t->nNodes; i++) {
            printf("%i ", t->preOrder[i]);
        }
        printf("\n");
    }

    for(nNum = 0; nNum < t->nNodes; nNum++) {
        if(t->nodes[nNum]->isLeaf) {
            terminalNodeCount++;
        }
    }
 
	
    // Now tell each node, in preorder, to p4_calculateExpectedComp.

    for(nNum = 0; nNum < t->nNodes; nNum++) {
        if(t->preOrder[nNum] != NO_ORDER) {
            //printf("about to calculate expected comp of node %i\n", t->preOrder[nNum]);
            p4_calculateExpectedComp(t->nodes[t->preOrder[nNum]]);
        }
    }

    bigTup = PyTuple_New(t->model->nParts);
    for(partNum = 0; partNum < t->model->nParts; partNum++) {
        mp = t->model->parts[partNum];
        //dp = t->data->parts[partNum];
        expPi = pdmatrix(t->nNodes, mp->dim);
        for(nNum = 0; nNum < t->nNodes; nNum++) {
            for(j = 0; j < mp->dim; j++) {
                expPi[nNum][j] = 0.0;
            }
        }

        // Add up the charFreqs in the rate categories ...
        for(i = 0; i < t->nNodes; i++) {
            if(t->preOrder[i] != NO_ORDER) {
                nNum = t->preOrder[i];
                for(j = 0; j < mp->nCat; j++) {
                    for(k = 0; k < mp->dim; k++) {
                        expPi[nNum][k] += t->nodes[nNum]->expectedComp[partNum][j][k];
                    }
                }
            }
        }

        // ... and divide by the number of rate categories to get mean.
        for(i = 0; i < t->nNodes; i++) {
            if(t->preOrder[i] != NO_ORDER) {
                nNum = t->preOrder[i];
                for(k = 0; k < mp->dim; k++) {
                    expPi[nNum][k] /= (double)mp->nCat;
                }
            }
        }

        mediumTup = PyTuple_New(terminalNodeCount);
        for(nNum = 0; nNum < t->nNodes; nNum++) {
            if(t->nodes[nNum]->isLeaf) {

                smallTup = PyTuple_New(mp->dim);
                for(j = 0; j < mp->dim; j++) {
                    PyTuple_SetItem(smallTup, j, PyFloat_FromDouble(expPi[nNum][j]));
                }
                PyTuple_SetItem(mediumTup, t->nodes[nNum]->seqNum, smallTup);
            }
        }
        free_pdmatrix(expPi);
        expPi = NULL;
        PyTuple_SetItem(bigTup, partNum, mediumTup);
    }
	
    return bigTup;
}


