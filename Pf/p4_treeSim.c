#include "pftypes.h"
#include "p4_tree.h"
#include "pmatrices.h"
#include "defines.h"
#include "util.h"
#include "p4_node.h"
#include <stdio.h>
#include <stdlib.h>





void p4_simulate(p4_tree *t)
{

	double  *picker = NULL; 
	double  theRanDoubleUpToOne;
	int i, j, k, l;
	int pNum, nNum, cNum, rNum;
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
			exit(0);
		}
		for(pNum = 0; pNum < d->nParts; pNum++) {
			t->simSequences[pNum] = malloc(t->nNodes * sizeof(int *));
			if(!t->simSequences[pNum]) {
				printf("Problem allocating memory for simSequences (2)\n");
				exit(0);
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
			exit(0);
		}
		for(pNum = 0; pNum < d->nParts; pNum++) {
			dp = d->parts[pNum];
			t->internalSequences[pNum] = malloc(internSeqCounter * sizeof(int *));
			if(!t->internalSequences[pNum]) {
				printf("Problem allocating memory for internalSequences (2)\n");
				exit(0);
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

	// make a picker, that is big enough for the biggest dim or nCat
	i = 0;  // the biggest dim or nCat so far
	for(pNum = 0; pNum < d->nParts; pNum++) {
		mp = t->model->parts[pNum];
		if(mp->dim > i) {  // if covarion, model dim will be bigger than data dim
			i = mp->dim;
		}		
		if(mp->nCat > i) {
			i = mp->nCat;
		}
	}
	//printf("malloc'ing a picker %i long.\n", i);
	picker = malloc(i * sizeof(double));   // it is freed at the end of this function
	if(!picker) {
		printf("Problem allocating memory for picker.\n");
		exit(0);
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


	// Fill data->parts[pNum]->simCats, a vector of ints, nChar long.
	// We are going to assume here that the frequencies of each cat
	// are all equal.  If that is not the case, I would need to refer
	// to gdasrv->freqs.
	for(pNum = 0; pNum < d->nParts; pNum++) {
		dp = d->parts[pNum];
		mp = t->model->parts[pNum];
		if(mp->isMixture) {
			// fill a picker based on mp->mixture->freqs
			picker[0] = mp->mixture->freqs[0];
			for(i = 1; i < mp->nCat - 1; i++) {
				picker[i] = picker[i - 1] + mp->mixture->freqs[i];
			}
			picker[mp->nCat - 1] = 1.0;
			if(0) {
				printf("Part %i: picker for simCats is: ", pNum);
				for(j = 0; j < mp->nCat; j++) {
					printf(" %7.4f", picker[j]);
				}
				printf("\n");
			}
			for(j = 0; j < dp->nChar; j++) {
				theRanDoubleUpToOne = ranDoubleUpToOne();
				for(k = 0; k < mp->nCat; k++) {
					if(theRanDoubleUpToOne < picker[k]) {
						dp->simCats[j] = k;
						break;
					}
				}
			}
			
		}
		else {
			if(mp->nCat > 1) {
				for(i = 0; i < dp->nChar; i++) {
					dp->simCats[i] = (int)floor(((double)mp->nCat) * ranDoubleUpToOne());
				}
			} else {
				for(i = 0; i < dp->nChar; i++) {
					dp->simCats[i] = 0;
				}
			}
		}
		if(0) {
			printf("part %i dp->simCats.  nCat=%i\n", pNum, mp->nCat);
			for(i = 0; i < dp->nChar; i++) {
				printf("%1i", dp->simCats[i]);
			}
			printf("\n");
		}
	}

#if 0
	mp = t->model->parts[0];
	dp = d->parts[0];
	for(j = 0; j < mp->nCat; j++) {
		picker[j] = 0.0;
	}
	for(cNum = 0; cNum < dp->nChar; cNum++) {
		//printf("%1i", dp->simCats[cNum]);
		picker[dp->simCats[cNum]] += 1.0;
	}
	for(j = 0; j < mp->nCat; j++) {
		picker[j] /= (double)d->parts[0]->nChar;
	}
	printf("Part 0 simCats comp  is:");
	for(j = 0; j < mp->nCat; j++) {
		printf(" %6.3f", picker[j]);
	}
	printf("\n");
	//exit(1);
#endif
	

	// fill picker based on the composition of part[pNum], then fill the root node simSequences
	for(pNum = 0; pNum < d->nParts; pNum++) {
		dp = d->parts[pNum];
		mp = t->model->parts[pNum];
		if(mp->isMixture) {
			//printf("simulate()  isMixture.\n");
			i = 0; // counter for "nCat"
			for(cNum = 0; cNum < mp->nComps; cNum++) {
				// fill a comp picker for cats with this comp
				picker[0] = mp->comps[cNum]->val[0];
				for(j = 1; j < mp->dim - 1; j++) {
					picker[j] = picker[j - 1] + mp->comps[cNum]->val[j];
				}
				picker[mp->dim -1] = 1.0;
				if(0) {
					printf("Part %i: picker for cat %i is: ", pNum, i);
					for(j = 0; j < mp->dim; j++) {
						printf(" %7.4f", picker[j]);
					}
					printf("\n");
				}
				for(rNum = 0; rNum < mp->nRMatrices; rNum++) {  // use the same picker for all cats with the same comp
					// simulate root sequence only for sites where dp->simCats[j] == i
					for(j = 0; j < d->parts[pNum]->nChar; j++) {
						if(dp->simCats[j] == i) {
							theRanDoubleUpToOne = ranDoubleUpToOne();
							for(k = 0; k < mp->dim; k++) {
								if(theRanDoubleUpToOne < picker[k]) {
									t->simSequences[pNum][t->root->nodeNum][j] = k;
									break;
								}
							}
						}
					}
					i++;
				}
			}
		}
		else {
			picker[0] = mp->comps[t->root->compNums[pNum]]->val[0];
			for(j = 1; j < mp->dim - 1; j++) {
				picker[j] = picker[j - 1] + mp->comps[t->root->compNums[pNum]]->val[j];
			}
			picker[mp->dim -1] = 1.0;

			if(0) {
				printf("Part %i picker is\n", pNum);
				for(j = 0; j < mp->dim; j++) {
					printf("%7.4f", picker[j]);
				}
				printf("\n");
			}

			// simulate root sequence
			for(j = 0; j < d->parts[pNum]->nChar; j++) {
				theRanDoubleUpToOne = ranDoubleUpToOne();
				for(k = 0; k < mp->dim; k++) {
					if(theRanDoubleUpToOne < picker[k]) {
						t->simSequences[pNum][t->root->nodeNum][j] = k;
						break;
					}
				}
			}
		}
	}


#if 0
	printf("Root node is nodeNum %i\n", t->root->nodeNum);
	for(pNum = 0; pNum < d->nParts; pNum++) {
		printf("Part %i root node sequence is (%li)\n", 
			   pNum, (long int)t->simSequences[pNum][t->root->nodeNum]);
		for(j = 0; j < d->parts[pNum]->nChar; j++) {
			printf("%1i", t->simSequences[pNum][t->root->nodeNum][j]);
		}
		printf("\n");
	}
	//exit(1);
#endif

#if 0
	mp = t->model->parts[0];
	for(j = 0; j < mp->dim; j++) {
		picker[j] = 0.0;
	}
	for(cNum = 0; cNum < d->parts[0]->nChar; cNum++) {
		//printf("%1i", t->simSequences[0][t->root->nodeNum][cNum]);
		picker[t->simSequences[0][t->root->nodeNum][cNum]] += 1.0;
	}
	for(j = 0; j < mp->dim; j++) {
		picker[j] /= (double)d->parts[0]->nChar;
	}
	printf("Part 0 root comp is:");
	for(j = 0; j < mp->dim; j++) {
		printf(" %6.3f", picker[j]);
	}
	printf("\n");
	exit(1);
#endif
	

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
				theRanDoubleUpToOne = ranDoubleUpToOne();
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
		//for(i = 0; i < dp->nChar; i++) {
		//	printf("%1i", dp->globalInvarSitesVec[i]);
		//}
		//printf("\n");
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
							theRanDoubleUpToOne = ranDoubleUpToOne();
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

	// Collapse the covarion dim symbols into dim/2
	for(pNum = 0; pNum < d->nParts; pNum++) {
		mp = t->model->parts[pNum];
		dp = d->parts[pNum];
		if(mp->doTSCovarion) {
			for(nNum = 0; nNum < t->nNodes; nNum++) {
				for(j = 0; j < dp->nChar; j++) {
					if(t->simSequences[pNum][nNum][j] >= dp->dim) {
						t->simSequences[pNum][nNum][j] -= dp->dim;
					}
				}
			}
		}
	}


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
	printf("\n");
	printf("part[0] is %li\n", (long int)d->parts[0]);
	printf("Part 0: sequence for first taxon\n");
	for(k = 0; k < d->parts[0]->nChar; k++) {
		printf("%1i", d->parts[0]->sequences[0][k]);
	}
	printf("\n");
#endif

	free(picker);
	//printf("p4_treeSim.c: finished simulate\n");
}

//==================================================
//==================================================


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

	if(0) {
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
	part          *dp = NULL;
	
	// We want the expected comp of the terminal nodes only.  They
	// should be in the order of the sequences.  Ie the first set
	// should be from node->seqNum=0, the second set should be
	// from node->seqNum=1, and so on.

	//printf("p4_treeSim.c: p4_expectedCompositionC here\n");

	if(0) {
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
		dp = t->data->parts[partNum];
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


