#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Python.h>
#include "pmatrices.h"
#include "defines.h"
#include "pftypes.h"
#include "part.h"
#include "eig.h"
#include "util.h"
#include <gsl/gsl_linalg.h>

#define DTOLERANCE	1e-37  // uncritically stolen from PAUP

part *newPart(int nTax, int nChar, char *equateSymbols, int nEquates, char *symbols, int dim)
{
  part *thePart;
  int i, j;

  //printf("got nTax %i, nChar %i, equateSymbols %s, nEquates %i, symbols %s, dim %i\n", 
  //        nTax, nChar, equateSymbols, nEquates, symbols, dim);
  thePart = (part *)malloc(sizeof(part));
	
  thePart->dim = dim;
  thePart->nTax = nTax;
  thePart->nChar = nChar;

  thePart->symbols = calloc((dim + 1), sizeof(char));
  for(i = 0; i < dim; i++) {
    thePart->symbols[i] = symbols[i];
  }

  thePart->equateSymbols = calloc((nEquates + 1), sizeof(char));
  for(i = 0; i < nEquates; i++) {
    thePart->equateSymbols[i] = equateSymbols[i];
  }

  //thePart->symbols = symbols;
  //thePart->equateSymbols = equateSymbols;
  //printf("part.c: newPart: equateSymbols are %s\n", thePart->equateSymbols);
  thePart->patterns = pimatrix(nTax, nChar);
  thePart->patternCounts = malloc(nChar * sizeof(int));
  thePart->sequencePositionPatternIndex = malloc(nChar * sizeof(int));
  thePart->nPatterns = 0;
  thePart->sequences = pimatrix(nTax, nChar);
  thePart->nEquates = nEquates;
  if(nEquates > 0) {
    thePart->equates = pimatrix(nEquates, dim);
  } else {
    thePart->equates = NULL;
  }
  for(i = 0; i < nTax; i++) {
    for(j = 0; j < nChar; j++) {
      thePart->sequences[i][j] = 0;
    }
  }
  if(nEquates > 0) {
    for(i = 0; i < nEquates; i++) {
      for(j = 0; j < dim; j++) {
	thePart->equates[i][j] = 0;
      }
    }
  }
  thePart->globalInvarSitesVec = NULL;
  thePart->globalInvarSitesArray = NULL;
	
  //thePart->siteLikes = malloc(nChar * sizeof(double));
  //for(i = 0; i < nChar; i++) {
  //	thePart->siteLikes[i] = 0.0;
  //}
  thePart->siteLikes = NULL;

  thePart->taxList = malloc(thePart->nTax * sizeof(int));
  thePart->simCats = NULL;
  return thePart;

}

void freePart(part *thePart)
{

  //printf("part.c: freePart\n"); fflush(stdout);
  //checkpoint("  about to free part");
  if(thePart->patterns) free_pimatrix((int **)thePart->patterns);
  if(thePart->patternCounts) free((int *)thePart->patternCounts);
  if(thePart->sequencePositionPatternIndex) free((int *)thePart->sequencePositionPatternIndex);
  if(thePart->sequences) free_pimatrix((int **)thePart->sequences);
  if(thePart->equates) free_pimatrix((int **)thePart->equates);
  if(thePart->equateSymbols) free(thePart->equateSymbols);
  if(thePart->symbols) free(thePart->symbols);

  thePart->patterns = NULL;
  thePart->patternCounts = NULL;
  thePart->sequencePositionPatternIndex = NULL;
  thePart->sequences = NULL;
  thePart->equates = NULL;
  thePart->equateSymbols = NULL;
  thePart->symbols = NULL;

  if(thePart->globalInvarSitesVec) free((int *)thePart->globalInvarSitesVec);
  thePart->globalInvarSitesVec = NULL;
  if(thePart->globalInvarSitesArray) free_pimatrix((int **)thePart->globalInvarSitesArray);
  thePart->globalInvarSitesArray = NULL;
  if(thePart->siteLikes) free(thePart->siteLikes);
  thePart->siteLikes = NULL;
	
  free(thePart->taxList);
  thePart->taxList = NULL;

  if(thePart->simCats) {
    free(thePart->simCats);
    thePart->simCats = NULL;
  }
  //printf("part.c: free-ing part %li\n", (long int)thePart); fflush(stdout);
  free(thePart);
  thePart = NULL;
  //checkpoint("      finished freeing part");
  //printf("part.c: finished freePart\n"); fflush(stdout);


}


void pokeSequences(part *thePart, char *theString)
{
  int i, j, k, m, isCoded;
  char	theChar;
  char   *charsLikeN;
  int     nCharsLikeN = 0;
	
  // Here we want a list of all characters that are fully ambiguous,
  // ie are N_LIKE.  That would include gaps, '?', and (for DNA)
  // 'n'.  They can generally be treated the same.
  charsLikeN = (char *)malloc((thePart->nEquates + 2) * sizeof(char));
  if(!charsLikeN) {
    printf("Failed to malloc charsLikeN.\n");
    exit(0);
  }
  for(i = 0; i < thePart->nEquates + 2; i++) {
    charsLikeN[i] = 'x';
  }
  for(i = 0; i < thePart->nEquates; i++) {
    m = 1; // m will mean it is N_LIKE.  Assume true to start.
    for(j = 0; j < thePart->dim; j++) {
      if(!thePart->equates[i][j]) {
	m = 0;
	break;
      }
    }
    if(m) {
      charsLikeN[nCharsLikeN] = thePart->equateSymbols[i];
      nCharsLikeN++;
    }
  }
  charsLikeN[nCharsLikeN] = '-';
  nCharsLikeN++;
  charsLikeN[nCharsLikeN] = '?';
  nCharsLikeN++;

#if 0
  printf("%i chars are n-like: ", nCharsLikeN);
  for(i = 0; i < thePart->nEquates + 2; i++) {
    printf("%c ", charsLikeN[i]);
  }
  printf("\n");
#endif


  //printf("pokeSequences here\n");
  //printf("the string is:\n%s\n", theString);
  //dumpPart(thePart);
  if(!thePart->sequences) {
    printf("pokeSequences: no thePart->sequences!  Bad!\n");
    return;
  }
  k = 0;
  for(i = 0; i < thePart->nTax; i++) {
    for(j = 0; j < thePart->nChar; j++) {
      theChar = theString[k];
      //printf("theString[%i] = %c\n", k, theString[k]);
      isCoded = 0;
      // is it one of the symbols?
      for(m = 0; m < thePart->dim; m++) {
	//printf("m is %i, testing %c\n", m, thePart->symbols[m]);
	if(theChar == thePart->symbols[m]) {
	  thePart->sequences[i][j] = m;
	  isCoded = 1;
	  k++;
	  break;
	}
      }

#if 0
      if(!isCoded) {
	// is it N-Like?
	for(m = 0; m < nCharsLikeN; m++) {
	  if(theChar == charsLikeN[m]) {
	    thePart->sequences[i][j] = N_LIKE;
	    k++;
	    isCoded = 1;
	    break;
	  }
	}
      }
      if(!isCoded) {
	// is it an equate?
	for(m = 0; m < thePart->nEquates; m++) {
	  //printf("m is %i, testing %c\n", m, thePart->equateSymbols[m]);
	  if(theChar == thePart->equateSymbols[m]) {
	    thePart->sequences[i][j] = EQUATES_BASE + m;
	    isCoded = 1;
	    k++;
	    break;
	  }
	}
      }
      if(!isCoded) {
	printf("part.c pokeSequences.\n");
	printf("    Got character %c.  This should not happen.\n", theChar);
	printf("    It is not in the symbols nor equates.\n");
	exit(0);
      }
#endif	

#if 1
      // if it has not been coded yet, look into the alternatives...
      if(!isCoded && theChar == '-') {
	thePart->sequences[i][j] = GAP_CODE;
	k++;
	isCoded = 1;
      }
      else if(!isCoded && theChar == '?') {
	thePart->sequences[i][j] = QMARK_CODE;
	k++;
	isCoded = 1;
      } else if(!isCoded) {
	//printf("checking whether %c is an equate.\n", theChar);
	for(m = 0; m < thePart->nEquates; m++) {
	  //printf("m is %i, testing %c\n", m, thePart->equateSymbols[m]);
	  if(theChar == thePart->equateSymbols[m]) {
	    thePart->sequences[i][j] = EQUATES_BASE + m;
	    isCoded = 1;
	    k++;
	    break;
	  }
	}
      }
      if(!isCoded) {
	printf("part.c pokeSequences.\n");
	printf("    Got character %c.  This should not happen.\n", theChar);
	printf("    It is not in the symbols nor in the equates.\n");
	exit(0);
      }	
#endif
    }
  }
					

	
#if 0
  for(i = 0; i < thePart->nTax; i++) {
    for(j = 0; j < thePart->nChar; j++) {
      printf("%4i", thePart->sequences[i][j]);
    }
    printf("\n");
  }
#endif

  free(charsLikeN);
  charsLikeN = NULL;

}



void pokeEquatesTable(part *thePart, char *theString)
{
  int i, j, k;
  char	theChar;
	
	
  //printf("pokeEquatesTable here\n");
  //printf("the string is:\n%s\n", theString);
  //dumpPart(thePart);
  if(!thePart->equates) {
    //printf("pokeEquatesTable: no thePart->equates!  Bad!\n");
    return;
  }

  k = 0;
  for(i = 0; i < thePart->nEquates; i++) {
    for(j = 0; j < thePart->dim; j++) {
      theChar = theString[k];
      //printf("theString[%i] = %c\n", k, theString[k]);
      if(theChar == '1') {
	thePart->equates[i][j] = 1;
      }
      k++;
    }
  }
	
  /*
    for(i = 0; i < thePart->nEquates; i++) {
    for(j = 0; j < thePart->dim; j++) {
    printf("%1i", thePart->equates[i][j]);
    }
    printf("\n");
    }
  */
  //dumpPart(thePart);

}

void makePatterns(part *thePart)
{
  int i, j;
  int	*slice;
  int appendPos;
  int	already, i2;
	

  // #####################
  // i for nChar; j for nTax
  // the arrays are eg sequences[j][i], slice[j], patternCounts[i]
  // #####################

  //printf("part.c: starting makePatterns:\n");

  // check for sequences and patterns
  if(!thePart->sequences) {
    printf("makePatterns: no sequences!  bad!\n");
    return;
  }
  if(!thePart->patterns) {
    printf("makePatterns: no memory allocated for patterns!  bad!\n");
    return;
  }
  if(!thePart->patternCounts) {
    printf("makePatterns: no memory allocated for patternCounts!  bad!\n");
    return;
  }
  if(!thePart->sequencePositionPatternIndex) {
    printf("makePatterns: no memory allocated for sequencePositionPatternIndex!  bad!\n");
    return;
  }
		
  // alloc and init
  for(i = 0; i < thePart->nChar; i++) {
    thePart->patternCounts[i] = 0;
    thePart->sequencePositionPatternIndex[i] = 0;
    //for(j = 0; j < thePart->nTax; j++) {
    //	thePart->sequences[j][i] = 0;
    //}
  }
  slice = malloc(thePart->nTax * sizeof(int));
	
  // add the first slice of the alignment to patterns, at sequence position 0
  for(j = 0; j < thePart->nTax; j++) {
    thePart->patterns[j][0] = thePart->sequences[j][0];
  }
  thePart->patternCounts[0] = 1;
  thePart->sequencePositionPatternIndex[0] = 0;

  // main loop
  appendPos = 1;  // cuz I've already got one slice
  for(i = 1; i < thePart->nChar; i++) {
    for(j = 0; j < thePart->nTax; j++) {
      slice[j] = thePart->sequences[j][i];
      //if(slice[j] >= thePart->dim) {
      //	printf("part.makePatterns(): Got bad character %i\n", slice[j]);
      //	printf("charPos = %i, seqNum = %i", i, j);
      //	exit(0);
      //}
    }
    // check to see if already exists
    already = -1;
    for(i2 = 0; i2 < appendPos; i2++) {
      for(j = 0; j < thePart->nTax; j++) {
	if(slice[j] == thePart->patterns[j][i2]) {
	  // if it goes to the end of the taxa without 'break'-ing
	  if(j == thePart->nTax - 1) {
	    already = i2;
	  }
	} else {			// if there is a difference, go to the next j2
	  break;
	}
      }
      // if I have found a column that already exists, stop looking
      if(already != -1) break;
    }

    if(already != -1) {	// if I have found a column that already exists, increment it's count
      thePart->patternCounts[already]++;
      thePart->sequencePositionPatternIndex[i] = already;
    } else {	// otherwize its a new pattern, so add it to the list
      for(j = 0; j < thePart->nTax; j++) {
	thePart->patterns[j][appendPos] = slice[j];
      }
      thePart->patternCounts[appendPos]++;
      thePart->sequencePositionPatternIndex[i] = appendPos;
      appendPos++;
    }
  }

  thePart->nPatterns = appendPos;
  free(slice);

  //printf("\npatterns (%i):\n", thePart->nPatterns);
  for(j = 0; j < thePart->nTax; j++) {
    for(i = 0; i < thePart->nPatterns; i++) {
      //printf("%3i", thePart->patterns[j][i]);
      if(thePart->patterns[j][i] >= thePart->dim) {
	printf("bad character %i\n", thePart->patterns[j][i]);
	exit(0);
      }
    }
    //printf("\n");
    //printf(".");
  }
#if 0
  printf("patternCounts:\n");
  already = 0; // a counter for patternCounts
  for(i = 0; i < thePart->nPatterns; i++) {
    printf("%5i", thePart->patternCounts[i]);
    already = already + thePart->patternCounts[i];
  }
  printf("\n");
  printf("Sum of patternCounts is %i\n", already);
	
  if(thePart->nPatterns > 92) {
    printf("Position 92:\n");
    for(j = 0; j < thePart->nTax; j++) {
      printf("%3i", thePart->patterns[j][92]);
    }
  }
  //printf("sequencePositionPatternIndex:\n");
  //for(i = 0; i < thePart->nChar; i++) {
  //	printf("%3i:%3i,  ", i, thePart->sequencePositionPatternIndex[i]);
  //}
  printf("\n");
#endif

  //printf("part.c: finished makePatterns:  nPatterns = %i\n", thePart->nPatterns);

}



void dumpPart(part *thePart)
{
  int max, i, j;
	
  printf("Dump (c) part\n");
  printf("    dim = %i\n", thePart->dim);
  printf("    nTax = %i\n", thePart->nTax);
  printf("    nChar = %i\n", thePart->nChar);
  printf("    symbols = %s\n", thePart->symbols);
  printf("    patterns = %li\n", (long int)thePart->patterns);
  printf("    patternCounts = %li\n", (long int)thePart->patternCounts);
  printf("    sequencePositionPatternIndex = %li\n", (long int)thePart->sequencePositionPatternIndex);
  printf("    nPatterns = %i\n", thePart->nPatterns);
  printf("    nEquates = %i\n", thePart->nEquates);
  printf("    equateSymbols = %s\n", thePart->equateSymbols);
  printf("    sequences = %li\n", (long int)thePart->sequences);
  if(thePart->sequences) {
    if(thePart->nChar >= 20) {
      max = 20;
    } 
    else {
      max = thePart->nChar;
    }
    printf("    sequences:\n");
    for(i = 0; i < thePart->nTax; i++) {
      printf("      sequences[%i]:", i);
      for(j = 0; j < max; j++) {
	printf("%4i", thePart->sequences[i][j]);
      }
      if(max < thePart->nChar) {
	printf("...");
      }
      printf("\n");
    }
  }
  if(thePart->nPatterns) {
    if(thePart->nPatterns >= 20) {
      max = 20;
    } 
    else {
      max = thePart->nPatterns;
    }
    printf("    patterns:\n");
    for(i = 0; i < thePart->nTax; i++) {
      printf("      patterns[%i]: ", i);
      for(j = 0; j < max; j++) {
	printf("%4i", thePart->patterns[i][j]);
      }
      if(max < thePart->nPatterns) {
	printf("...");
      }
      printf("\n");
    }
  }


  if(thePart->nEquates) {
    printf("    equates:\n");
    for(i = 0; i < thePart->nEquates; i++) {
      printf("      equates[%i]: ", i);
      for(j = 0; j < thePart->dim; j++) {
	printf("%1i", thePart->equates[i][j]);
      }
      printf("\n");
    }
  }


  if(thePart->globalInvarSitesVec) {
    printf("    globalInvarSitesVec:\n");
    printf("      ");
    for(i = 0; i < thePart->nPatterns; i++) {
      printf("%4i", thePart->globalInvarSitesVec[i]);
    }
    printf("\n");
  }

  if(thePart->globalInvarSitesArray) {
    if(thePart->nPatterns >= 20) {
      max = 20;
    } 
    else {
      max = thePart->nPatterns;
    }
    printf("    globalInvarSitesArray:\n");
    for(i = 0; i < thePart->dim; i++) {
      printf("      symbols[%i] = %c : ", i,  thePart->symbols[i]);
      for(j = 0; j < max; j++) {
	printf("%1i", thePart->globalInvarSitesArray[i][j]);
      }
      if(max < thePart->nPatterns) {
	printf("...");
      }
      printf("\n");
    }
  }


}





PyObject *singleSequenceBaseCounts(part *thePart, int sequenceNum)
{
  int i;
  int	*totals = NULL;
  PyObject	*thePyList;
	
  totals = pivector(thePart->dim);
  for(i = 0; i < thePart->dim; i++) {
    totals[i] = 0;
  }
	


  if(thePart->nPatterns) {
    for(i = 0; i < thePart->nPatterns; i++) {
      if(thePart->patterns[sequenceNum][i] >= 0) {
	totals[thePart->patterns[sequenceNum][i]] = 
	  totals[thePart->patterns[sequenceNum][i]] + thePart->patternCounts[i];
      }
    }
  }

  else if(thePart->nChar) {
    for(i = 0; i < thePart->nChar; i++) {
      if(thePart->sequences[sequenceNum][i] >= 0) {
	totals[thePart->sequences[sequenceNum][i]] = 
	  totals[thePart->sequences[sequenceNum][i]] + 1;
      }
    }
  }
  else {
    printf("part: singleSequenceBaseCounts: no sequences.\n");
    exit(1);
  }
			
  thePyList = PyList_New(thePart->dim);
  for(i = 0; i < thePart->dim; i++) {
    //PyList_SetItem(thePyList, i, PyFloat_FromDouble(totals[i]));
    PyList_SetItem(thePyList, i, PyInt_FromLong((long int)totals[i]));
  }
	
  free(totals);
  totals = NULL;
  return thePyList;

}



PyObject *symbolSequences(part *thePart)
{
  int i, j;
  char *oneLongSeq;
  int theCharNum;
  int spot;
  PyObject *theOneLongSeqPy;

#if 0
  printf("part: symbolSequences: \n");
  printf("part: %li \n", (long int)thePart);
  //printf("length = %i\n", thePart->nChar);
  printf("first bit of sequence: ");
  for(i = 0; i < 15; i++) {
    printf("%3i", thePart->sequences[0][i]);
  }
  printf("...\n");
#endif

  oneLongSeq = calloc(((thePart->nTax * thePart->nChar) + 1), sizeof(char));
  if(!oneLongSeq) {
    printf("part: symbolSequences: failed to malloc totals array\n");
    exit(1);
  }

  //oneLongSeq[thePart->nTax * thePart->nChar] = '\0';

#if 0
  if(thePart->nPatterns) {
    for(i = 0; i < thePart->nPatterns; i++) {
      if(thePart->patterns[sequenceNum][i] >= 0) {
	totals[thePart->patterns[sequenceNum][i]] = 
	  totals[thePart->patterns[sequenceNum][i]] + thePart->patternCounts[i];
      }
      else {
	// put ambiguity code stuff here
      }
    }
  }
  else if(thePart->nChar) {
#endif

    if(thePart->nChar) {
      spot = 0;
      for(i = 0; i < thePart->nTax; i++) {
	for(j = 0; j < thePart->nChar; j++) {
	  theCharNum = thePart->sequences[i][j];
	  //printf("%3i", theCharNum);
	  if(theCharNum >= 0) {
	    oneLongSeq[spot] = thePart->symbols[theCharNum];
	  } else if (theCharNum == GAP_CODE) {
	    oneLongSeq[spot] = '-';
	  } else if (theCharNum == QMARK_CODE) {
	    oneLongSeq[spot] = '?';
	  } else if (theCharNum >= EQUATES_BASE && 
		     theCharNum < EQUATES_BASE + thePart->nEquates) {
	    oneLongSeq[spot] = thePart->equateSymbols[theCharNum - EQUATES_BASE];
	  } else {
	    printf("part.c: symbolSequences: character number %i is not recognized.\n",
		   theCharNum);
	    exit(1);
	  }
	  spot++;
	}
      }
    }

    else {
      printf("part.c: symbolSequences: no sequences.\n");
      exit(1);
    }
    //printf("\n");
#if 0
    printf("oneLongSeq = %s\n", oneLongSeq);
#endif

    theOneLongSeqPy = Py_BuildValue("s", oneLongSeq);
    free(oneLongSeq);	
    return theOneLongSeqPy;

  }



  double unconstrainedLogLike(part *thePart)
  {
    int i, j;
    double dsum;

    // check for silliness...
    if(!thePart->nPatterns) {
      printf("part.c: unconstrainedLogLike: no patterns.\n");
      printf("   Needs both sequences and patterns.\n");
      exit(1);
    }
	
    // check for bad characters...
    for(i = 0; i < thePart->nTax; i++) {
      for(j = 0; j < thePart->nChar; j++) {
	if(thePart->sequences[i][j] < 0) {
	  printf("part.c: unconstrainedLogLike: bad character.\n");
	  printf("   Can\'t do this calculation if there are any\n");
	  printf("   gaps, unknowns, ambiguities, or equates, ok?\n");
	  exit(1);
	}
      }
    }

    dsum = 0.0;
    for(i = 0; i < thePart->nPatterns; i++) {
      dsum = dsum + (thePart->patternCounts[i] * log((double)(thePart->patternCounts[i])));
    }
    dsum = dsum - (thePart->nChar * log((double)(thePart->nChar)));
	
    return dsum;

  }

  void setGlobalInvarSitesVec(part *thePart)
  {
    // use the globalInvarSitesVec to hold part constant sites mapping
    // If its a small positive int value then it is invar,  If its variable, its Zero.
    // Note that it is for the patterns, not the sequences.
    // The Array holds more info on how the sites are invar.  Eg for DNA
    // if it is 1000 then it is invar A.  If it is 1010 it is invar a and g. (ie r).

    int p, t;  // pos, tax
    int i, j;
    int tester[thePart->nTax][thePart->dim];
    int result[thePart->dim];
    int theSum;

    // create and initialize the vector and array
    if(!thePart->globalInvarSitesVec) {
      thePart->globalInvarSitesVec = (int *)malloc(sizeof(int) * thePart->nChar);
      for(p = 0; p < thePart->nChar; p++) {
	thePart->globalInvarSitesVec[p] = 0;  // to initialize
      }
    }
    if(!thePart->globalInvarSitesArray) {
      thePart->globalInvarSitesArray = pimatrix(thePart->dim, thePart->nChar);
      for(i = 0; i < thePart->dim; i++) {
	for(j = 0; j < thePart->nChar; j++) {
	  thePart->globalInvarSitesArray[i][j] = 0;
	}
      }
    }

    // Heres the scoop.  If all the sites are the same and >= zero,
    // then of course it is an invar site.
    // However, if there are gaps or equates, then it may be an invar site
    // even though all the sites are not the same, if they are compatible
    // equates.  Eg the site raa is invar, but yaa is not.
    // Here's a tough one: the site hra (where h is a, c, or t) is invar.
    // And of course hrg is a varied site.
    //
    // So what I do is make a tester, an array of ints, which I fill for each site.
    // For example, for hra, the tester would be
    //      1101   for h
    //      1010   for r
    //      1000   for a
    //
    // Then I can get a result vector, in this case 1000.  
    // This result vector goes in the globalInvarSitesArray.
    // And the total of the 1's in that result (in this case 1) goes
    // in the globalInvarSitesVector.
    // So if its a zero, then its a varied site,

    for(p = 0; p < thePart->nPatterns; p++) {
      for(t = 0; t < thePart->nTax; t++) {
	//printf("part = %i, tax = %i, patterns[tax][pat] = %i\n", 
	//	   p, t, thePart->patterns[t][p]);
	// if its a regular symbol character...
	if(thePart->patterns[t][p] >= 0) {
	  for(j = 0; j < thePart->dim; j++) {
	    tester[t][j] = 0;
	  }
	  tester[t][thePart->patterns[t][p]] = 1;
	} else if(thePart->patterns[t][p] == N_LIKE || 
		  thePart->patterns[t][p] == GAP_CODE || 
		  thePart->patterns[t][p] == QMARK_CODE) {
	  for(j = 0; j < thePart->dim; j++) {
	    tester[t][j] = 1;
	  }
	} else {
	  for(j = 0; j < thePart->dim; j++) {
	    tester[t][j] = thePart->equates[thePart->patterns[t][p] - EQUATES_BASE][j];
	  }
	}
      }
      //for(i = 0; i < thePart->nTax; i++) {
      //	for(j = 0; j < thePart->dim; j++) {
      //		printf("%4i", tester[i][j]);
      //	}
      //	printf("\n");
      //}
      for(j = 0; j < thePart->dim; j++) {
	result[j] = 1;
	for(i = 0; i < thePart->nTax; i++) {
	  if(!tester[i][j]) {
	    result[j] = 0;
	    break;
	  }
	}
      }
      //printf("result:  ");
      //for(j = 0; j < thePart->dim; j++) {
      //	printf("%1i", result[j]);
      //}
      //printf("\n");
      theSum = 0;
      for(j = 0; j < thePart->dim; j++) {
	theSum = theSum + result[j];
      }
      //printf("The sum: %i\n", theSum);
      thePart->globalInvarSitesVec[p] = theSum;
      for(j = 0; j < thePart->dim; j++) {
	thePart->globalInvarSitesArray[j][p] = result[j];
      }
    }
	


    if(0){
      printf("globalInvarSitesVec.\n");
      printf("    A small positive int means possibly constant.\n");  
      printf("    Zero means varied\n");
      for(p = 0; p < thePart->nPatterns; p++) {
	printf("%1i", thePart->globalInvarSitesVec[p]);
      }
      printf("\n\nglobalInvarSitesArray:\n");
      for(i = 0; i < thePart->dim; i++) {
	for(j = 0; j < thePart->nPatterns; j++) {
	  printf("%1i", thePart->globalInvarSitesArray[i][j]);
	}
	printf("\n");
      }

    }
    if(0){
      int nInvarSites = 0;

      for(p = 0; p < thePart->nPatterns; p++) {
	if(thePart->globalInvarSitesVec[p]) {
	  nInvarSites += thePart->patternCounts[p];
	}
      }
      printf("part.c setGlobalInvarSitesVec().  nInvarSites=%i\n", nInvarSites);
    }

  }

  void partCompositionC(part *thePart, double *theComp)

  // This can get comps of single, or multiple, sequences, as well
  // as the whole lot.  Which sequences are measured is determined
  // by the thePart->taxList vector.

  {
    int      i, j, k;
    int      hasEquates = 0;
    double  *symbolFreq = NULL;
    double  *equateFreq = NULL;
    double	*comp = NULL;
    double	*symbSum = NULL;
    double   epsilon = 1.0e-12;
    int      maxIterations = 1000;
    double   x = 0.0;
    double   diff = 0.0;
    double   oldComp = 0.0;
    double  *results = NULL;
    int      nSites=0;
    int      grandNSites = 0;
    int      nGapsMissings = 0;
    int      dummy = 0;
    int      seqNum = 0;

#if 0
    printf("part.c: partCompositionC()");
    for(i = 0; i < thePart->nTax; i++) {
      printf(" %i ", thePart->taxList[i]);
    }
    printf("\n");
#endif

    symbolFreq = pdvector(thePart->dim);
    comp = pdvector(thePart->dim);
    symbSum = pdvector(thePart->dim);
    results = pdvector(thePart->dim);
    for(i = 0; i < thePart->dim; i++) {
      symbolFreq[i] = 0.0;
      comp[i] = 0.0;
      symbSum[i] = 0.0;
      results[i] = 0.0;
    }
    equateFreq = pdvector(thePart->nEquates);
    for(i = 0; i < thePart->nEquates; i++) {
      equateFreq[i] = 0.0;
    }

    //printf("equates: %s\n", thePart->equateSymbols);
    //printf("taxList: ");
    //for(i = 0; i < thePart->nTax; i++) {
    //	printf(" %2i", thePart->taxList[i]);
    //}
    //printf("\n");

    for(seqNum = 0; seqNum < thePart->nTax; seqNum++) {
      if(thePart->taxList[seqNum]) {
	nGapsMissings = 0;

	for(k = 0; k < thePart->dim; k++) {
	  symbolFreq[k] = 0.0;
	  comp[k] = 0.0;
	  symbSum[k] = 0.0;
	}
	for(k = 0; k < thePart->nEquates; k++) {
	  equateFreq[k] = 0.0;
	}

	for(j = 0; j < thePart->nChar; j++) {
	  //printf("%4i", thePart->sequences[i][j]);
	  if(thePart->sequences[seqNum][j] >= 0) {
	    symbolFreq[thePart->sequences[seqNum][j]] = 
	      symbolFreq[thePart->sequences[seqNum][j]] + 1.0;
	  }
	  else if(thePart->sequences[seqNum][j] == QMARK_CODE || 
		  thePart->sequences[seqNum][j] ==  GAP_CODE) {
	    nGapsMissings += 1;
	  }
	  else { // an equate
	    //printf("(%i)", thePart->sequences[i][j] - EQUATES_BASE);
	    equateFreq[thePart->sequences[seqNum][j] - EQUATES_BASE] += 1.0; 
	  }
	}
	nSites = thePart->nChar - nGapsMissings;
	grandNSites += nSites;
	//printf("nSites = %i, grandNSites=%i\n", nSites, grandNSites);


	for(i = 0; i < thePart->nEquates; i++) {
	  if(equateFreq[i] > 0.0) {
	    hasEquates = 1;
	    break;
	  }
	}
	
#if 0
	printf("symbolFreq, ");
	for(j = 0; j < thePart->dim; j++) {
	  printf(" %8.5f", symbolFreq[j]);
	}
	printf("\n");
	printf("equateFreq, ");
	for(j = 0; j < thePart->nEquates; j++) {
	  printf(" %8.5f", equateFreq[j]);
	}
	printf("\n");
#endif

	// There is no point unless we have nSites
	if(nSites) {

	  x = 1.0 / ((float) thePart->dim);
	  for(i = 0; i < thePart->dim; i++) {
	    comp[i] = x;
	  }

	  for(dummy = 0; dummy < maxIterations; dummy++) {
	    for(j = 0; j < thePart->dim; j++) {
	      symbSum[j] = symbolFreq[j];
	    }
	    if(hasEquates) {
	      for(j = 0; j < thePart->nEquates; j++) {
		if(equateFreq[j] > 0.0) {
		  x = 0.0;
		  for(k = 0; k < thePart->dim; k++) {
		    if(thePart->equates[j][k]) {
		      x = x + comp[k];
		    }
		  }
		  for(k = 0; k < thePart->dim; k++) {
		    if(thePart->equates[j][k]) {
		      symbSum[k] = symbSum[k] + (equateFreq[j] * (comp[k] / x));
		    }
		  }
		}
	      }
	    }
	    x = 0.0;
	    for(j = 0; j < thePart->dim; j++) {
	      x = x + symbSum[j];
	    }
	    diff = 0.0;
	    for(j = 0; j < thePart->dim; j++) {
	      oldComp = comp[j];
	      comp[j] = symbSum[j] / x;
	      diff = diff + fabs(comp[j] - oldComp);
	    }
	    //printf("iteration %i, diff=%8.5f, ", i, diff);
	    //for(j = 0; j < thePart->dim; j++) {
	    //	printf(" %8.5f", comp[j]);
	    //}
	    //printf("\n");
	    if(diff < epsilon) {
	      break;
	    }
	  }
	  for(j = 0; j < thePart->dim; j++) {
	    results[j] = results[j] + (comp[j] * (double)nSites);
	  }
	}
      }
    }

#if 0
    printf("grandNSites = %i\n", grandNSites);
    printf("results: ");
    for(j = 0; j < thePart->dim; j++) {
      printf(" %12.5f", results[j] / grandNSites);
    }
    printf("\n");
#endif

    if(grandNSites) {
      for(i = 0; i < thePart->dim; i++) {
	theComp[i] = results[i] / (double)grandNSites;
      }
    }
    else {
      // its a missing sequence.  All comps are zero.
      for(i = 0; i < thePart->dim; i++) {
	theComp[i] = 0.0;
      }
    }

    free(symbolFreq);
    symbolFreq = NULL;
    free(equateFreq);
    equateFreq = NULL;
    free(comp);
    comp = NULL;
    free(symbSum);
    symbSum = NULL;
    free(results);
    results = NULL;
  }



  PyObject *partCompositionP(part *thePart)
  {
    int      i;
    double	*comp = NULL;
    PyObject	*thePyList;

    comp = pdvector(thePart->dim);

    partCompositionC(thePart, comp);
    thePyList = PyList_New(thePart->dim);
    for(i = 0; i < thePart->dim; i++) {
      PyList_SetItem(thePyList, i, PyFloat_FromDouble(comp[i]));
    }

    free(comp);
    comp = NULL;
    return thePyList;
  }


  int partSequenceSitesCount(part *thePart, int sequenceNum)  // ie without gaps or missings
  {

    int nGapsMissings=0;
    int i;

    for(i = 0; i < thePart->nChar; i++) {
      if(thePart->sequences[sequenceNum][i] == GAP_CODE ||
	 thePart->sequences[sequenceNum][i] == QMARK_CODE) {
	nGapsMissings++;
      }
    }
    return thePart->nChar - nGapsMissings;
  }

	
  void calcEmpiricalRMatrixViaMatrixLog(part *thePart)
  {
    int       i, j, k;
    int       ch1, ch2;
    int       total = 0;
    double    factor = 0.0;
    int     **bigF = NULL;
    double  **bigP = NULL;
    double  **bigQR = NULL;
    double   *comp = NULL;
    eig      *theEig = NULL;

    bigF = psimatrix(thePart->dim);
    bigP = psdmatrix(thePart->dim);
    bigQR = psdmatrix(thePart->dim);
    comp = pdvector(thePart->dim);

    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	bigF[i][j] = 0;
	bigP[i][j] = 0.0;
	bigQR[i][j] = 0.0;
      }
      comp[i] = 0.0;
    }

    // I want the composition for all sequences
    for(i = 0; i < thePart->nTax; i++) {
      thePart->taxList[i] = 1;
    }

    // Get the composition
    partCompositionC(thePart, comp);
    printf("Composition:\n");
    for(i = 0; i < thePart->dim; i++) {
      printf("    %g", comp[i]);
    }
    printf("\n");

    // fill the bigF matrix, the matrix of changes
    for(i = 0; i < (thePart->nTax - 1); i++) {
      for(j = i + 1; j < thePart->nTax; j++) {
	for(k = 0; k < thePart->nChar; k++) {  // it would be more efficient if this was the outside loop
	  //printf("%i %i %i:  %4i %4i\n", i, j, k, thePart->sequences[i][k], thePart->sequences[j][k]);
	  ch1 = thePart->sequences[i][k];
	  ch2 = thePart->sequences[j][k];
	  if(ch1 < QMARK_CODE || ch2 < QMARK_CODE) {
	    printf("calcEmpiricalRMatrixViaMatrixLog()\n");
	    printf("    This doesn\'t work if there are ambiguities.\n");
	    exit(0);
	  }
	  if(ch1 < 0 || ch2 < 0) {
	    // pass
	  } else {
	    bigF[ch1][ch2] = bigF[ch1][ch2] + 1;
	    bigF[ch2][ch1] = bigF[ch2][ch1] + 1;
	  }
	}
      }
    }

    printf("bigF:\n");
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	printf("  %5i", bigF[i][j]);
      }
      printf("\n");
    }
			
    for(i = 0; i < thePart->dim; i++) {
      total = 0;
      for(j = 0; j < thePart->dim; j++) {
	total = total + bigF[i][j];
      }
      for(j = 0; j < thePart->dim; j++) {
	bigP[i][j] = (double)bigF[i][j] / (double)total;
      }
    }

    printf("bigP (ced):\n");
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	printf("  %10.7f", bigP[i][j]);
      }
      printf("\n");
    }

    printf("Sums of the rows\n");
    for(i = 0; i < thePart->dim; i++) {
      factor = 0.0;
      for(j = 0; j < thePart->dim; j++) {
	factor = factor + bigP[i][j];
      }
      printf("  %10.7f", factor);
    }
    printf("\n");



    theEig = allocEig(thePart->dim, bigP);
    matrixLog(theEig, bigQR);
	
    printf("loggified bigP (before normalizing):\n");
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	printf("  %10.7f", bigQR[i][j]);
      }
      printf("\n");
    }


    // sum the off diagonal elements * pi(row)
    factor = 0.0;
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	if (i != j) {
	  factor = factor + (bigQR[i][j] * comp[i]);
	}
      }
    }
    printf("The sum of the off diagonal elements * comp(row) is %f\n", factor);
	
    normalizeBigQ(bigQR, comp, thePart->dim);

    printf("bigQ (ie loggified bigP after normalizing):\n");
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	printf("  %10.7f", bigQR[i][j]);
      }
      printf("\n");
    }

    // sum the off diagonal elements (again)
    factor = 0.0;
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	if (i != j) {
	  factor = factor + (bigQR[i][j] * comp[i]);
	}
      }
    }
    printf("The sum of the off diagonal elements * comp(row) is %f\n", factor);

    printf("Sums of the rows\n");
    for(i = 0; i < thePart->dim; i++) {
      factor = 0.0;
      for(j = 0; j < thePart->dim; j++) {
	factor = factor + bigQR[i][j];
      }
      printf("  %10.7f", factor);
    }
    printf("\n");


    // calculate bigR, by dividing the columns by the composition
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	bigQR[i][j] = bigQR[i][j] / comp[j];
      }
    }

    printf("bigR:\n");
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	if(i == j) {
	  printf("   -        ");
	} else {
	  printf("  %10.7f", bigQR[i][j]);
	}
      }
      printf("\n");
    }


    factor = 1.0 / bigQR[thePart->dim - 2][thePart->dim - 1];
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	bigQR[i][j] = bigQR[i][j] * factor;
      }
    }

    printf("bigR (scaled):\n");
    for(i = 0; i < thePart->dim; i++) {
      for(j = 0; j < thePart->dim; j++) {
	if(i == j) {
	  printf("   -        ");
	} else {
	  printf("  %10.7f", bigQR[i][j]);
	}
      }
      printf("\n");
    }

	
    free_psimatrix(bigF);
    bigF = NULL;
    free_psdmatrix(bigP);
    bigP = NULL;
    free_psdmatrix(bigQR);
    bigQR = NULL;
    free(comp);
    comp = NULL;
    freeEig(theEig);
    theEig = NULL;
  }





  double steelCRCInvariants(part *thePart)
  {
    // From Thollesson's LDDist.  Thanks Mikael!  Estimates the
    // fraction of constant sites that are invariant by
    // Capture-Recapture according to Steel et al.  Steel, M., Huson,
    // D. and Lockhart, P.J. (2000) Invariable sites models and their
    // use in phylogeny reconstruction. Syst. Biol., 49, 225-232.
    // Returns the fraction of the constant sites that are also
    // estimated to be invariant



    const int pick = 4;
    const int scale = 2;
    int TaxNumCutOff = 50;
    int RepNums = TaxNumCutOff*(TaxNumCutOff-1)*(TaxNumCutOff-2)*(TaxNumCutOff-3)/8;
    int RepsToDo = RepNums / (scale * (thePart->nTax / pick));
    int replicates = 0;
    double pSt = 0;
    int theTaxa[pick];
    int *theArr;
    int i,j,k,l;
    int tx, CurrentSize,Rep,choice;

    //theArr = new int[TaxNum];

    if(thePart->nTax < 4) {
      printf("part steelCRCInvariants().  Need a minimum of 4 sequences.  Got %i\n", thePart->nTax);
      exit(0);
    }
    theArr = malloc(thePart->nTax * sizeof(int));
    if(!theArr) {
      printf("Failed to allocate memory for theArr\n");
      exit(0);
    }
	
    if (thePart->nTax > TaxNumCutOff) {
      for (Rep=0;Rep<RepsToDo;Rep++) {
	// Initialize the list
	for (i = 0; i< thePart->nTax; i++) {
	  theArr[i] = i;			
	}
	CurrentSize = thePart->nTax - 1;
	
	while (CurrentSize>pick) {
	  // Choose four
	  for (tx=0;tx<pick;tx++) {
	    choice = rand()%(CurrentSize);
	    theTaxa[tx] = theArr[choice];
	    theArr[choice]=theArr[CurrentSize];
	    CurrentSize--;
	  }
	  pSt += qEval(thePart, theTaxa[0],theTaxa[1],theTaxa[2],theTaxa[3]);
	  replicates++;
			
	}
      }
    } else {
      for (l=3;l<thePart->nTax;l++) {
	for (k=2;k<l;k++) {
	  for (j=1;j<k;j++) {
	    for (i=0;i<j;i++) {
	      pSt += qEval(thePart,i,j,k,l);
	      replicates++;
	    }
	  }
	}
      }
    }

    free(theArr);
    theArr = NULL;
    //return vSum/replicates;
    return pSt/replicates;

  }

  double qEval(part *thePart, int i, int j, int k, int l)
  {
    // From Thollesson's LDDist.  Thanks Mikael!  
    int lmij,lmkl,lmik,lmjl,lmil,lmjk; //,lmijkl,lmikjl,lmiljk;
    int total;
    int unvaried;
    int site;
    int Term1, Term2, Term3;
    int MaxTerm;
    double TempSum;
	
    int mij=0;
    int mik=0;
    int mil=0;
    int mjk=0;
    int mjl=0;
    int mkl=0;
	
    int mijkl=0;
    int mikjl=0;
    int miljk=0;

    unvaried=0;
    total=0;
	
    for (site=0;site<thePart->nChar;site++) {
      //if (FullFour(i,j,k,l,site)) {
      if (thePart->sequences[i][site] >= 0 &&
	  thePart->sequences[j][site] >= 0 &&
	  thePart->sequences[k][site] >= 0 &&
	  thePart->sequences[l][site] >= 0) {
	total++;
			
	lmij = (int)(!(thePart->sequences[i][site] == thePart->sequences[j][site]));
	lmik = (int)(!(thePart->sequences[i][site] == thePart->sequences[k][site]));
	lmil = (int)(!(thePart->sequences[i][site] == thePart->sequences[l][site]));
	lmjk = (int)(!(thePart->sequences[j][site] == thePart->sequences[k][site]));
	lmjl = (int)(!(thePart->sequences[j][site] == thePart->sequences[l][site]));
	lmkl = (int)(!(thePart->sequences[k][site] == thePart->sequences[l][site]));

	mij+=lmij;
	mik+=lmik;
	mil+=lmil;
	mjk+=lmjk;
	mjl+=lmjl;
	mkl+=lmkl;
			
	mijkl+= (int)(lmij && lmkl);
	mikjl+= (int)(lmik && lmjl);
	miljk+= (int)(lmil && lmjk);
	unvaried+= (!(lmij||lmik||lmil||lmjk||lmjl||lmkl));
      }
    }

    Term1 = mij*mkl/mijkl;
    Term2 = mik*mjl/mikjl;
    Term3 = mil*mjk/miljk;


	
    MaxTerm = (Term1>Term2 ? (Term1>Term3 ? Term1 : (Term3>Term2 ? Term3 : Term2)) : (Term2>Term3 ? Term2 : Term3));

    TempSum = (MaxTerm>total ? 1.0 : ((MaxTerm>(total-unvaried)) ? (double) ((double) MaxTerm/total) : (double) (1.0-(double)unvaried/total)));
    return ((1.0-TempSum)*total/unvaried);

  }

  void recodeNLike(part *thePart)
  {
    int     i, j, k, isNLike;
    int    *codesLikeN;
    int     nCodesLikeN = 0;
	
    // Here we want a list of all character codes that are fully
    // ambiguous, ie are N_LIKE.  That would include gaps, '?', and
    // (for DNA) 'n'.  They can generally be treated the same.

    codesLikeN = (int *)malloc((thePart->nEquates + 3) * sizeof(int));
    if(!codesLikeN) {
      printf("Failed to malloc codesLikeN.\n");
      exit(0);
    }
    for(i = 0; i < thePart->nEquates + 2; i++) {
      codesLikeN[i] = N_LIKE;
    }
    for(i = 0; i < thePart->nEquates; i++) {
      isNLike = 1; 
      for(j = 0; j < thePart->dim; j++) {
	if(!thePart->equates[i][j]) {
	  isNLike = 0;
	  break;
	}
      }
      if(isNLike) {
	codesLikeN[nCodesLikeN] = EQUATES_BASE + i;
	nCodesLikeN++;
      }
    }
    codesLikeN[nCodesLikeN] = GAP_CODE;
    nCodesLikeN++;
    codesLikeN[nCodesLikeN] = QMARK_CODE;
    nCodesLikeN++;
    codesLikeN[nCodesLikeN] = N_LIKE;
    nCodesLikeN++;

#if 0
    printf("%i n-like codes are: ", nCodesLikeN);
    for(i = 0; i < thePart->nEquates + 2; i++) {
      printf("%3i ", codesLikeN[i]);
    }
    printf("\n");
#endif

    for(i = 0; i < thePart->nTax; i++) {
      for(j = 0; j < thePart->nChar; j++) {
	for(k = 0; k < nCodesLikeN; k++) {
	  if(thePart->sequences[i][j] == codesLikeN[k]) {
	    thePart->sequences[i][j] = N_LIKE;
	    break;
	  }
	}
      }
    }

    free(codesLikeN);
  }


  double partMeanNCharsPerSite(part *thePart)
  {
    int charNum,taxNum,symbNum,nSites,totalThisSite,theTotal;
    int *counts;
	
    counts = (int *)malloc((thePart->dim) * sizeof(int));
    if(!counts) {
      printf("Failed to malloc counts.\n");
      exit(1);
    }

    for(symbNum = 0; symbNum < thePart->dim; symbNum++) {
      counts[symbNum] = 0;
    }
    nSites = 0;
    theTotal = 0;
    for(charNum = 0; charNum < thePart->nChar; charNum++) {
      for(taxNum = 0; taxNum < thePart->nTax; taxNum++) {
	symbNum = thePart->sequences[taxNum][charNum];
	if((symbNum >= 0) && (symbNum < thePart->dim)) {
	  counts[symbNum] += 1;
	}
      }
      totalThisSite = 0;
      for(symbNum = 0; symbNum < thePart->dim; symbNum++) {
	if(counts[symbNum]) {
	  totalThisSite += 1;
	  counts[symbNum] = 0;
	}
      }
      if(totalThisSite > 1) {
	theTotal += totalThisSite;
	nSites += 1;
      }
    }

    free(counts);
    return (double)theTotal / (double)nSites;
  }

  int partSimpleConstantSitesCount(part *thePart)
  {
    int charNum,taxNum,symbNum,topCharNum,taxNum2,isConstant;
    int counts = 0;
	

    for(charNum = 0; charNum < thePart->nChar; charNum++) {
      topCharNum = -1;
      isConstant = 1;
      for(taxNum = 0; taxNum < thePart->nTax; taxNum++) {
	symbNum = thePart->sequences[taxNum][charNum];
	if((symbNum >= 0) && (symbNum < thePart->dim)) {
	  topCharNum = symbNum;
	  break;
	}
      }
      for(taxNum2 = taxNum; taxNum2 < thePart->nTax; taxNum2++) {
	symbNum = thePart->sequences[taxNum2][charNum];
	if((symbNum >= 0) && (symbNum < thePart->dim) && symbNum != topCharNum) {
	  isConstant = 0;
	  break;
	}
      }
      if(topCharNum != -1 && isConstant) {
	counts += 1;
      }
    }

    return counts;
  }

  double partBigXSquared(part *thePart)
  {
    int charNum,taxNum,symbNum;
    double **obs, *exp, *soc;
    double oneOverNTax, xSq;
	
    obs = pdmatrix(thePart->nTax, thePart->dim);
    if(!obs) {
      printf("Failed to malloc obs.\n");
      exit(1);
    }
    exp = pdvector(thePart->dim);
    if(!exp) {
      printf("Failed to malloc exp.\n");
      exit(1);
    }
    soc = pdvector(thePart->dim);
    if(!soc) {
      printf("Failed to malloc soc.\n");
      exit(1);
    }

    // zero them
    for(symbNum = 0; symbNum < thePart->dim; symbNum++) {
      exp[symbNum] = 0.0;
      for(taxNum = 0; taxNum < thePart->nTax; taxNum++) {
	obs[taxNum][symbNum] = 0.0;
      }
    }

    oneOverNTax = 1.0 / thePart->nTax;

    // fill obs
    for(charNum = 0; charNum < thePart->nChar; charNum++) {
      for(taxNum = 0; taxNum < thePart->nTax; taxNum++) {
	symbNum = thePart->sequences[taxNum][charNum];
	if((symbNum >= 0) && (symbNum < thePart->dim)) {
	  obs[taxNum][symbNum] += 1;
	} else {
	  printf("This function cannot handle gaps and ambiguities.  Returning -2.0\n");
	  free(obs);
	  free(exp);
	  return -2.0;
	}
      }
    }
    // fill soc
    for(symbNum = 0; symbNum < thePart->dim; symbNum++) {
      soc[symbNum] = 0.0;
      for(taxNum = 0; taxNum < thePart->nTax; taxNum++) {
	soc[symbNum] += obs[taxNum][symbNum];
      }
      //if(soc[symbNum] == 0.0) {
      //	printf("This function cannot handle zeros in column sums.  Returning -1.0\n");
      //	free(obs);
      //	free(exp);
      //	free(soc);
      //	return -1.0;
      //}
			
    }

    // fill exp
    for(symbNum = 0; symbNum < thePart->dim; symbNum++) {
      exp[symbNum] = oneOverNTax * soc[symbNum] ;
    }
	

    xSq = 0.0;
    for(taxNum = 0; taxNum < thePart->nTax; taxNum++) {
      for(symbNum = 0; symbNum < thePart->dim; symbNum++) {
	if(exp[symbNum]) { // silently skip column zeros
	  xSq += ((obs[taxNum][symbNum] - exp[symbNum]) * (obs[taxNum][symbNum] - exp[symbNum])) /exp[symbNum];
	}
      }
    }

    free(obs);
    free(exp);
    free(soc);
    return xSq;
  }
