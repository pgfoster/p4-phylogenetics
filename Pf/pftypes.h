#include <gsl/gsl_rng.h>

typedef struct dataStruct data;
typedef struct partStruct part;
typedef struct eigStruct eig;
typedef struct pcomplexStruct pcomplex;
typedef struct rellStuffStruct rellStuff;

typedef struct p4_treeStruct p4_tree;
typedef struct p4_nodeStruct p4_node;
typedef struct p4_modelStruct p4_model;
typedef struct p4_modelPartStruct p4_modelPart;
typedef struct p4_compStruct p4_comp;
typedef struct p4_rMatrixStruct p4_rMatrix;
typedef struct p4_gdasrvStruct p4_gdasrv;
typedef struct p4_pInvarStruct p4_pInvar;
typedef struct p4_mixtureStruct p4_mixture;
typedef struct p4_tSCovarionStruct p4_tSCovarion;
typedef struct p4_bigQAndEigStruct p4_bigQAndEig;

typedef struct nexusTokenStruct nexusToken;


struct dataStruct {
	int 	nTax;
	int     nParts;
	part  **parts;
	double  unconstrainedLogLike;
};

struct partStruct {
	data   *data;
	int 	dim;
	int     nTax;
	int		nChar;
	char   *symbols;
	int	  **patterns;
	int    *patternCounts;
	int    *sequencePositionPatternIndex;
	int		nPatterns;
	int	  **sequences;
	int     nEquates;
	int   **equates;
	char   *equateSymbols;
	int    *globalInvarSitesVec;
	int   **globalInvarSitesArray;
	double *siteLikes;
	int    *taxList;
	int    *simCats;
	//double  logLike;
};

struct eigStruct {
	int	dim;	// the dimension of the matrix, ie dim x dim
	double	**mat;	// the matrix
	double	**mwork;	// a utility matrix
	double	**eigvecs;	// the eigenvectors
	double	**inverseEigvecs;	// the inverse of the eigenvector matrix
	double	*eigvals;	// the eigenvalues
	double	*eigvalsImag;	// the imaginary part of the eigvals
	double	*dwork;		// a utility vector
	int		*iwork;		// a utility vector
	int      complexEig;
	pcomplex	**Ceigvecs;
	pcomplex	**CinverseEigvecs;
	pcomplex	**Cwork;
	pcomplex	*Ccol;
	int      complexStuffHasBeenAllocated;
};


struct pcomplexStruct{
	double re;
	double im;
};


struct rellStuffStruct{
	int             nTrees;
	int             nChar;
	double        **mat;
	const gsl_rng  *gsl_rng;
};



//===========================================
//===========================================
//===========================================


struct p4_treeStruct {
	int        nNodes;
    int        nLeaves;
	p4_node  **nodes;
	p4_node   *root;
    data      *data;
    p4_model  *model;
	int        nParts;
	int       *preOrder;
	int       *postOrder;
	int       *ints;  // utility vector, nNodes long.
	p4_node  **stack;
    double     logLike;
	double    *partLikes;
	int     ***simSequences;
	int     ***internalSequences;
};


struct p4_nodeStruct {
	int		       nodeNum;
	p4_tree       *tree;
	p4_node       *parent;
	p4_node       *leftChild;
	p4_node       *sibling;
	int            seqNum;
	int		       isLeaf;
	int            nParts;
	double        *brLen;
	double         savedBrLen;
	int           *compNums;
	int           *rMatrixNums;
	int           *gdasrvNums;
	double     ****bigPDecks;
	double     ****bigPDecks_1stD;
	double     ****bigPDecks_2ndD;
	double     ****cl;  // conditionalLikelihoods[partNum][nCats][dim][nChar]
	double     ****cl2;
	double     ****pickerDecks;
	int		       clNeedsUpdating;
	int            cl2NeedsUpdating;
	//int            brLenChanged;
	double	 ***expectedComp;  // need one for each part, one for each rate category
	double	  **observedComp;  // need one for each part
#if 0
	int         bFlag;
	int         uFlag;
	int       **catVec;         // for simulations.  replaced by using part->simRateCats
#endif
};


struct p4_modelStruct {
	int            nParts;
	p4_modelPart **parts;
	int            doRelRates;
	int            relRatesAreFree;
	int            nFreePrams;
	int            isHet;
	int           *rMatrixNormalizeTo1;
};

struct p4_modelPartStruct {
	int                dim;
  //int                clNeedsUpdating;
	int                nComps;
	p4_comp          **comps;
	int                nRMatrices;
	p4_rMatrix       **rMatrices;
	int                nGdasrvs;
	p4_gdasrv        **gdasrvs;
	int                nCat;
	p4_pInvar         *pInvar;
	int                isMixture;
	p4_mixture        *mixture;
	double            *relRate;
	p4_bigQAndEig   ***bigQAndEigThing; // nComps * nRMatrices
	int                doTSCovarion;
	p4_tSCovarion     *tSCov;
	double            *freqsTimesOneMinusPInvar;
	int               *bQETneedsReset;//a numpy int vector, nComps * nRMatrices, accessed by [(cNum * nRMatrices) + rNum]
};

struct p4_compStruct {
	//int        dim;
	int        free;
	double    *val;
};

struct p4_rMatrixStruct {
	//int        dim;
	int        free;
	int        spec;
	double   **bigR;
    double    *kappa;
};

struct p4_gdasrvStruct {
	int        free;
	double    *val;
	double    *freqs;
	double    *rates;
	int        nCat;
};

struct p4_pInvarStruct {
	int        free;
	double    *val;
};

struct p4_mixtureStruct {
	int        free;
	double    *freqs; // length nCat
	double    *rates;
};



struct p4_tSCovarionStruct {
	int        free;
	double    *s1; // rate on->off
	double    *s2; // rate off->on
    //int        halfDim;
    double    *halfComp;
	double     pOn;
	double     pOff;
};

struct p4_bigQAndEigStruct {
	double    **bigQ;
	eig        *qEig;
	//int         *needsReset;
};

// ================================

//struct np_gdasrvStruct {    // not used
//	int        free;
//	double    *val;
//	double    *freqs;
//	double    *rates;
//};

struct nexusTokenStruct {
	int       *writeVisibleComments;
	int       *getP4CommandComments;
	int       *getWeightCommandComments;
	int       *getAllCommandComments;
	int       *getLineEndings;
	int       *max;
	int       *tokLen;
	char      *tok;
	int       *embeddedCommentLen;
	char      *embeddedComment;
	int       *savedCommentLen;
	FILE      *filePtr;
};

