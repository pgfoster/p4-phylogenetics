#include <Python.h>
#include <numpy/arrayobject.h>

#include "pftypes.h"
#include "data.h"
#include "part.h"
#include "defines.h"

#include "dlsgamma.h"
#include "gamma.h"
#include "util.h"
#include "pmatrices.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics_double.h>
//#include <gsl/gsl_combination.h>
#include "eig.h"

#include "p4_tree.h"
#include "p4_node.h"
#include "p4_model.h"
#include "p4_treeCopyVerify.h"
#include "logDet.h"
#include "proteinModels.h"


static PyObject *
pf_newData(PyObject *self, PyObject *args)
{
    int nTax;
    int nParts;
	
    if(!PyArg_ParseTuple(args, "ii", &nTax, &nParts)) {
        printf("Error pf_newData: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("l", (long int)newData(nTax, nParts));
	
}

static PyObject *
pf_freeData(PyObject *self, PyObject *args)
{
    data	*theData;
	
    if(!PyArg_ParseTuple(args, "l", &theData)) {
        printf("Error pf_freeData: couldn't parse tuple\n");
        return NULL;
    }

    freeData(theData);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_dumpData(PyObject *self, PyObject *args)
{
    data	*theData;
	
    if(!PyArg_ParseTuple(args, "l", &theData)) {
        printf("Error pf_dumpData: couldn't parse tuple\n");
        return NULL;
    }

    dumpData(theData);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_pokePartInData(PyObject *self, PyObject *args)
{
    part    *thePart;
    data	*theData;
    int     i;
	
    if(!PyArg_ParseTuple(args, "lli", &thePart, &theData, &i)) {
        printf("Error pf_pokePartInData: couldn't parse tuple\n");
        return NULL;
    }

    theData->parts[i] = thePart;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_bootstrapData(PyObject *self, PyObject *args)
{
    data    *referenceData;
    data	*toFillData;
    const gsl_rng  *gsl_rng;
	
    if(!PyArg_ParseTuple(args, "lll", &referenceData, &toFillData, &gsl_rng)) {
        printf("Error pf_bootstrapData: couldn't parse tuple\n");
        return NULL;
    }

    bootstrapData(referenceData, toFillData, gsl_rng);
    Py_INCREF(Py_None);
    return Py_None;
}




//#####################################################################

static PyObject *
pf_newPart(PyObject *self, PyObject *args)
{
    int nTax;
    int nChar;
    int nEquates;
    char *symbols;
    char *equateSymbols;
    int	dim;
	
    if(!PyArg_ParseTuple(args, "iisisi", &nTax, &nChar, &equateSymbols, 
                         &nEquates, &symbols, &dim)) {
        printf("Error pf_newPart: couldn't parse tuple\n");
        return NULL;
    }

    //printf("got nTax %i, nChar %i, symbols %s, dim %i\n", nTax, nChar, symbols, dim);
    return Py_BuildValue("l", (long int)newPart(nTax, nChar, equateSymbols, nEquates, symbols, dim));
	
}


static PyObject *
pf_freePart(PyObject *self, PyObject *args)
{
    part	*thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_freePart: couldn't parse tuple\n");
        return NULL;
    }

    freePart(thePart);
    Py_INCREF(Py_None);
    return Py_None;
}



static PyObject *
pf_pokeEquatesTable(PyObject *self, PyObject *args)
{

    long int	thePart;
    char	*theString;
	
    if(!PyArg_ParseTuple(args, "ls", &thePart, &theString)) {
        printf("Error pf_pokeEquatesTable: couldn't parse tuple\n");
        return NULL;
    }

    pokeEquatesTable((part *)thePart, theString);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_pokeSequences(PyObject *self, PyObject *args)
{

    long int	thePart;
    char	*theString;
	
    if(!PyArg_ParseTuple(args, "ls", &thePart, &theString)) {
        printf("Error pf_pokeSequences: couldn't parse tuple\n");
        return NULL;
    }

    pokeSequences((part *)thePart, theString);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_makePatterns(PyObject *self, PyObject *args)
{

    long int	thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_makePatterns: couldn't parse tuple\n");
        return NULL;
    }

    makePatterns((part *)thePart);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_dumpPart(PyObject *self, PyObject *args)
{

    long int	thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_dumpPart: couldn't parse tuple\n");
        return NULL;
    }

    dumpPart((part *)thePart);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_singleSequenceBaseCounts(PyObject *self, PyObject *args)
{
    long int thePart;
    int sequenceNum;
	
    if(!PyArg_ParseTuple(args, "li", &thePart, &sequenceNum)) {
        printf("Error pf_singleSequenceBaseCounts: couldn't parse tuple\n");
        return NULL;
    }

    return singleSequenceBaseCounts((part *)thePart, sequenceNum);  // its a PyObject
}

static PyObject *
pf_symbolSequences(PyObject *self, PyObject *args)
{
    part *thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_symbolSequences: couldn't parse tuple\n");
        return NULL;
    }

    return symbolSequences(thePart);  // its a PyObject
	
}

static PyObject *
pf_partPatternCount(PyObject *self, PyObject *args)
{
    long int thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_partPatternCount: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("i", ((part *)thePart)->nPatterns);
	
}

static PyObject *
pf_partMeanNCharsPerSite(PyObject *self, PyObject *args)
{
    long int thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_partMeanNCharsPerSite: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", partMeanNCharsPerSite((part *)thePart));
	
}

static PyObject *
pf_partSimpleConstantSitesCount(PyObject *self, PyObject *args)
{
    long int thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_partSimpleConstantSitesCount: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("i", partSimpleConstantSitesCount((part *)thePart));
	
}
PyDoc_STRVAR(doc_partSimpleConstantSitesCount,
             "pf.partSimpleConstantSitesCount(part *p) -> int\n\
Of the sites that are not all gaps+ambigs,\n\
if the chars that are not gaps or ambigs are all the same \n\
then it is considered here to be a constant site.\n\
This function returns the count of such sites in the part.");


static PyObject *
pf_partBigXSquared(PyObject *self, PyObject *args)
{
    long int thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_partBigXSquared: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", partBigXSquared((part *)thePart));
	
}

static PyObject *
pf_setGlobalInvarSitesVec(PyObject *self, PyObject *args)
{
    long int	thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_setGlobalInvarSitesVec: couldn't parse tuple\n");
        return NULL;
    }

    setGlobalInvarSitesVec((part *)thePart);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_pokePartTaxListAtIndex(PyObject *self, PyObject *args)
{
    part	*thePart;
    int      val;
    int      indx;
	
    if(!PyArg_ParseTuple(args, "lii", &thePart, &val, &indx)) {
        printf("Error pf_pokePartTaxListAtIndex: couldn't parse tuple\n");
        return NULL;
    }

    if(thePart->taxList) {
        thePart->taxList[indx] = val;
    }
    else {
        printf("pf_pokePartTaxListAtIndex: no taxList\n");
        exit(1);
    }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_partComposition(PyObject *self, PyObject *args)
{
    part	*thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_partComposition: couldn't parse tuple\n");
        return NULL;
    }

    return partCompositionP(thePart);
}

static PyObject *
pf_partSequenceSitesCount(PyObject *self, PyObject *args)   // ie length of sequence, without gaps or missings.
{
    part	*thePart;
    int      sequenceNum;
	
    if(!PyArg_ParseTuple(args, "li", &thePart, &sequenceNum)) {
        printf("Error pf_partSequenceSitesCount: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("i", partSequenceSitesCount(thePart, sequenceNum));
}


static PyObject *
pf_getSiteLikes(PyObject *self, PyObject *args)
{
    part	 *thePart;
    PyObject *thePyList;
    int       i;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_getSiteLikes: couldn't parse tuple\n");
        return NULL;
    }

    thePyList = PyList_New(thePart->nChar);
    for(i = 0; i < thePart->nChar; i++) {
        PyList_SetItem(thePyList, i, PyFloat_FromDouble(thePart->siteLikes[i]));
    }
    return thePyList;
}

/*
static PyObject *
pf_getSiteRates(PyObject *self, PyObject *args)
{
    p4_tree  *theTree;
    part	 *thePart;
    int       pNum;
    PyArrayObject *oSiteRates;
    double       *siteRates;
    PyArrayObject *oGammaCats;
    int       *gammaCats;
    PyArrayObject *oWork;
    double    *work;
    //int       i;
    double    theLogLike;
	
    if(!PyArg_ParseTuple(args, "lliOOO", &theTree, &thePart, &pNum, &oSiteRates, &oGammaCats, &oWork)) {
        printf("Error pf_getSiteRates: couldn't parse tuple\n");
        return NULL;
    }
    //printf("pfmodule.c pf_getSiteRates here.\n");
    siteRates = (double *)(oSiteRates->data);
    gammaCats = (int *)(oGammaCats->data);
    work = (double *)(oWork->data);
    //for(i = 0; i < thePart->nChar; i++) {
    //	printf("%i   %i\n", i, gammaCats[i]);
    //}
    theLogLike = p4_partLogLikeSiteRates(theTree, thePart, pNum, 0, siteRates, gammaCats, work);
    //printf("part %i  logLike = %f\n", pNum, theLogLike);
    Py_INCREF(Py_None);
    return Py_None;
	
}
*/

static PyObject *
pf_getUnconstrainedLogLike(PyObject *self, PyObject *args)
{
    part	 *thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_getUnconstrainedLogLike: couldn't parse tuple\n");
        return NULL;
    }

    return PyFloat_FromDouble(unconstrainedLogLike(thePart));
}


static PyObject *
pf_calcEmpiricalRMatrixViaMatrixLog(PyObject *self, PyObject *args)
{
    part	 *thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_getUnconstrainedLogLike: couldn't parse tuple\n");
        return NULL;
    }

    calcEmpiricalRMatrixViaMatrixLog(thePart);
    Py_INCREF(Py_None);
    return Py_None;
	
}




// ######################################################



// #######################################################


static PyObject *
pf_reseedCRandomizer(PyObject *self, PyObject *args)
{
    int theSeed;
	
    if(!PyArg_ParseTuple(args, "i", &theSeed)) {
        printf("Error pf_reseedCRandomizer: couldn't parse tuple\n");
        return NULL;
    }
    srandom(theSeed);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_chiSquaredProb(PyObject *self, PyObject *args)
{
    double chiSq;
    int	df;
	
    if(!PyArg_ParseTuple(args, "di", &chiSq, &df)) {
        printf("Error pf_chiSquaredProb: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("d", Chi2P(chiSq, df));
}

static PyObject *
pf_studentsTProb(PyObject *self, PyObject *args)
{
    double t;
    int	df;
	
    if(!PyArg_ParseTuple(args, "di", &t, &df)) {
        printf("Error pf_studentsTProb: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("d", StudentTP(t, df));
}

static PyObject *
pf_normProb(PyObject *self, PyObject *args)
{
    double z; // Z = no. of standard deviations from the mean
	
    if(!PyArg_ParseTuple(args, "d", &z)) {
        printf("Error pf_normProb: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("d", NormP(z, NULL));
}

static PyObject *
pf_normPDF(PyObject *self, PyObject *args)
{
    double z; // Z = no. of standard deviations from the mean
    double pdf;
    //double nrmP;
    //PyObject	*thePyList;
	
    if(!PyArg_ParseTuple(args, "d", &z)) {
        printf("Error pf_normPDF: couldn't parse tuple\n");
        return NULL;
    }
    //thePyList = PyList_New(2);
    //nrmP = NormP(z, &pdf);
    //PyList_SetItem(thePyList, 0, PyFloat_FromDouble(nrmP));
    //PyList_SetItem(thePyList, 1, PyFloat_FromDouble(pdf));
    //return thePyList;
    NormP(z, &pdf);
    return Py_BuildValue("d", pdf);
}

static PyObject *
pf_pointChi2(PyObject *self, PyObject *args)
{
    double prob;
    double v; // degrees of freedom. Oddly, a double.
	
    if(!PyArg_ParseTuple(args, "dd", &prob, &v)) {
        printf("Error pf_pointChi2: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("d", PointChi2(prob, v));
}


// -----------------------------------
// rellStuff
// -----------------------------------


// static PyObject *
// pf_setRellMemory(PyObject *self, PyObject *args)
// {
//     int   nTrees;
//     int   nChar;
//     int   gslrng_seed;
//     rellStuff *rStuff;
//     const gsl_rng_type * T;
	
//     if(!PyArg_ParseTuple(args, "iii", &nTrees, &nChar, &gslrng_seed)) {
//         printf("Error pf_setRellMemory: couldn't parse tuple\n");
//         return NULL;
//     }
//     rStuff = malloc(sizeof(rellStuff));
//     if(!rStuff) {
//         printf("Failed to make rellStuff\n");
//         exit(1);
//     }
//     rStuff->nTrees = nTrees;
//     rStuff->nChar = nChar;
//     rStuff->mat = pdmatrix(nTrees, nChar);
//     gsl_rng_env_setup();  // make a generator.  See gsl docs.
//     T = gsl_rng_default;
//     rStuff->gsl_rng = gsl_rng_alloc(T);
//     if(!rStuff->gsl_rng) {
//         printf("setRellMemory.  Failed to make a gsl_rng.\n");
//         exit(1);
//     }
//     // set the seed for the random number generator
//     gsl_rng_set(rStuff->gsl_rng, (unsigned long int)gslrng_seed);

//     return Py_BuildValue("l", (long int)rStuff);
// }

// static PyObject *
// pf_pokeRellMemory(PyObject *self, PyObject *args)
// {
//     int        treeNum;
//     int        charNum;
//     double     val;
//     rellStuff *rStuff;
	
//     if(!PyArg_ParseTuple(args, "iidl", &treeNum, &charNum, &val, &rStuff)) {
//         printf("Error pf_pokeRellMemory: couldn't parse tuple\n");
//         return NULL;
//     }
//     rStuff->mat[treeNum][charNum] = val;
//     Py_INCREF(Py_None);
//     return Py_None;
// }


// static PyObject *
// pf_rell(PyObject *self, PyObject *args)
// {
//     int   nBoots;
//     rellStuff *rStuff;
	
//     if(!PyArg_ParseTuple(args, "il", &nBoots, &rStuff)) {
//         printf("Error pf_rell: couldn't parse tuple\n");
//         return NULL;
//     }
//     return rell(nBoots, rStuff);
// }


// static PyObject *
// pf_freeRellMemory(PyObject *self, PyObject *args)
// {
//     rellStuff *rStuff;
	
//     if(!PyArg_ParseTuple(args, "l", &rStuff)) {
//         printf("Error pf_freeRellMemory: couldn't parse tuple\n");
//         return NULL;
//     }
//     if(rStuff->gsl_rng) {
//         gsl_rng_free((gsl_rng *)rStuff->gsl_rng);
//         rStuff->gsl_rng = NULL;
//     }
//     free_pdmatrix(rStuff->mat);
//     free(rStuff);
//     Py_INCREF(Py_None);
//     return Py_None;
// }


// -----------------------------------
// gsl stuff
// -----------------------------------


// static PyObject *
// pf_get_gsl_rng(PyObject *self, PyObject *args)
// {
//     const gsl_rng_type * T;
//     const gsl_rng  *gsl_rng;
	
	
//     gsl_rng_env_setup();  // make a generator.  See gsl docs.
//     T = gsl_rng_default;
//     gsl_rng = gsl_rng_alloc(T);
//     if(!gsl_rng) {
//         printf("pf_get_gsl_rng.  Failed to make a gsl_rng.\n");
//         exit(1);
//     }
//     return Py_BuildValue("l", (long int)gsl_rng);
// }


static PyObject *
pf_gsl_rng_get(PyObject *self, PyObject *args)
{
    const gsl_rng_type * T;
    const gsl_rng  *g;
	
	
    gsl_rng_env_setup();  // make a generator.  See gsl docs.
    T = gsl_rng_default;
    g = gsl_rng_alloc(T);
    if(!g) {
        printf("pf_gsl_rng_get.  Failed to make a gsl_rng.\n");
        exit(1);
    }
    return Py_BuildValue("l", (long int)g);
}

static PyObject *
pf_gsl_rng_free(PyObject *self, PyObject *args)
{
    gsl_rng  *g;   // not const
	
    if(!PyArg_ParseTuple(args, "l", &g)) {
        printf("Error pf_gsl_rng_free: couldn't parse tuple\n");
        return NULL;
    }

    gsl_rng_free(g);
    Py_INCREF(Py_None);
    return Py_None;
	
}


static PyObject *
pf_gsl_rng_set(PyObject *self, PyObject *args)
{
    const gsl_rng  *g;
    long int        newSeed;
	
    if(!PyArg_ParseTuple(args, "ll", &g, &newSeed)) {
        printf("Error pf_gsl_rng_set: couldn't parse tuple\n");
        return NULL;
    }

    gsl_rng_set(g, (unsigned long int)newSeed);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_gsl_rng_size(PyObject *self, PyObject *args)
{
    const gsl_rng  *g;
    size_t          theSize;
	
    if(!PyArg_ParseTuple(args, "l", &g)) {
        printf("Error pf_gsl_rng_size: couldn't parse tuple\n");
        return NULL;
    }

    theSize = gsl_rng_size(g);
    return Py_BuildValue("i", (int)theSize);
}

static PyObject *
pf_gsl_rng_uniform(PyObject *self, PyObject *args)
{
    const gsl_rng  *g;	
	
    if(!PyArg_ParseTuple(args, "l", &g)) {
        printf("Error pf_gsl_rng_uniform: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", (double)gsl_rng_uniform(g));
}

static PyObject *
pf_gsl_rng_getstate(PyObject *self, PyObject *args)
{
    const gsl_rng * r;	
    PyArrayObject * oByteArray;
    void * state;
    size_t n;
    
    if(!PyArg_ParseTuple(args, "lO", &r, &oByteArray)) {
        printf("Error pf_gsl_rng_getstate: couldn't parse tuple\n");
        return NULL;
    }

    // The C library function 
    // void *memcpy(void *str1, const void *str2, size_t n) 
    // copies n characters from memory area str2 to
    // memory area str1.

    state = gsl_rng_state(r);
    n = gsl_rng_size(r);
    memcpy(oByteArray->data, state, n);

    Py_INCREF(Py_None);
    return Py_None;
    
}

static PyObject *
pf_gsl_rng_setstate(PyObject *self, PyObject *args)
{
    const gsl_rng * r;	
    PyArrayObject * oByteArray;
    void * state;
    size_t n;
	
    if(!PyArg_ParseTuple(args, "lO", &r, &oByteArray)) {
        printf("Error pf_gsl_rng_setstate: couldn't parse tuple\n");
        return NULL;
    }

    // The C library function 
    // void *memcpy(void *str1, const void *str2, size_t n) 
    // copies n characters from memory area str2 to
    // memory area str1.

    state = gsl_rng_state(r);
    n = gsl_rng_size(r);
    memcpy(state, oByteArray->data, n);

    Py_INCREF(Py_None);
    return Py_None;
    
}


static PyObject *
pf_gsl_ran_chisq_pdf(PyObject *self, PyObject *args)
{
    double x, nu;
	
    if(!PyArg_ParseTuple(args, "dd", &x, &nu)) {
        printf("Error pf_gsl_ran_chisq_pdf: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", gsl_ran_chisq_pdf(x, nu));
	
}
PyDoc_STRVAR(doc_gsl_ran_chisq_pdf,
             "pf.gsl_ran_chisq_pdf(double x, double dof) -> double\n\
Wrapper on the GSL function.\n\
Return the probability density function at a point x for a chi square\n\
with dof degrees of freedom.");


static PyObject *
pf_gsl_ran_gamma(PyObject *self, PyObject *args)
{
    double a, b;
    const gsl_rng  *R;
	
    if(!PyArg_ParseTuple(args, "ldd", &R, &a, &b)) {
        printf("Error pf_gsl_ran_gamma: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", gsl_ran_gamma(R, a, b));
	
}
PyDoc_STRVAR(doc_gsl_ran_gamma,
             "pf.gsl_ran_gamma(const gsl_rng *R, double a, double b) -> double\n\
Wrapper on the GSL function.\n\
Return a number drawn from the gamma distribution.\n\
Only use this with func.gsl_ran_gamma()");


static PyObject *
pf_gsl_ran_gamma_pdf(PyObject *self, PyObject *args)
{
    double x, a, b;
	
    if(!PyArg_ParseTuple(args, "ddd", &x, &a, &b)) {
        printf("Error pf_gsl_ran_gamma_pdf: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", gsl_ran_gamma_pdf(x, a, b));
	
}
PyDoc_STRVAR(doc_gsl_ran_gamma_pdf,
             "pf.gsl_ran_gamma_pdf(double x, double a, double b) -> double\n\
Wrapper on the GSL function.\n\
This function computes the probability density p(x) at x for a gamma distribution with parameters a and b.");

//static PyObject *
//pf_gsl_ran_exponential_pdf(PyObject *self, PyObject *args)
//{
//	double x, mu;
//	
//	if(!PyArg_ParseTuple(args, "dd", &x, &mu)) {
//		printf("Error pf_gsl_ran_exponential_pdf: couldn't parse tuple\n");
//		return NULL;
//	}
//
//	return Py_BuildValue("d", gsl_ran_exponential_pdf(x, mu));
//	
//}
/*
  PyDoc_STRVAR(doc_gsl_ran_exponential_pdf,
  "pf.gsl_ran_exponential_pdf(double x, double mu) -> double\n	\
  Wrapper on the GSL function.\n\
  This is not the same as the Exp() dist used in branch length priors!\n\
  This function computes the probability density p(x) at x for an exponential distribution with mean mu");
*/

static PyObject *
pf_gsl_sf_lngamma(PyObject *self, PyObject *args)
{
    double a;
	
    if(!PyArg_ParseTuple(args, "d", &a)) {
        printf("Error pf_gsl_sf_lngamma: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", gsl_sf_lngamma(a));
	
}
PyDoc_STRVAR(doc_gsl_sf_lngamma,
             "pf.gsl_sf_lngamma(double a) -> double\n\
Wrapper on the GSL function.\n");

static PyObject *
pf_gsl_sf_gamma(PyObject *self, PyObject *args)
{
    double a;
	
    if(!PyArg_ParseTuple(args, "d", &a)) {
        printf("Error pf_gsl_sf_gamma: couldn't parse tuple\n");
        return NULL;
    }

    return Py_BuildValue("d", gsl_sf_gamma(a));
	
}
PyDoc_STRVAR(doc_gsl_sf_gamma,
             "pf.gsl_sf_gamma(double a) -> double\n\
Wrapper on the GSL function.\n");


static PyObject *
pf_gsl_sf_beta(PyObject *self, PyObject *args)
{
    double a, b;
  
    if(!PyArg_ParseTuple(args, "dd", &a, &b)) {
        printf("Error pf_gsl_sf_beta: couldn't parse tuple\n");
        return NULL;
    }
  
    return Py_BuildValue("d", gsl_sf_beta(a, b));
	
}
PyDoc_STRVAR(doc_gsl_sf_beta,
             "pf.gsl_sf_beta(double a, double b) -> double\n\
Wrapper on the GSL function.\n");


static PyObject *
pf_gsl_ran_dirichlet(PyObject *self, PyObject *args)
{
    PyArrayObject *oAlpha, *oTheta;
    double *alpha, *theta;
    int k;
    const gsl_rng  *R;
	
    if(!PyArg_ParseTuple(args, "liOO", &R, &k, &oAlpha, &oTheta)) {
        printf("Error pf_gsl_ran_dirichlet: couldn't parse tuple\n");
        return NULL;
    }
    alpha = (double *)(oAlpha->data);
    theta = (double *)(oTheta->data);

    // theta gets overwritten with the new new numbers.
    gsl_ran_dirichlet(R, k, alpha, theta);
	
    Py_INCREF(Py_None);
    return Py_None;
	
}
PyDoc_STRVAR(doc_gsl_ran_dirichlet,
             "pf.gsl_ran_dirichlet(const gsl_rng *R, int k, Numpy_array alpha, Numpy_array theta) -> void\n\
Wrapper on the GSL function.\n\
Return a draw from the dirichlet distribution.\n\
Only use this with func.gsl_ran_dirichlet()");



static PyObject *
pf_gsl_ran_dirichlet_pdf(PyObject *self, PyObject *args)
{
    PyArrayObject *oAlpha, *oTheta;
    double *alpha, *theta;
    int k;
	
    if(!PyArg_ParseTuple(args, "iOO", &k, &oAlpha, &oTheta)) {
        printf("Error pf_gsl_dirichlet_pdf: couldn't parse tuple\n");
        return NULL;
    }
    alpha = (double *)(oAlpha->data);
    theta = (double *)(oTheta->data);
    return Py_BuildValue("d", gsl_ran_dirichlet_pdf(k, alpha, theta));
}
PyDoc_STRVAR(doc_gsl_ran_dirichlet_pdf,
             "pf.gsl_ran_dirichlet_pdf(int k, NumPy vector of floats alpha, NumPy vector of floats theta) -> double\n\n");


static PyObject *
pf_gsl_ran_dirichlet_lnpdf(PyObject *self, PyObject *args)
{
    PyArrayObject *oAlpha, *oTheta;
    double *alpha, *theta;
    int k;
	
    if(!PyArg_ParseTuple(args, "iOO", &k, &oAlpha, &oTheta)) {
        printf("Error pf_gsl_dirichlet_lnpdf: couldn't parse tuple\n");
        return NULL;
    }
    alpha = (double *)(oAlpha->data);
    theta = (double *)(oTheta->data);
    return Py_BuildValue("d", gsl_ran_dirichlet_lnpdf(k, alpha, theta));
}
PyDoc_STRVAR(doc_gsl_ran_dirichlet_lnpdf,
             "pf.gsl_ran_dirichlet_lnpdf(int k, NumPy vector of floats alpha, NumPy vector of floats theta) -> double\n\n");


static PyObject *
pf_gsl_ran_lognormal_pdf(PyObject *self, PyObject *args)
{
    double x, zeta, sigma;
	
    if(!PyArg_ParseTuple(args, "ddd", &x, &zeta, &sigma)) {
        printf("Error pf_gsl_lognormal_pdf: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("d", gsl_ran_lognormal_pdf(x, zeta, sigma));
}
PyDoc_STRVAR(doc_gsl_ran_lognormal_pdf,
             "pf.gsl_ran_lognormal_pdf(float x, float zeta, float sigma) -> double\n\n");




static PyObject *
pf_gsl_meanVariance(PyObject *self, PyObject *args)
{
    PyArrayObject *oSeq, *oMean, *oVariance;
    double *seq, *mean, *variance;
    int seqLen;
    //int i;
	
    if(!PyArg_ParseTuple(args, "OiOO", &oSeq, &seqLen, &oMean, &oVariance)) {
        printf("Error pf_gsl_meanVariance: couldn't parse tuple\n");
        return NULL;
    }
    seq = (double *)(oSeq->data);
    mean = (double *)(oMean->data);
    variance = (double *)(oVariance->data);
    //printf("mean->data is at %li\n", (long int)mean);
    //printf("variance->data is at %li\n", (long int)variance);
    //for(i = 0; i < seqLen; i++) {
    //	printf("%f ", seq[i]);
    //}
    //printf("\n");
    mean[0] = gsl_stats_mean(seq, 1, seqLen);
    //printf("xxx variance is %f\n", gsl_stats_variance_m(seq, 1, seqLen, mean[0]));
    variance[0] = gsl_stats_variance_m(seq, 1, seqLen, mean[0]);
    //printf("%f %f\n", mean[0], variance[0]);

    Py_INCREF(Py_None);
    return Py_None;
	
}
PyDoc_STRVAR(doc_gsl_meanVariance,
             "pf.gsl_meanVariance(NumPy vector of floats, N, NumPy mean, NumPy variance)\n\n");


#if 0
static PyObject *
pf_gsl_combination_alloc(PyObject *self, PyObject *args)
{
    gsl_combination  *c;
    int n, k;
	
    if(!PyArg_ParseTuple(args, "ii", &n, &k)) {
        printf("Error pf_gsl_combination_alloc: couldn't parse tuple\n");
        return NULL;
    }
    c = gsl_combination_alloc(n, k);
    if(!c) {
        printf("pf_gsl_combination_alloc.  Failed to alloc.\n");
        exit(1);
    }
    return Py_BuildValue("l", (long int)c);
}

static PyObject *
pf_gsl_combination_free(PyObject *self, PyObject *args)
{
    gsl_combination  *c;

    if(!PyArg_ParseTuple(args, "l", &c)) {
        printf("Error pf_gsl_combination_free: couldn't parse tuple\n");
        return NULL;
    }
    gsl_combination_free(c);
    Py_INCREF(Py_None);
    return Py_None;
}
#endif






// -----------------------------------
// Eig stuff
// -----------------------------------


static PyObject *
pf_getSquareDoubleMatrix(PyObject *self, PyObject *args)
{
    int dim;
	
    if(!PyArg_ParseTuple(args, "i", &dim)) {
        printf("Error pf_getSquareDoubleMatrix: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("l", psdmatrix(dim));
}

static PyObject *
pf_pokeValInSquareDoubleMatrix(PyObject *self, PyObject *args)
{
    int i, j;
    double val;
    long int  mat;
	
    if(!PyArg_ParseTuple(args, "ldii", &mat, &val, &i, &j)) {
        printf("Error pf_pokeValInSquareDoubleMatrix: couldn't parse tuple\n");
        return NULL;
    }
    ((double **)mat)[i][j] = val;

    Py_INCREF(Py_None);
    return Py_None;
	
}

static PyObject *
pf_getValFromSquareDoubleMatrix(PyObject *self, PyObject *args)
{
    int i, j;
    long int  mat;
	
    if(!PyArg_ParseTuple(args, "lii", &mat, &i, &j)) {
        printf("Error pf_getValFromSquareDoubleMatrix: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("d", ((double **)mat)[i][j]);
}

static PyObject *
pf_freeSquareDoubleMatrix(PyObject *self, PyObject *args)
{
    long int  mat;
	
    if(!PyArg_ParseTuple(args, "l", &mat)) {
        printf("Error pf_freeSquareDoubleMatrix: couldn't parse tuple\n");
        return NULL;
    }
    free_psdmatrix((double **)mat);

    Py_INCREF(Py_None);
    return Py_None;
	
}

static PyObject *
pf_newEig(PyObject *self, PyObject *args)
{
    int dim;
    long int mat;
	
    if(!PyArg_ParseTuple(args, "il", &dim, &mat)) {
        printf("Error pf_newEig: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("l", allocEig(dim, (double **)mat));
}

static PyObject *
pf_freeEig(PyObject *self, PyObject *args)
{
    long int theEig;
	
    if(!PyArg_ParseTuple(args, "l", &theEig)) {
        printf("Error pf_freeEig: couldn't parse tuple\n");
        return NULL;
    }
    freeEig((eig *)theEig);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_matrixExp(PyObject *self, PyObject *args)
{
    long int theEig;
    long int result;
    double factor;
    int i;
	
    if(!PyArg_ParseTuple(args, "ldl", &theEig, &factor, &result)) {
        printf("Error pf_matrixExp: couldn't parse tuple\n");
        return NULL;
    }

    //printf("result = %li\n", result);

    // void matrixExpTimesBranchLength(eig *anEig, double branchLength, double **result)
    for(i = 0; i < 100000; i++) {
	eigensystem((eig *)theEig);
        matrixExpTimesBranchLength((eig *)theEig, factor, (double **)result);
    }
    //dump_psdmatrix(((eig *)theEig)->mat, 4);
    //printf("factor = %f\n", factor);
    //dump_psdmatrix((double **)result, 4);

    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_matrixLog(PyObject *self, PyObject *args)
{
    long int theEig;
    long int result;
	
    if(!PyArg_ParseTuple(args, "ll", &theEig, &result)) {
        printf("Error pf_matrixExp: couldn't parse tuple\n");
        return NULL;
    }

    eigensystem((eig *)theEig);
    // void matrixLog(eig *anEig, double **result)
    matrixLog((eig *)theEig, (double **)result);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_getBigQ(PyObject *self, PyObject *args)
{
    p4_model *aModel;
    int theDim;
    int partNum;
    int compNum;
    int rMatrixNum;
    PyArrayObject *oNumpyBigQ;
    int i, j;
    double **eigBigQ;
    double *numpyData;
	
    if(!PyArg_ParseTuple(args, "liiiiO", &aModel, &theDim, &partNum, &compNum, &rMatrixNum, &oNumpyBigQ)) {
        printf("Error pf_getBigQ: couldn't parse tuple\n");
        return NULL;
    }

    //printf("getBigQ here.\n");
    eigBigQ = aModel->parts[partNum]->bigQAndEigThing[compNum][rMatrixNum]->bigQ;
    numpyData = (double *)(oNumpyBigQ->data);
    //dump_psdmatrix(eigBigQ, theDim);
    for(i = 0; i < theDim; i++) {
        for(j = 0; j < theDim; j++) {
	    numpyData[(i * theDim) + j] = eigBigQ[i][j];
        }
    } 
	

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_getBigR(PyObject *self, PyObject *args)
{
    int spec;
    PyArrayObject *oNumpyBigR;
    int i, j;
    double **bigR;
    double *numpyData;
    int theDim = 20;
	
    if(!PyArg_ParseTuple(args, "iO", &spec, &oNumpyBigR)) {
        printf("Error pf_getBigR: couldn't parse tuple\n");
        return NULL;
    }

    bigR = (double **)psdmatrix(theDim);
    // fill the bigR
    if(spec == RMATRIX_CPREV) {
        cpREVRMatrix(bigR);
    }
    else if(spec == RMATRIX_D78) {
        d78RMatrix(bigR);
    }
    else if(spec == RMATRIX_JTT) {
        jttRMatrix(bigR);
    }
    else if (spec == RMATRIX_MTREV24) {
        mtREV24RMatrix(bigR);
    }
    else if(spec == RMATRIX_MTMAM) {
        mtmamRMatrix(bigR);
    }
    else if(spec == RMATRIX_WAG) {
        wagRMatrix(bigR);
    }
    else if(spec == RMATRIX_RTREV) {
        rtRevRMatrix(bigR);
    }

    else if(spec == RMATRIX_TMJTT94) {
        tmjtt94RMatrix(bigR);
    }
    else if(spec == RMATRIX_TMLG99) {
        tmlg99RMatrix(bigR);
    }
    else if(spec == RMATRIX_LG) {
        lgRMatrix(bigR);
    }
    else if(spec == RMATRIX_BLOSUM62_A) {
        blosum62RMatrix(bigR);
    }
    else if(spec == RMATRIX_HIVB) {
        hivbRMatrix(bigR);
    }
    else if(spec == RMATRIX_MTART) {
        mtartRMatrix(bigR);
    }
    else if(spec == RMATRIX_MTZOA) {
        mtzoaRMatrix(bigR);
    }
    else if(spec == RMATRIX_GCPREV) {
        gcpREVRMatrix(bigR);
    }
    else if(spec == RMATRIX_STMTREV) {
        stmtREVRMatrix(bigR);
    }
    else if(spec == RMATRIX_VT) {
        vtRMatrix(bigR);
    }
    else if(spec == RMATRIX_PMB) {
        pmbRMatrix(bigR);
    }

    numpyData = (double *)(oNumpyBigR->data);
    //dump_psdmatrix(eigBigQ, theDim);
    for(i = 0; i < theDim; i++) {
        for(j = 0; j < theDim; j++) {
	    numpyData[(i * theDim) + j] = bigR[i][j];
        }
    } 
    free_psdmatrix(bigR);
    bigR = NULL;

    Py_INCREF(Py_None);
    return Py_None;
}

// -----------------------------------
// 
// -----------------------------------


static PyObject *
pf_recodeNLike(PyObject *self, PyObject *args)
{
    part *thePart;
	
    if(!PyArg_ParseTuple(args, "l", &thePart)) {
        printf("Error pf_recodeNLike: couldn't parse tuple\n");
        return NULL;
    }
    recodeNLike(thePart);
    Py_INCREF(Py_None);
    return Py_None;
}



// -----------------------------------
// ======== p4_functions
// -----------------------------------

static PyObject *
pf_p4_newTree(PyObject *self, PyObject *args)
{
    int nNodes;
    int nLeaves;
    data *aData;
    p4_model *aModel;
    PyArrayObject *oPreOrder, *oPostOrder, *oPartLikes, *oNewtAndBrentPowellOptPassLimit;
    int  *preOrder, *postOrder, *newtAndBrentPowellOptPassLimit;
    double *partLikes;
    //int i;
	
    if(!PyArg_ParseTuple(args, "iiOOOOll",&nNodes, &nLeaves, &oPreOrder, &oPostOrder, 
                         &oNewtAndBrentPowellOptPassLimit, &oPartLikes, &aData, &aModel)) {
        printf("Error pf_p4_newTree: couldn't parse tuple\n");
        return NULL;
    }

    preOrder = (int *)oPreOrder->data;
    postOrder = (int *)(oPostOrder->data);
    partLikes = (double *)(oPartLikes->data);
    newtAndBrentPowellOptPassLimit = (int *)(oNewtAndBrentPowellOptPassLimit->data);
    //printf("pfmodule.c pf_p4_newTree. newtAndBrentPowellOptPassLimit is %i\n", newtAndBrentPowellOptPassLimit[0]);
    return Py_BuildValue("l", p4_newTree(nNodes, nLeaves, preOrder, postOrder, partLikes, aData, aModel, newtAndBrentPowellOptPassLimit));
}


static PyObject *
pf_p4_freeTree(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;
	
    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_freeTree: couldn't parse tuple\n");
        return NULL;
    }
    p4_freeTree(aTree);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_dumpTree(PyObject *self, PyObject *args)
{
    p4_tree   *aTree;
	
    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_dumpTree: couldn't parse tuple\n");
        return NULL;
    }
    p4_dumpTree(aTree);
    Py_INCREF(Py_None);
    return Py_None;
}



// -----------------------------------

static PyObject *
pf_p4_newNode(PyObject *self, PyObject *args)
{
    int nodeNum;
    p4_tree   *aTree;
    int        seqNum;
    int        isLeaf;
    int        inTree;
	
    if(!PyArg_ParseTuple(args, "iliii", &nodeNum, &aTree, &seqNum, &isLeaf, &inTree)) {
        printf("Error pf_p4_newNode: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("l", p4_newNode(nodeNum, aTree, seqNum, isLeaf, inTree));
}


static PyObject *
pf_p4_freeNode(PyObject *self, PyObject *args)
{
    p4_node  *aNode;
	
    if(!PyArg_ParseTuple(args, "l", &aNode)) {
        printf("Error pf_p4_freeNode: couldn't parse tuple\n");
        return NULL;
    }

    if(!aNode->tree->model) {
        printf("p4_freeNode().  A model is required, and this has no model\n");
        exit(1);
    }
    p4_freeNode(aNode);
    Py_INCREF(Py_None);
    return Py_None;
}




// -----------------------------------

static PyObject *
pf_p4_newModel(PyObject *self, PyObject *args)
{
    int nParts;
    int doRelRates;
    int relRatesAreFree;
    int nFreePrams;
    int isHet;
    PyArrayObject *oRMatrixNormalizeTo1;
    PyArrayObject *oPINVAR_MIN;
    PyArrayObject *oPINVAR_MAX;
    PyArrayObject *oKAPPA_MIN;
    PyArrayObject *oKAPPA_MAX;
    PyArrayObject *oGAMMA_SHAPE_MIN;
    PyArrayObject *oGAMMA_SHAPE_MAX;
    PyArrayObject *oPIVEC_MIN;
    PyArrayObject *oPIVEC_MAX;
    PyArrayObject *oRATE_MIN;
    PyArrayObject *oRATE_MAX;
    PyArrayObject *oRELRATE_MIN;
    PyArrayObject *oRELRATE_MAX;
    PyArrayObject *oBRLEN_MIN;
    PyArrayObject *oBRLEN_MAX;
	
    if(!PyArg_ParseTuple(args, "iiiiiOOOOOOOOOOOOOOO", &nParts, 
                         &doRelRates, 
                         &relRatesAreFree, 
                         &nFreePrams, 
                         &isHet, 
                         &oRMatrixNormalizeTo1,
                         &oPINVAR_MIN,
                         &oPINVAR_MAX,
                         &oKAPPA_MIN,
                         &oKAPPA_MAX,
                         &oGAMMA_SHAPE_MIN,
                         &oGAMMA_SHAPE_MAX,
                         &oPIVEC_MIN,
                         &oPIVEC_MAX,
                         &oRATE_MIN,
                         &oRATE_MAX,
                         &oRELRATE_MIN,
                         &oRELRATE_MAX,
                         &oBRLEN_MIN,
                         &oBRLEN_MAX )) {
        printf("Error pf_p4_newModel: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("l", p4_newModel(nParts, 
                                          doRelRates, 
                                          relRatesAreFree, 
                                          nFreePrams, 
                                          isHet, 
                                          (int *)(oRMatrixNormalizeTo1->data),
                                          (double *)(oPINVAR_MIN->data),
                                          (double *)(oPINVAR_MAX->data),
                                          (double *)(oKAPPA_MIN->data),
                                          (double *)(oKAPPA_MAX->data),
                                          (double *)(oGAMMA_SHAPE_MIN->data),
                                          (double *)(oGAMMA_SHAPE_MAX->data),
                                          (double *)(oPIVEC_MIN->data),
                                          (double *)(oPIVEC_MAX->data),
                                          (double *)(oRATE_MIN->data),
                                          (double *)(oRATE_MAX->data),
                                          (double *)(oRELRATE_MIN->data),
                                          (double *)(oRELRATE_MAX->data),
                                          (double *)(oBRLEN_MIN->data),
                                          (double *)(oBRLEN_MAX->data)));
}


static PyObject *
pf_p4_freeModel(PyObject *self, PyObject *args)
{
    p4_model  *aModel;
	
    if(!PyArg_ParseTuple(args, "l", &aModel)) {
        printf("Error pf_p4_freeModel: couldn't parse tuple\n");
        return NULL;
    }
    if(aModel) {
        p4_freeModel(aModel);
    }
    else {
        printf("p4_freeModel().  aModel is NULL.\n");
    }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_dumpModel(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
	
    if(!PyArg_ParseTuple(args, "l", &aModel)) {
        printf("Error pf_p4_dumpModel: couldn't parse tuple\n");
        return NULL;
    }
    p4_dumpModel(aModel);
    Py_INCREF(Py_None);
    return Py_None;
}

// -----------------------------------

static PyObject *
pf_p4_newModelPart(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         dim;
    int         nComps;
    int         nRMatrices;
    int         nGdasrvs;
    int         nCat;
    int         pInvarFree;
    PyArrayObject *bQETneedsResetO;
	

    if(!PyArg_ParseTuple(args, "liiiiiiiO", &aModel, &pNum, &dim, &nComps, &nRMatrices, &nGdasrvs, &nCat, &pInvarFree, &bQETneedsResetO)) {
        printf("Error pf_p4_newModelPart: couldn't parse tuple\n");
        return NULL;
    }
    p4_newModelPart(aModel, pNum, dim, nComps, nRMatrices, nGdasrvs, nCat, pInvarFree, (int *)(bQETneedsResetO->data));
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_resetBQET(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         compNum;
    int         rMatrixNum;
	

    if(!PyArg_ParseTuple(args, "liii", &aModel, &pNum, &compNum, &rMatrixNum)) {
        printf("Error pf_p4_resetBQET: couldn't parse tuple\n");
        return NULL;
    }
    p4_resetBQET(aModel, pNum, compNum, rMatrixNum);
    Py_INCREF(Py_None);
    return Py_None;
}




static PyObject *
pf_p4_newComp(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         mNum;
    int         free;
    PyArrayObject  *val;

    if(!PyArg_ParseTuple(args, "liiiO", &aModel, &pNum, &mNum, &free, &val)) {
        printf("Error pf_p4_newComp: couldn't parse tuple\n");
        return NULL;
    }
    p4_newComp(aModel, pNum, mNum, free, val);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_newRMatrix(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         mNum;
    int         free;
    int         spec;

    if(!PyArg_ParseTuple(args, "liiii", &aModel, &pNum, &mNum, &free, &spec)) {
        printf("Error pf_p4_newRMatrix: couldn't parse tuple\n");
        return NULL;
    }
    p4_newRMatrix(aModel, pNum, mNum, free, spec);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_newGdasrv(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         mNum;
    int         nCat;
    int         free;
    PyArrayObject  *val, *freqs, *rates;

    if(!PyArg_ParseTuple(args, "liiiiOOO", &aModel, &pNum, &mNum, &nCat, &free, &val, &freqs, &rates)) {
        printf("Error pf_p4_newGdasrv: couldn't parse tuple\n");
        return NULL;
    }
    //printf("p4_newGdasrv. nCat=%i\n", nCat);
    return Py_BuildValue("l", (long int)p4_newGdasrv(aModel, pNum, mNum, nCat, free, val, freqs, rates));
}





//---------------------------------------

#if 0
// no longer needed as comp.val is a numpy array
static PyObject *
pf_p4_setCompVal(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         cNum;
    int         vNum;
    double      val;
	
    if(!PyArg_ParseTuple(args, "liiid", &aModel, &pNum, &cNum, &vNum, &val)) {
        printf("Error pf_p4_setCompVal: couldn't parse tuple\n");
        return NULL;
    }

    aModel->parts[pNum]->comps[cNum]->val[vNum] = val;

    Py_INCREF(Py_None);
    return Py_None;
}
#endif


static PyObject *
pf_p4_setRMatrixBigR(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         cNum;
    int         i;
    int         j;
    double      val;
	
    if(!PyArg_ParseTuple(args, "liiiid", &aModel, &pNum, &cNum, &i, &j, &val)) {
        printf("Error pf_p4_setRMatrixBigR: couldn't parse tuple\n");
        return NULL;
    }
    aModel->parts[pNum]->rMatrices[cNum]->bigR[i][j] = val;
    aModel->parts[pNum]->rMatrices[cNum]->bigR[j][i] = val;
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_setKappa(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         cNum;
    double      val;
    double      alpha = 1.0/3.0;
    double      beta;
    p4_rMatrix *r;
	
    if(!PyArg_ParseTuple(args, "liid", &aModel, &pNum, &cNum, &val)) {
        printf("Error pf_p4_setKappa: couldn't parse tuple\n");
        return NULL;
    }
    r = aModel->parts[pNum]->rMatrices[cNum];
    r->kappa[0] = val;
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

    Py_INCREF(Py_None);
    return Py_None;
}


#if 0
static PyObject *
pf_p4_setGdasrvVal(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    int         cNum;
    double      val;
    p4_gdasrv  *g;
	
    if(!PyArg_ParseTuple(args, "liid", &aModel, &pNum, &cNum, &val)) {
        printf("Error pf_p4_setGdasrvVal: couldn't parse tuple\n");
        return NULL;
    }
    g = aModel->parts[pNum]->gdasrvs[cNum];
    g->val[0] = val;

    // set the rates
    DiscreteGamma(g->freqs, g->rates, g->val[0], g->val[0], aModel->parts[pNum]->nCat, 0);
    if(1){
        int i;
        printf("pfmodule: setGdasrvVal: gdasrv[part %i] shape=%f, free=%i\n", pNum, g->val[0], g->free);
        for(i = 0; i < aModel->parts[pNum]->nCat; i++) {
            printf("  %i  %f   %f\n", i, g->freqs[i], g->rates[i]);
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}
#endif


static PyObject *
pf_gdasrvCalcRates(PyObject *self, PyObject *args)
{
    p4_gdasrv  *g;
	
    if(!PyArg_ParseTuple(args, "l", &g)) {
        printf("Error pf_gdasrvCalcRates: couldn't parse tuple\n");
        return NULL;
    }

    //printf("about to DiscreteGamma. g->val[0]=%f\n", g->val[0]);
    DiscreteGamma(g->freqs, g->rates, g->val[0], g->val[0], g->nCat, 0);
    if((0)){
        int i;
        printf("pfmodule: pf_gdasrvCalcRates: shape=%f, free=%i\n", g->val[0], g->free);
        for(i = 0; i < g->nCat; i++) {
            printf("  %i  %f   %f\n", i, g->freqs[i], g->rates[i]);
        }
    }
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_gdasrvCalcRates_np(PyObject *self, PyObject *args)
{
    int  nGammaCat;
    double val;
    PyArrayObject *oFreqs;
    double *freqs;
    PyArrayObject *oRates;
    double *rates;
	
	
    if(!PyArg_ParseTuple(args, "idOO", &nGammaCat, &val, &oFreqs, &oRates)) {
        printf("Error pf_gdasrvCalcRates_np: couldn't parse tuple\n");
        return NULL;
    }
    //printf("pfmodule: pf_gdasrvCalcRates_np() here.\n");
    freqs = (double *)oFreqs->data;
    rates = (double *)oRates->data;
    DiscreteGamma(freqs, rates, val, val, nGammaCat, 0);
#if 0
    {
        int i;
        printf("pfmodule: pf_gdasrvCalcRates_np: shape=%f, free=? nGammaCat=%i\n", val, nGammaCat);
        for(i = 0; i < nGammaCat; i++) {
            printf("  %i  %f   %f\n", i, freqs[i], rates[i]);
        }
    }
#endif
    Py_INCREF(Py_None);
    return Py_None;
}




static PyObject *
pf_p4_setPInvarVal(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    double      val;
	
    if(!PyArg_ParseTuple(args, "lid", &aModel, &pNum, &val)) {
        printf("Error pf_p4_setPInvarVal: couldn't parse tuple\n");
        return NULL;
    }
    aModel->parts[pNum]->pInvar->val[0] = val;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_setRelRateVal(PyObject *self, PyObject *args)
{
    p4_model   *aModel;
    int         pNum;
    double      val;
	
    if(!PyArg_ParseTuple(args, "lid", &aModel, &pNum, &val)) {
        printf("Error pf_p4_setRelRateVal: couldn't parse tuple\n");
        return NULL;
    }
    aModel->parts[pNum]->relRate[0] = val;
    Py_INCREF(Py_None);
    return Py_None;
}


//--------------------------------------------------


static PyObject *
pf_p4_setNodeRelation(PyObject *self, PyObject *args)
{
    p4_node *aNode;
    int relation;
    int relNum;
    p4_node *theRelative = NULL;
	
    if(!PyArg_ParseTuple(args, "lii", &aNode, &relation, &relNum)) {
        printf("Error pf_p4_setNodeRelation: couldn't parse tuple\n");
        return NULL;
    }

    if(relNum >= 0) {
        theRelative = aNode->tree->nodes[relNum];
    } else {
        theRelative = NULL;
    }
    if(relation == 0) {
        aNode->parent = theRelative;
    }
    else if(relation == 1) {
        aNode->leftChild = theRelative;
    }
    else if(relation == 2) {
        aNode->sibling = theRelative;
    }
    else {
        printf("Error in p4_setNodeRelation: \"relation\" is out of range\n");
        exit(1);
    }
	
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_setTreeRoot(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;
    p4_node  *aNode;

    if(!PyArg_ParseTuple(args, "ll", &aTree, &aNode)) {
        printf("Error pf_p4_setTreeRoot: couldn't parse tuple\n");
        return NULL;
    }
    aTree->root = aNode;
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_setBrLen(PyObject *self, PyObject *args)
{
    p4_node  *aNode;
    double    theBrLen;
    //int       brLenChanged;

    //if(!PyArg_ParseTuple(args, "ldi", &aNode, &theBrLen, &brLenChanged)) {
    if(!PyArg_ParseTuple(args, "ld", &aNode, &theBrLen)) {
        printf("Error pf_p4_setBrLen: couldn't parse tuple\n");
        return NULL;
    }
    aNode->brLen[0] = theBrLen;
    //aNode->brLenChanged = brLenChanged;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_getTreeLen(PyObject *self, PyObject *args)
{
    p4_tree *aTree;
    p4_node *aNode;
    double leng;
    int i, j;
  
    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_getTreeLen: couldn't parse tuple\n");
        return NULL;
    }
  
    leng = 0.0;
    for(i = 0; i < aTree->nNodes; i++) {
        j = aTree->postOrder[i];
        if(j != NO_ORDER) {
            aNode = aTree->nodes[j];
            if(aNode != aTree->root) {
                leng += aNode->brLen[0];
            }
        }
    }
  
    return Py_BuildValue("d", leng);
	
}
PyDoc_STRVAR(doc_p4_getTreeLen,
             "pf.p4_getTreeLen(cTree) -> double\n\
Get sum of br lens.  If the node is in theTree.nodes but not in the tree, it is skipped.\n");



//-----------------------

static PyObject *
pf_p4_setCompNum(PyObject *self, PyObject *args)
{
    p4_node  *aNode;
    int       pNum;
    int       val;

    if(!PyArg_ParseTuple(args, "lii", &aNode, &pNum, &val)) {
        printf("Error pf_p4_setCompNum: couldn't parse tuple\n");
        return NULL;
    }
    aNode->compNums[pNum] = val;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_setRMatrixNum(PyObject *self, PyObject *args)
{
    p4_node  *aNode;
    int       pNum;
    int       val;

    if(!PyArg_ParseTuple(args, "lii", &aNode, &pNum, &val)) {
        printf("Error pf_p4_setRMatrixNum: couldn't parse tuple\n");
        return NULL;
    }
    aNode->rMatrixNums[pNum] = val;
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_setGdasrvNum(PyObject *self, PyObject *args)
{
    p4_node  *aNode;
    int       pNum;
    int       val;

    if(!PyArg_ParseTuple(args, "lii", &aNode, &pNum, &val)) {
        printf("Error pf_p4_setGdasrvNum: couldn't parse tuple\n");
        return NULL;
    }
    aNode->gdasrvNums[pNum] = val;
    Py_INCREF(Py_None);
    return Py_None;
}


//static PyObject *
//pf_p4_setPreAndPostOrder(PyObject *self, PyObject *args)
//{
//	p4_tree  *aTree;
//	int       pos;
//	int       pre;
//	int       post;
//
//	if(!PyArg_ParseTuple(args, "liii", &aTree, &pos, &pre, &post)) {
//		printf("Error pf_p4_setPreAndPostOrder: couldn't parse tuple\n");
//		return NULL;
//	}
//	aTree->preOrder[pos] = pre;
//	aTree->postOrder[pos] = post;
//    Py_INCREF(Py_None);
//	return Py_None;
//}



//------------------------------



static PyObject *
pf_p4_setPrams(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;
    int       pNum;

    if(!PyArg_ParseTuple(args, "li", &aTree, &pNum)) {
        printf("Error pf_p4_setPrams: couldn't parse tuple\n");
        return NULL;
    }
    if(pNum == -1) {
        p4_setPrams(aTree);
    }
    else {
        p4_setPramsPart(aTree, pNum);
        //aTree->model->parts[pNum]->clNeedsUpdating = 1;
    }
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_setPramsTest(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;
    int       pNum;

    if(!PyArg_ParseTuple(args, "li", &aTree, &pNum)) {
        printf("Error pf_p4_setPramsTest: couldn't parse tuple\n");
        return NULL;
    }
    if(pNum == -1) {
        p4_setPramsTest(aTree);
    }
    else {
        p4_setPramsPartTest(aTree, pNum);
        //aTree->model->parts[pNum]->clNeedsUpdating = 1;
    }
    
    Py_INCREF(Py_None);
    return Py_None;
}



static PyObject *
pf_p4_treeLogLike(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;
    int       getSiteLikes;

    if(!PyArg_ParseTuple(args, "li", &aTree, &getSiteLikes)) {
        printf("Error pf_p4_treeLogLike: couldn't parse tuple\n");
        return NULL;
    }

    //printf("wxyz parts=%li\n", (long int)aTree->model->parts); 
    return Py_BuildValue("d", p4_treeLogLike(aTree, getSiteLikes));
}

static PyObject *
pf_p4_partLogLike(PyObject *self, PyObject *args)
{
    p4_tree *treePtr;
    part *part;
    int   partNum;
    int   doSiteLikes;
	
    if(!PyArg_ParseTuple(args, "llii", &treePtr, &part, &partNum, &doSiteLikes)) {
        printf("Error pf_p4_partLogLike: couldn't parse tuple\n");
        return NULL;
    }
    // double p4_partLogLike(p4_tree *aTree, part *dp, int pNum, int getSiteLikes)
    return Py_BuildValue("d", p4_partLogLike(treePtr, part, partNum, doSiteLikes));
}



static PyObject *
pf_p4_calculateBigPDecks(PyObject *self, PyObject *args)
{
    p4_node  *aNode;

    if(!PyArg_ParseTuple(args, "l", &aNode)) {
        printf("Error pf_p4_calculateBigPDecks: couldn't parse tuple\n");
        return NULL;
    }

    p4_calculateBigPDecks(aNode);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_calculateAllBigPDecksAllParts(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;

    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_calculateAllBigPDecksAllParts: couldn't parse tuple\n");
        return NULL;
    }

    p4_calculateAllBigPDecksAllParts(aTree);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_setConditionalLikelihoodsOfInternalNodePart(PyObject *self, PyObject *args)
{
    p4_node  *aNode;
    int       pNum;

    if(!PyArg_ParseTuple(args, "li", &aNode, &pNum)) {
        printf("Error pf_p4_setConditionalLikelihoodsOfInternalNodePart: couldn't parse tuple\n");
        return NULL;
    }

    p4_setConditionalLikelihoodsOfInternalNodePart(aNode, pNum);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_allBOBYQAOptimize(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;
    int doBrLens;

    if(!PyArg_ParseTuple(args, "li", &aTree, &doBrLens)) {
        printf("Error pf_p4_allBOBYQAOptimize: couldn't parse tuple\n");
        return NULL;
    }
    p4_allBOBYQAOptimize(aTree, doBrLens);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_newtAndBOBYQAOpt(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;

    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_newtAndBOBYQAOpt: couldn't parse tuple\n");
        return NULL;
    }
    p4_newtAndBOBYQAOpt(aTree);
    Py_INCREF(Py_None);
    return Py_None;
}




static PyObject *
pf_p4_allBrentPowellOptimize(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;

    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_allBrentPowellOptimize: couldn't parse tuple\n");
        return NULL;
    }
    p4_allBrentPowellOptimize(aTree);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_newtAndBrentPowellOpt(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;

    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_newtAndBrentPowellOpt: couldn't parse tuple\n");
        return NULL;
    }
    p4_newtAndBrentPowellOpt(aTree);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_newtSetup(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;

    if(!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_p4_newtSetup: couldn't parse tuple\n");
        return NULL;
    }
    p4_newtSetup(aTree);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_getBrLens(PyObject *self, PyObject *args)
{
    p4_tree *theTree;
 
    if(!PyArg_ParseTuple(args, "l", &theTree)) {
        printf("Error pf_p4_getBrLens: couldn't parse tuple\n");
        return NULL;
    }
 
    return p4_getBrLens(theTree);  // its a PyObject! -- a list!
}


static PyObject *
pf_p4_getFreePrams(PyObject *self, PyObject *args)
{
    p4_tree *theTree;
 
    if(!PyArg_ParseTuple(args, "l", &theTree)) {
        printf("Error pf_p4_getFreePrams: couldn't parse tuple\n");
        return NULL;
    }
 
    return p4_getFreePrams(theTree);  // its a PyObject! -- a list!
}

static PyObject *
pf_p4_getRelRate(PyObject *self, PyObject *args)
{
    p4_model *theModel;
    int       partNum;
 
    if(!PyArg_ParseTuple(args, "li", &theModel, &partNum)) {
        printf("Error pf_p4_getRelRate: couldn't parse tuple\n");
        return NULL;
    }
 
    return Py_BuildValue("d", theModel->parts[partNum]->relRate[0]);
	
}

//----------------------------------------

static PyObject *
pf_p4_simulate(PyObject *self, PyObject *args)
{
    p4_tree *aTree;
    p4_tree *refTree;
    const gsl_rng  *g;
	
    if(!PyArg_ParseTuple(args, "lll", &aTree, &refTree, &g)) {
        printf("Error pf_simulate: couldn't parse tuple\n");
        return NULL;
    }
    // printf("pf_p4_simulate(), got refTree %li\n", (long int)refTree);
    // if(refTree == NULL) {
    //     printf("pf_p4_simulate() refTree is NULL\n");
    // }
    p4_simulate(aTree, refTree, g);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_p4_drawAncState(PyObject *self, PyObject *args)
{
    p4_tree *aTree;
    int partNum;
    int seqPos;
    PyArrayObject *oResultsTuple;
	
    if(!PyArg_ParseTuple(args, "liiO", &aTree, &partNum, &seqPos, &oResultsTuple)) {
        printf("Error pf_simulate: couldn't parse tuple\n");
        return NULL;
    }
    // Note this calls p4_drawAncStateP() not p4_drawAncState()
    p4_drawAncStateP(aTree, partNum, seqPos, (int *)oResultsTuple->data);
    Py_INCREF(Py_None);
    return Py_None;
}
PyDoc_STRVAR(doc_p4_drawAncState,
             "pf.p4_drawAncState(p4_tree *t, int partNum, int seqPos, numpy.ndarray draw) -> None\n\
Arg t is a eg myPythonTree.cTree\n\
Arg seqPos is the sequence position in the part.\n\
Arg draw is a numpy array of length 4, dtype numpy.int32\n\
eg myDraw = numpy.zeros(4, numpy.int32)\n\
\n\
Make a draw from the tree ancestral state (at the root).  Results are in draw.\n\
draw[0] = chStNum;      # char state number if variable\n\
draw[1] = catNum;       # gamma category if variable\n\
draw[2] = isInvar;\n\
draw[3] = invarChNum;   # char state num if invariable");


static PyObject *
pf_p4_expectedCompositionCounts(PyObject *self, PyObject *args)
{
    p4_tree *theTree;
    int      partNum;
 
    if(!PyArg_ParseTuple(args, "li", &theTree, &partNum)) {
        printf("Error pf_p4_expectedCompositionCounts: couldn't parse tuple\n");
        return NULL;
    }
 
    return p4_expectedCompositionCounts(theTree, partNum);  // in p4_treeSim.c
}
PyDoc_STRVAR(doc_p4_expectedCompositionCounts,
             "pf.p4_expectedCompositionCounts(p4_tree *t, int partNum) -> Py tuple\n\
Returns a tuple of tuples.  Each inner tuple is a composition.\n\
Terminal nodes only.  The comps are in the order of the data.\n\
Ie the first is from the node that has seqNum=0.");


static PyObject *
pf_p4_expectedComposition(PyObject *self, PyObject *args)
{
    p4_tree *theTree;
    //int      partNum;
 
    if(!PyArg_ParseTuple(args, "l", &theTree)) {
        printf("Error pf_p4_expectedComposition: couldn't parse tuple\n");
        return NULL;
    }
 
    return p4_expectedComposition(theTree);  // in p4_treeSim.c
	
}
PyDoc_STRVAR(doc_p4_expectedComposition,
             "pf.p4_expectedComposition(p4_tree *t) -> Py tuple\n\
Returns a tuple of tuples of tuples.  Each inner tuple is a composition\n\
of a sequence.  Each middle tuple is a data partition.\n\
Terminal nodes only.  The comps are in the order of the data.\n\
Ie the first is from the node that has seqNum=0.");




//------------------------------


static PyObject *
pf_p4_verifyIdentityOfTwoTrees(PyObject *self, PyObject *args)
{
    p4_tree *aTree;
    p4_tree *bTree;
	
    if(!PyArg_ParseTuple(args, "ll", &aTree, &bTree)) {
        printf("Error pf_verifyIdentityOfTwoTrees: couldn't parse tuple\n");
        return NULL;
    }
    return Py_BuildValue("i", p4_verifyIdentityOfTwoTrees(aTree, bTree));
}


static PyObject *
pf_p4_copyCondLikes(PyObject *self, PyObject *args)
{
    p4_tree *aTree;
    p4_tree *bTree;
    int doAll;
	
    if(!PyArg_ParseTuple(args, "lli", &aTree, &bTree, &doAll)) {
        printf("Error pf_copyCondLikes: couldn't parse tuple\n");
        return NULL;
    }
    p4_copyCondLikes(aTree, bTree, doAll);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_copyBigPDecks(PyObject *self, PyObject *args)
{
    p4_tree *aTree;
    p4_tree *bTree;
    int doAll;
	
    if(!PyArg_ParseTuple(args, "lli", &aTree, &bTree, &doAll)) {
        printf("Error pf_copyBigPDecks: couldn't parse tuple\n");
        return NULL;
    }
    p4_copyBigPDecks(aTree, bTree, doAll);
    Py_INCREF(Py_None);
    return Py_None;
}


static PyObject *
pf_p4_copyModelPrams(PyObject *self, PyObject *args)
{
    p4_tree *aTree;
    p4_tree *bTree;
	
    if(!PyArg_ParseTuple(args, "ll", &aTree, &bTree)) {
        printf("Error pf_copyModelPrams: couldn't parse tuple\n");
        return NULL;
    }
    p4_copyModelPrams(aTree, bTree);
    Py_INCREF(Py_None);
    return Py_None;
}



static PyObject *
pf_zeroNumPyInts(PyObject *self, PyObject *args)
{
    PyArrayObject *theNumPyArray;
    int theLen;
    int  *theData;
    int   i;
	
    if(!PyArg_ParseTuple(args, "Oi", &theNumPyArray, &theLen)) {
        printf("Error pf_zeroNumPyInts: couldn't parse tuple\n");
        return NULL;
    }
    theData = (int *)theNumPyArray->data;
    for(i = 0; i < theLen; i++) {
        theData[i] = 0;
    }
    Py_INCREF(Py_None);
    return Py_None;
}

// ------------------------------------------

static PyObject *
pf_logDetFillFxy(PyObject *self, PyObject *args)
{
    PyArrayObject *nUnambig,
        *nAmbig,
        *nDoubleGap,
        *seq,
        *refUnambigCountMatrix,
        *refAmbigCountMatrix,
        *allSymbolNums,
        *normUnambig,
        *bigFxy,
        *equatesArray;

    int sNum1, sNum2, nTax, nChar, dim, bigDim;
	
    if(!PyArg_ParseTuple(args, "OOOOiiOOOiiiiOOO", &nUnambig,
                         &nAmbig,
                         &nDoubleGap,
                         &seq,
                         &sNum1,
                         &sNum2,
                         &refUnambigCountMatrix,
                         &refAmbigCountMatrix,
                         &allSymbolNums,
                         &nTax,
                         &nChar,
                         &dim,
                         &bigDim,
                         &normUnambig,
                         &bigFxy,
                         &equatesArray)) {
        printf("Error pf_logDetFillFxy: couldn't parse tuple\n");
        return NULL;
    }

    logDetFillFxy((int *)(nUnambig->data), 
                  (int *)(nAmbig->data), 
                  (int *)(nDoubleGap->data), 
                  (char *)(seq->data), 
                  sNum1, 
                  sNum2, 
                  (int *)(refUnambigCountMatrix->data), 
                  (int *)(refAmbigCountMatrix->data), 
                  (int *)(allSymbolNums->data),
                  nTax,
                  nChar,
                  dim,
                  bigDim,
                  (double *)(normUnambig->data),
                  (double *)(bigFxy->data),
                  (int *)(equatesArray->data));
	
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_effectiveSampleSize(PyObject *self, PyObject *args)
{
    PyArrayObject *valsO,
        *meanO,
        *gammaStatAtPreviousLagO,
        *gammaStatO,
        *varStatO,
        *gammaStatAtLagZeroO;

    int nSamples, maxLag, j, lag;
    double *vals,
        *mean,
        *gammaStatAtPreviousLag,
        *gammaStat,
        *varStat,
        *gammaStatAtLagZero;
	
    if(!PyArg_ParseTuple(args, "OOiiOOOO", &valsO,
                         &meanO,
                         &nSamples,
                         &maxLag,
                         &gammaStatAtPreviousLagO,
                         &gammaStatO,
                         &varStatO,
                         &gammaStatAtLagZeroO)) {
        printf("Error pf_effectiveSampleSize: couldn't parse tuple\n");
        return NULL;
    }

    vals = (double *)(valsO->data);
    mean = (double *)(meanO->data);
    gammaStatAtPreviousLag = (double *)(gammaStatAtPreviousLagO->data);
    gammaStat = (double *)(gammaStatO->data);
    varStat = (double *)(varStatO->data);
    gammaStatAtLagZero = (double *)(gammaStatAtLagZeroO->data);
    //printf("c here 1, maxLag is %i, nSamples is %i, mean is %f\n", maxLag, nSamples, mean[0]);
    lag = 0;
    while(lag < maxLag) {
        gammaStat[0] = 0.0;
        for(j = 0; j < (nSamples - lag); j++) {
            gammaStat[0] += (vals[j] - mean[0]) * (vals[j + lag] - mean[0]);
        }

        gammaStat[0] /= (double)(nSamples - lag);

        if(lag == 0) {
            varStat[0] = gammaStat[0];
            gammaStatAtLagZero[0] = gammaStat[0];
            //printf("c got gammaStatAtLagZero = %f\n", gammaStatAtLagZero[0]);
        }
        else if((lag % 2) == 0) {
            if(gammaStatAtPreviousLag[0] + gammaStat[0] > 0.0) {
                varStat[0] += 2.0 * (gammaStatAtPreviousLag[0] + gammaStat[0]);
            }
            else {
                break;
            }
        }

        lag += 1;
        gammaStatAtPreviousLag[0] = gammaStat[0];
        //gammaStat[0] = 0.0;
    }

    //printf("c: lag is %i, gammaStatAtLagZero is %g, varStat is %f\n", lag, gammaStatAtLagZero[0], varStat[0]);
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
pf_newtonRaftery94_eqn16(PyObject *self, PyObject *args)
{
    PyArrayObject *logLikesO;
    double harmMean, 
        delta,
        *logLikes;
    int verbose, len;
	
    if(!PyArg_ParseTuple(args, "Oiddi", &logLikesO, &len, &harmMean, &delta, &verbose)) {
        printf("Error pf_newtonRaftery94_eqn16: couldn't parse tuple\n");
        return NULL;
    }

    logLikes = (double *)(logLikesO->data);

    //return Py_BuildValue("d", p4_partLogLike(treePtr, part, partNum, doSiteLikes));
    return Py_BuildValue("d", newtonRaftery94_eqn16(logLikes, len, harmMean, delta, verbose));
}

// ==================
// callback
// ==================

// static PyObject *mcmcTreeCallback = NULL;

static PyObject *
pf_setMcmcTreeCallback(PyObject *dummy, PyObject *args)
{
    p4_tree  *aTree;
    PyObject *result = NULL;
    PyObject *temp;

    if (!PyArg_ParseTuple(args, "lO:setMcmcTreeCallback", &aTree, &temp)) {
        printf("Error pf_setMcmcTreeCallback: couldn't parse tuple\n");
        return NULL;
    }
    if (!PyCallable_Check(temp)) {
            PyErr_SetString(PyExc_TypeError, "pf_setMcmcTreeCallback(): parameter must be callable");
            return NULL;
    }
    Py_XINCREF(temp);              /* Add a reference to new callback */
    Py_XDECREF(aTree->mcmcTreeCallback);  /* Dispose of previous callback */
    aTree->mcmcTreeCallback = temp;       /* Remember new callback */
    /* Boilerplate to return "None" */
    Py_INCREF(Py_None);
    result = Py_None;
    
    return result;
}

static PyObject *
pf_unsetMcmcTreeCallback(PyObject *dummy, PyObject *args)
{
    p4_tree  *aTree;

    if (!PyArg_ParseTuple(args, "l", &aTree)) {
        printf("Error pf_unsetMcmcTreeCallback: couldn't parse tuple\n");
        return NULL;
    }
    Py_XDECREF(aTree->mcmcTreeCallback);  /* Dispose of previous callback */
    aTree->mcmcTreeCallback = NULL; 
    Py_INCREF(Py_None);
    return Py_None;
}


#if 1
static PyObject *
pf_test(PyObject *self, PyObject *args)
{

    printf("pf_test here.\n");
	
    Py_INCREF(Py_None);
    return Py_None;
}
#endif

#if 0
static PyObject *
pf_test(PyObject *self, PyObject *args)
{
    PyArrayObject *oTheArray;
    int *theIntPtr;
    int theInt;
	
    if(!PyArg_ParseTuple(args, "O", &oTheArray)) {
        printf("Error pf_test: couldn't parse tuple\n");
        return NULL;
    }

    printf("xxxx %li\n", (long int)oTheArray->data);
    printf("xyxx %i\n", (int)oTheArray->data[0]);            // no workee
    printf("xyyx %i\n", (int)((int *)(oTheArray->data))[0]); // works
    theIntPtr = (int *)oTheArray->data;
    theInt = theIntPtr[0];
    printf("xzxx %i\n", theInt); 
    Py_INCREF(Py_None);
    return Py_None;
}
#endif

#if 0
static PyObject *
pf_test(PyObject *self, PyObject *args)
{
    p4_tree  *aTree;
    //p4_node  *aNode;
    int pNum, compNum, chNum;
    p4_modelPart   *mp;
    double myIt;


    if(!PyArg_ParseTuple(args, "liii", &aTree, &pNum, &compNum, &chNum)) {
        printf("Error pf_test: couldn't parse tuple\n");
        return NULL;
    }
    
    //p4_setPramsPart(aTree, pNum);
    //void p4_resetBigQAndEig(p4_tree *aTree, int pNum)
    //p4_resetGdasrv(aTree, pNum);
    //printf("\n");

    // {
    //     p4_modelPart   *mp;
    //     double sum=0.0;
    //     int i, j;

    //     mp = aTree->model->parts[pNum];
    //     //printf("    part %i, nComps=%i, nRMatrices=%i, nCat=%i\n", pNum, mp->nComps, mp->nRMatrices, mp->nCat);
    //     for(i = 0; i < mp->nComps; i++) {
    //         for(j = 0; j < mp->dim; j++) {
    //             //if(mp->comps[i]->val[j] < (0.5 * PIVEC_MIN)) {
    //             if(mp->comps[i]->val[j] < PIVEC_MIN) {               // Was half PIVEC_MIN.  Why?
    //                 printf("p4_setPramsPart()  part %i, comp %i, value %i is %g   Bad.\n", pNum, i, j, mp->comps[i]->val[j]);
    //                 exit(1);
    //             }
    //         }
    //         // check
    //         sum = 0.0;
    //         for(j = 0; j < mp->dim; j++) {
    //             sum += mp->comps[i]->val[j];
    //         }
    //         // if its wrong this time, give up  (diff of 1e-16 was too small, raised to 1e-14)
    //         if(fabs(sum - 1.0) > 1e-14) {
    //             printf("****p4_setPramsPart()  part %i, comp %i, values do not sum to 1.0.  sum=%g, %f\n", pNum, i, sum, sum);
    //             printf("****  sum - 1.0 = %g\n", sum - 1.0);
    //             //exit(1);
    //         } else {
    //             printf("****p4_setPramsPart()  part %i, comp %i, ok\n", pNum, i);
    //         }
    //     }
    // }

    mp = aTree->model->parts[pNum];
    myIt = mp->comps[compNum]->val[chNum];
    
    return Py_BuildValue("d", myIt);
    // Py_INCREF(Py_None);
    // return Py_None;
}
#endif





// ======================================================
// ######################################################
// ======================================================


static PyMethodDef pfMethods[] = {
    {"newData", pf_newData, METH_VARARGS},
    {"freeData", pf_freeData, METH_VARARGS},
    {"dumpData", pf_dumpData, METH_VARARGS},
    {"pokePartInData", pf_pokePartInData, METH_VARARGS},
    {"bootstrapData", pf_bootstrapData, METH_VARARGS},

    {"newPart", pf_newPart, METH_VARARGS},
    {"freePart", pf_freePart, METH_VARARGS},
    {"pokeSequences", pf_pokeSequences, METH_VARARGS},
    {"pokeEquatesTable", pf_pokeEquatesTable, METH_VARARGS},
    {"makePatterns", pf_makePatterns, METH_VARARGS},
    {"dumpPart", pf_dumpPart, METH_VARARGS},
    {"singleSequenceBaseCounts", pf_singleSequenceBaseCounts, METH_VARARGS},
    {"symbolSequences", pf_symbolSequences, METH_VARARGS},
    {"partPatternCount", pf_partPatternCount, METH_VARARGS},
    {"partMeanNCharsPerSite", pf_partMeanNCharsPerSite, METH_VARARGS},
    {"partSimpleConstantSitesCount", pf_partSimpleConstantSitesCount, METH_VARARGS, doc_partSimpleConstantSitesCount},
    {"partBigXSquared", pf_partBigXSquared, METH_VARARGS},
    {"setGlobalInvarSitesVec", pf_setGlobalInvarSitesVec, METH_VARARGS},
    {"pokePartTaxListAtIndex", pf_pokePartTaxListAtIndex, METH_VARARGS},
    {"partComposition", pf_partComposition, METH_VARARGS},
    {"partSequenceSitesCount", pf_partSequenceSitesCount, METH_VARARGS},
    {"getSiteLikes", pf_getSiteLikes, METH_VARARGS},
    //{"getSiteRates", pf_getSiteRates, METH_VARARGS},
    {"getUnconstrainedLogLike", pf_getUnconstrainedLogLike, METH_VARARGS},
    {"calcEmpiricalRMatrixViaMatrixLog", pf_calcEmpiricalRMatrixViaMatrixLog, METH_VARARGS},

    {"reseedCRandomizer", pf_reseedCRandomizer, METH_VARARGS},
    {"chiSquaredProb", pf_chiSquaredProb, METH_VARARGS},
    {"studentsTProb", pf_studentsTProb, METH_VARARGS},
    {"normProb", pf_normProb, METH_VARARGS},
    {"normPDF", pf_normPDF, METH_VARARGS},
    {"pointChi2", pf_pointChi2, METH_VARARGS},

    //{"setRellMemory", pf_setRellMemory, METH_VARARGS},
    //{"pokeRellMemory", pf_pokeRellMemory, METH_VARARGS},
    //{"rell", pf_rell, METH_VARARGS},
    //{"freeRellMemory", pf_freeRellMemory, METH_VARARGS},

    //{"get_gsl_rng", pf_get_gsl_rng, METH_VARARGS},
    {"gsl_rng_get", pf_gsl_rng_get, METH_VARARGS},
    {"gsl_rng_free", pf_gsl_rng_free, METH_VARARGS},
    {"gsl_rng_set", pf_gsl_rng_set, METH_VARARGS},
    {"gsl_rng_size", pf_gsl_rng_size, METH_VARARGS},
    {"gsl_rng_uniform", pf_gsl_rng_uniform, METH_VARARGS},
    {"gsl_rng_getstate", pf_gsl_rng_getstate, METH_VARARGS},
    {"gsl_rng_setstate", pf_gsl_rng_setstate, METH_VARARGS},
    {"gsl_ran_chisq_pdf", pf_gsl_ran_chisq_pdf, METH_VARARGS, doc_gsl_ran_chisq_pdf},
    {"gsl_ran_gamma", pf_gsl_ran_gamma, METH_VARARGS, doc_gsl_ran_gamma},
    {"gsl_ran_gamma_pdf", pf_gsl_ran_gamma_pdf, METH_VARARGS, doc_gsl_ran_gamma_pdf},
    //{"gsl_ran_exponential_pdf", pf_gsl_ran_exponential_pdf, METH_VARARGS, doc_gsl_ran_exponential_pdf},
    {"gsl_sf_lngamma", pf_gsl_sf_lngamma, METH_VARARGS, doc_gsl_sf_lngamma},
    {"gsl_sf_gamma", pf_gsl_sf_gamma, METH_VARARGS, doc_gsl_sf_gamma},
    {"gsl_sf_beta", pf_gsl_sf_beta, METH_VARARGS, doc_gsl_sf_beta},
    {"gsl_ran_dirichlet", pf_gsl_ran_dirichlet, METH_VARARGS, doc_gsl_ran_dirichlet},
    {"gsl_ran_dirichlet_pdf", pf_gsl_ran_dirichlet_pdf, METH_VARARGS, doc_gsl_ran_dirichlet_pdf},
    {"gsl_ran_dirichlet_lnpdf", pf_gsl_ran_dirichlet_lnpdf, METH_VARARGS, doc_gsl_ran_dirichlet_lnpdf},
    {"gsl_ran_lognormal_pdf", pf_gsl_ran_lognormal_pdf, METH_VARARGS, doc_gsl_ran_lognormal_pdf},
    {"gsl_meanVariance", pf_gsl_meanVariance, METH_VARARGS, doc_gsl_meanVariance},
    //{"gsl_combination_alloc", pf_gsl_combination_alloc, METH_VARARGS},
    //{"gsl_combination_free", pf_gsl_combination_free, METH_VARARGS},

    {"getSquareDoubleMatrix", pf_getSquareDoubleMatrix, METH_VARARGS},
    {"pokeValInSquareDoubleMatrix", pf_pokeValInSquareDoubleMatrix, METH_VARARGS},
    {"getValFromSquareDoubleMatrix", pf_getValFromSquareDoubleMatrix, METH_VARARGS},
    {"freeSquareDoubleMatrix", pf_freeSquareDoubleMatrix, METH_VARARGS},
    {"newEig", pf_newEig, METH_VARARGS},
    {"freeEig", pf_freeEig, METH_VARARGS},
    {"matrixExp", pf_matrixExp, METH_VARARGS},
    {"matrixLog", pf_matrixLog, METH_VARARGS},
    {"getBigQ", pf_getBigQ, METH_VARARGS},
    {"getBigR", pf_getBigR, METH_VARARGS},

    {"recodeNLike", pf_recodeNLike, METH_VARARGS},

    {"p4_newTree", pf_p4_newTree, METH_VARARGS},
    {"p4_freeTree", pf_p4_freeTree, METH_VARARGS},
    {"p4_dumpTree", pf_p4_dumpTree, METH_VARARGS},

    {"p4_newNode", pf_p4_newNode, METH_VARARGS},
    {"p4_freeNode", pf_p4_freeNode, METH_VARARGS},

    {"p4_newModel", pf_p4_newModel, METH_VARARGS},
    {"p4_freeModel", pf_p4_freeModel, METH_VARARGS},
    {"p4_dumpModel", pf_p4_dumpModel, METH_VARARGS},

    {"p4_newModelPart", pf_p4_newModelPart, METH_VARARGS},
    {"p4_resetBQET", pf_p4_resetBQET, METH_VARARGS},
    {"p4_newComp", pf_p4_newComp, METH_VARARGS},
    {"p4_newRMatrix", pf_p4_newRMatrix, METH_VARARGS},
    {"p4_newGdasrv", pf_p4_newGdasrv, METH_VARARGS},

    //{"p4_setCompVal", pf_p4_setCompVal, METH_VARARGS},
    {"p4_setRMatrixBigR", pf_p4_setRMatrixBigR, METH_VARARGS},
    {"p4_setKappa", pf_p4_setKappa, METH_VARARGS},
    //{"p4_setGdasrvVal", pf_p4_setGdasrvVal, METH_VARARGS},
    {"gdasrvCalcRates", pf_gdasrvCalcRates, METH_VARARGS},
    {"gdasrvCalcRates_np", pf_gdasrvCalcRates_np, METH_VARARGS},
    {"p4_setPInvarVal", pf_p4_setPInvarVal, METH_VARARGS},
    {"p4_setRelRateVal", pf_p4_setRelRateVal, METH_VARARGS},

    {"p4_setNodeRelation", pf_p4_setNodeRelation, METH_VARARGS},
    {"p4_setTreeRoot", pf_p4_setTreeRoot, METH_VARARGS},
    {"p4_setBrLen", pf_p4_setBrLen, METH_VARARGS},
    {"p4_getTreeLen", pf_p4_getTreeLen, METH_VARARGS, doc_p4_getTreeLen},

    {"p4_setCompNum", pf_p4_setCompNum, METH_VARARGS},
    {"p4_setRMatrixNum", pf_p4_setRMatrixNum, METH_VARARGS},
    {"p4_setGdasrvNum", pf_p4_setGdasrvNum, METH_VARARGS},
    //{"p4_setPreAndPostOrder", pf_p4_setPreAndPostOrder, METH_VARARGS},

    {"p4_setPrams", pf_p4_setPrams, METH_VARARGS},
    {"p4_setPramsTest", pf_p4_setPramsTest, METH_VARARGS},
    {"p4_treeLogLike", pf_p4_treeLogLike, METH_VARARGS},
    {"p4_partLogLike", pf_p4_partLogLike, METH_VARARGS},
    {"p4_calculateBigPDecks", pf_p4_calculateBigPDecks, METH_VARARGS},
    {"p4_calculateAllBigPDecksAllParts", pf_p4_calculateAllBigPDecksAllParts, METH_VARARGS},
    {"p4_setConditionalLikelihoodsOfInternalNodePart", pf_p4_setConditionalLikelihoodsOfInternalNodePart, METH_VARARGS},
    {"p4_allBOBYQAOptimize", pf_p4_allBOBYQAOptimize, METH_VARARGS},
    {"p4_newtAndBOBYQAOpt", pf_p4_newtAndBOBYQAOpt, METH_VARARGS},
    {"p4_allBrentPowellOptimize", pf_p4_allBrentPowellOptimize, METH_VARARGS},
    {"p4_newtAndBrentPowellOpt", pf_p4_newtAndBrentPowellOpt, METH_VARARGS},
    {"p4_newtSetup", pf_p4_newtSetup, METH_VARARGS},
    {"p4_getBrLens", pf_p4_getBrLens, METH_VARARGS},
    {"p4_getFreePrams", pf_p4_getFreePrams, METH_VARARGS},
    {"p4_getRelRate", pf_p4_getRelRate, METH_VARARGS},
    {"p4_simulate", pf_p4_simulate, METH_VARARGS},
    {"p4_drawAncState", pf_p4_drawAncState, METH_VARARGS, doc_p4_drawAncState},
    {"p4_expectedCompositionCounts", pf_p4_expectedCompositionCounts, METH_VARARGS, doc_p4_expectedCompositionCounts},
    {"p4_expectedComposition", pf_p4_expectedComposition, METH_VARARGS, doc_p4_expectedComposition},

    {"p4_verifyIdentityOfTwoTrees", pf_p4_verifyIdentityOfTwoTrees, METH_VARARGS},
    {"p4_copyCondLikes", pf_p4_copyCondLikes, METH_VARARGS},
    {"p4_copyBigPDecks", pf_p4_copyBigPDecks, METH_VARARGS},
    {"p4_copyModelPrams", pf_p4_copyModelPrams, METH_VARARGS},

    {"zeroNumPyInts", pf_zeroNumPyInts, METH_VARARGS},
    {"logDetFillFxy", pf_logDetFillFxy, METH_VARARGS},
    {"effectiveSampleSize", pf_effectiveSampleSize, METH_VARARGS},
    {"newtonRaftery94_eqn16", pf_newtonRaftery94_eqn16, METH_VARARGS},
    {"setMcmcTreeCallback", pf_setMcmcTreeCallback, METH_VARARGS},
    {"unsetMcmcTreeCallback", pf_unsetMcmcTreeCallback, METH_VARARGS},

    {"test", pf_test, METH_VARARGS},

    {NULL, NULL}
};


static struct PyModuleDef pfmodule = {
   PyModuleDef_HEAD_INIT,
   "pf",   /* name of module */
   NULL, /* module documentation, may be NULL */
   -1,       /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
   pfMethods
};

PyMODINIT_FUNC
PyInit_pf(void)
{
    return PyModule_Create(&pfmodule);
}

