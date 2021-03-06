#include <Python.h>
#include <numpy/arrayobject.h>

p4_model *p4_newModel(int nParts, int doRelRates, int relRatesAreFree, int nFreePrams, int isHet, int *rMatrixNormalizeTo1, double *PINVAR_MIN, double *PINVAR_MAX, double *KAPPA_MIN, double *KAPPA_MAX, double *GAMMA_SHAPE_MIN, double *GAMMA_SHAPE_MAX, double *PIVEC_MIN, double *PIVEC_MAX, double *RATE_MIN, double *RATE_MAX, double *RELRATE_MIN, double *RELRATE_MAX, double *BRLEN_MIN, double *BRLEN_MAX);
void p4_freeModel(p4_model *aModel);
void p4_dumpModel(p4_model *aModel);
void p4_newModelPart(p4_model *aModel, int pNum, int dim, int nComps, int nRMatrices, int nGdasrvs, int nCat, int pInvarFree, int *bQETneedsReset);
void p4_freeModelPart(p4_modelPart *aModelPart);
void p4_newComp(p4_model *aModel, int pNum, int mNum, int free, PyArrayObject *val);
void p4_newRMatrix(p4_model *aModel, int pNum, int mNum, int free, int spec);
p4_gdasrv *p4_newGdasrv(p4_model *aModel, int pNum, int mNum, int nCat, int free, PyArrayObject *val, PyArrayObject *freqs, PyArrayObject *rates);
void p4_resetBQET(p4_model *aModel, int pNum, int compNum, int rMatrixNum);

