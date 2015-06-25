#include <Python.h>
#include <numpy/arrayobject.h>

p4_model *p4_newModel(int nParts, int doRelRates, int relRatesAreFree, int nFreePrams, int isHet, int *rMatrixNormalizeTo1);
void p4_freeModel(p4_model *aModel);
void p4_dumpModel(p4_model *aModel);
void p4_newModelPart(p4_model *aModel, int pNum, int dim, int nComps, int nRMatrices, int nGdasrvs, int nCat, int isMixture, int mixtureIsFree, int pInvarFree, int doTSCovarion, int tSCovIsFree, int *bQETneedsReset);
void p4_freeModelPart(p4_modelPart *aModelPart);
void p4_newComp(p4_model *aModel, int pNum, int mNum, int free);
void p4_newRMatrix(p4_model *aModel, int pNum, int mNum, int free, int spec);
p4_gdasrv *p4_newGdasrv(p4_model *aModel, int pNum, int mNum, int nCat, int free, PyArrayObject *val, PyArrayObject *freqs, PyArrayObject *rates);
void p4_resetBQET(p4_model *aModel, int pNum, int compNum, int rMatrixNum);

