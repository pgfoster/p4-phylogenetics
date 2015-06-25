#include "Python.h"

// p4_tree.c
p4_tree *p4_newTree(int nNodes, int nLeaves, int *preOrder, int *postOrder, double *partLikes, data *aData, p4_model *aModel);
void p4_freeTree(p4_tree *aTree);
void p4_dumpTree(p4_tree *aTree);
void p4_setPrams(p4_tree *aTree);
void p4_setPramsPart(p4_tree *aTree, int pNum);
void p4_calculateAllBigPDecksAllParts(p4_tree *aTree);
void p4_setPramsTest(p4_tree *aTree);
void p4_setPramsPartTest(p4_tree *aTree, int pNum);
double p4_treeLogLike(p4_tree *aTree, int getSiteLikes);
double p4_partLogLike(p4_tree *aTree, part *p, int partNum, int getSiteLikes);
double p4_partLogLikeSiteRates(p4_tree *aTree, part *p, int pNum, int getSiteLikes, double *siteRates, int *gammaCats, double *work);
void p4_getPreOrderNodeNumsAbove(p4_tree *aTree, p4_node *aNode);

// p4_treeOpt.c
void p4_windUpParameters(p4_tree *aTree, double *parameters, double *lBounds, double *uBounds, int *compStarts);
PyObject *p4_getBrLens(p4_tree *aTree);
PyObject *p4_getFreePrams(p4_tree *aTree);
void p4_unWindParameters(p4_tree *aTree, double *parameters);
double p4_minusLogLikeForBrent(double *parameters);
void p4_allBrentPowellOptimize(p4_tree *aTree);
void p4_simplexOptimize(p4_tree *aTree, PyObject *treeObject, PyObject *theMethod);
void p4_newtAndBrentPowellOpt(p4_tree *aTree);
void p4_newtAnd1DBrent(p4_tree *aTree);
void p4_bracketTheMinimum(p4_tree *aTree, double *guess1, double *guess2, double *lBound, double *uBound, double *pram, double downFactor);
typedef double (*MinimizeFxn)(double *);
//typedef double (*MinimizeFxn)(double *);
//double LocalMin(double a, double b, double eps, double t, MinimizeFxn f, double *px);
double LocalMin(double a, double b, double eps, double t, double (*f)(double *), double *px);
void BracketMinimum(double *pA, double *pB, double (*func)(double []), double minAllowed, double maxAllowed);

// p4_treeSim.c
void p4_simulate(p4_tree *t);
PyObject *p4_expectedCompositionCounts(p4_tree *t, int partNum);
PyObject *p4_expectedComposition(p4_tree *t);

// p4_treeNewt.c
void p4_newtSetup(p4_tree *aTree);
//void p4_newt(p4_tree *aTree);
void p4_newtAround(p4_tree *aTree, double epsilon, double likeDelta);
void p4_newtNode(p4_node *aNode, double epsilon);
void p4_setNodeCL2(p4_tree *aTree, p4_node *aNode);


