#include  "pftypes.h"
#include  "p4_tree.h"
#include  "pmatrices.h"
#include  "util.h"
#include "brent.h"
#include "defines.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <nlopt.h>

static  int	nPrams = 0;
static	p4_tree *thisTree = NULL;
static  int doBranchLengths = 0;
static  int likelihoodEvaluations = 0;



void p4_windUpParameters(p4_tree *aTree, double *parameters, double *lBounds, double *uBounds)
{
    int pos = 0;
    int i, j, mtNum, pNum;
    p4_modelPart  *mp;
    p4_comp       *c;
    p4_rMatrix    *r;
    p4_gdasrv     *g;
    p4_node       *n;

    //printf("p4_windUpParameters() here.\n");
    for(pNum = 0; pNum < aTree->nParts; pNum++) {
        mp = aTree->model->parts[pNum];

        // comps
        for(mtNum = 0; mtNum < mp->nComps; mtNum++) {
            c = mp->comps[mtNum];
            if(c->free) {
                for(i = 0; i < mp->dim - 1; i++) {
                    parameters[pos] = c->val[i];
                    if(lBounds) {
                        lBounds[pos] = aTree->model->PIVEC_MIN[0];
                        uBounds[pos] = aTree->model->PIVEC_MAX[0];
                    }
                    pos++;
                }
            }
        }

        // rMatrices
        for(mtNum = 0; mtNum < mp->nRMatrices; mtNum++) {
            r = mp->rMatrices[mtNum];
            if(r->free) {
                if(r->spec == RMATRIX_2P) {
                    parameters[pos] = r->kappa[0];
                    if(lBounds) {
                        lBounds[pos] = aTree->model->KAPPA_MIN[0];
                        uBounds[pos] = aTree->model->KAPPA_MAX[0];
                    }
                    pos++;
                }
                else {
                    for(i = 0; i < mp->dim - 2; i++) {
                        for(j = i + 1; j < mp->dim; j++) {
                            parameters[pos] = r->bigR[i][j];
                            if(lBounds) {
                                lBounds[pos] = aTree->model->RATE_MIN[0];
                                uBounds[pos] = aTree->model->RATE_MAX[0];
                            }
                            pos++;
                        }
                    }
                }
            }
        }

        // gdasrvs
        for(mtNum = 0; mtNum < mp->nGdasrvs; mtNum++) {
            g = mp->gdasrvs[mtNum];
            if(g->free) {
                parameters[pos] = g->val[0];
                if(lBounds) {
                    lBounds[pos] = aTree->model->GAMMA_SHAPE_MIN[0];
                    uBounds[pos] = aTree->model->GAMMA_SHAPE_MAX[0];
                }
                pos++;
            }
        }

        // pInvar
        if(mp->pInvar->free) {
            parameters[pos] = mp->pInvar->val[0];
            if(lBounds) {
                lBounds[pos] = aTree->model->PINVAR_MIN[0];
                uBounds[pos] = aTree->model->PINVAR_MAX[0];
            }
            pos++;
        }

    }

    // relRates are done after the stuff above.  Ie out of the 'parts' loop.
    if(aTree->model->relRatesAreFree) {
        for(pNum = 0; pNum < aTree->nParts - 1; pNum++) {
            mp = aTree->model->parts[pNum];
            parameters[pos] = mp->relRate[0];
            if(lBounds) {
                lBounds[pos] = aTree->model->RELRATE_MIN[0];
                uBounds[pos] = aTree->model->RELRATE_MAX[0];
            }
            pos++;
        }
    }

    // brLens
    if(doBranchLengths) {
        for(i = 0; i < aTree->nNodes; i++) {
            n = aTree->nodes[i];
            if(n != aTree->root) {
                parameters[pos] = n->brLen[0];
                if(lBounds) {
                    lBounds[pos] = aTree->model->BRLEN_MIN[0];
                    uBounds[pos] = aTree->model->BRLEN_MAX[0];
                }
                pos++;
            }
        }
    }
	
    // summary
    if((0)) {
        pos = 0;
        for(i = 0; i < aTree->model->nFreePrams; i++) {
            printf("parameters[%i] = %f\n", i, parameters[i]);
            pos++;
        }
        for(i = aTree->model->nFreePrams; i < aTree->model->nFreePrams + (aTree->nNodes - 1); i++) {
            printf("brLen pram[%i] = %f\n", i, parameters[i]);
            pos++;
        }
        //exit(1);
    }
}




PyObject *p4_getBrLens(p4_tree *aTree)
{
    int i;
    PyObject        *theList;
 
    theList = PyList_New(aTree->nNodes);
    for(i = 0; i < aTree->nNodes; i++){
        if(aTree->nodes[i] != aTree->root) {
            PyList_SetItem(theList, i, PyFloat_FromDouble(aTree->nodes[i]->brLen[0]));
        }
    }
    PyList_SetItem(theList, aTree->root->nodeNum, PyFloat_FromDouble(-1.0));  // for the root node.  A placeholder
    return theList;
}


PyObject *p4_getFreePrams(p4_tree *aTree)
{

    PyObject        *theList;
    double *prams;
    int i;
 
    // We want aTree->model->nFreePrams, not nPrams
    theList = PyList_New(aTree->model->nFreePrams);
    prams = malloc(aTree->model->nFreePrams * sizeof(double));
    doBranchLengths = 0;
 
    p4_windUpParameters(aTree, prams, NULL, NULL);
    for(i = 0; i < aTree->model->nFreePrams; i++){
        PyList_SetItem(theList, i, PyFloat_FromDouble(prams[i]));
    }

#if 0
    printf("p4_getFreePrams()\n");
    for(i = 0; i < aTree->model->nFreePrams; i++){
        printf("    %2i   %6.4f\n", i, prams[i]);
    }
#endif

    free(prams);
    return theList;
}



void p4_unWindParameters(p4_tree *aTree, double *parameters)
{
    int pos = 0;
    int i, j, mtNum, pNum, totLen;
    double     sum, sum2;
    double     par;
    double     reducedUpper;
    double     diff;
    double     factor;
    int hit_limit = 0;
    p4_modelPart  *mp = NULL;
    p4_comp       *c = NULL;
    p4_rMatrix    *r = NULL;
    p4_gdasrv     *g = NULL;
    p4_node       *n = NULL;

    //printf("p4_unWindParameters() here.\n");
    // Do relRates after this loop.
    for(pNum = 0; pNum < aTree->nParts; pNum++) {
        mp = aTree->model->parts[pNum];

        // comps
        for(mtNum = 0; mtNum < mp->nComps; mtNum++) {
            c = mp->comps[mtNum];
            if(c->free) {
                // Here we restore the dim-1 parameters, but
                // effectively we move the zero point to be PIVEC_MIN.

                // loop up to less than dim - 1, not dim
                sum = 0.0;
                for(i = 0; i < mp->dim - 1; i++) {
                    par = parameters[pos];
                    if(par <= aTree->model->PIVEC_MIN[0]) {
                        hit_limit = 1;
                        par = 0.0;
                    } else {
                        par = par - aTree->model->PIVEC_MIN[0];
                    }
                    c->val[i] = par;
                    sum += par;
                    pos++;
                }

                // Make a new target for the sum, that takes into account our temporary "zero".
                reducedUpper = 1.0 - (mp->dim * aTree->model->PIVEC_MIN[0]);
                // Assign a value to c->val[dim-1]
                diff = reducedUpper - sum;
                if(diff < 0.0) {
                    c->val[mp->dim - 1] = 0.0;
                } else {
                    c->val[mp->dim - 1] = diff;
                    sum += diff;
                }

                // Change all the values in proportion to the sum, and then add back PIVEC_MIN
                factor = reducedUpper / sum;
                
                for(i = 0; i < mp->dim; i++) {
                    c->val[i] *= factor;
                    c->val[i] += aTree->model->PIVEC_MIN[0];
                }
                
                if((1)) {
                    sum2 = 0.0;
                    for(i = 0; i < mp->dim; i++) {
                        sum2 += c->val[i];
                    }
                    if(fabs(1.0 - sum2) > 1.0e-14) {
                        printf("%%%% p4_unWindParameters()  part %i, comp %i, values do not sum to 1.0.  sum2=%g\n", 
                               pNum, mtNum, sum);
                        printf("%%%%  sum2 - 1.0 = %g\n", sum2 - 1.0);
                        printf("%%%% sum was %g  %f\n", sum, sum);
                        printf("%%%% fabs(sum - 1.0) = %g\n", fabs(sum - 1.0));
                        //exit(1);
                    }
                }
                if((1)) {
                    for(i = 0; i < mp->dim; i++) {
                        if(c->val[i] < aTree->model->PIVEC_MIN[0]) {
                            printf("%%%% p4_unWindParameters()  part %i, comp %i, val[%i] %g smaller than PIVEC_MIN %g\n", 
                               pNum, mtNum, i, c->val[i], aTree->model->PIVEC_MIN[0]);
                        }
                    }
                }
            }

#if 0
            if(c->free) {
                sum = 0.0;
                for(i = 0; i < mp->dim - 1; i++) {
                    par = parameters[pos];
                    if(par < aTree->model->PIVEC_MIN[0]) {
                        //hit_limit = 1;
                        par = (aTree->model->PIVEC_MIN[0] * (mp->dim - 1)) + (aTree->model->PIVEC_MIN[0] * ranDoubleUpToOne());
                    } else if(par > aTree->model->PIVEC_MAX[0]) {
                        //hit_limit = 1;
                        par = aTree->model->PIVEC_MAX[0];
                    }
                    c->val[i] = par;
                    sum += par;
                    pos++;
                }
                // At this point, sum might be only a few times PIVEC_MIN, or it might be dim-1 * PIVEC_MAX
                if(sum < 1.0) {
                    c->val[mp->dim - 1] = 1.0 - sum;
                    if(c->val[mp->dim - 1] < aTree->model->PIVEC_MIN[0]) {
                        hit_limit = 1;
                        c->val[mp->dim - 1] = aTree->model->PIVEC_MIN[0] + (aTree->model->PIVEC_MIN[0] * ranDoubleUpToOne());
                    } else if(c->val[mp->dim - 1] > aTree->model->PIVEC_MAX[0]) {
                        hit_limit = 1;
                        c->val[mp->dim - 1] = aTree->model->PIVEC_MAX[0];
                    }
                    sum += c->val[mp->dim - 1];
                    //printf("    a  %f  %f  %f  %f  %f\n", c->val[0], c->val[1], c->val[2], c->val[3], sum);
                    if(sum != 1.0) {  // This will happen if c->val[mp->dim - 1] == PIVEC_MAX
                        //printf("xx unwind. sum = %f  Don't panic.\n", sum);
                        for(i = 0; i < mp->dim; i++) {
                            c->val[i] = c->val[i] / sum;
                        }
                        //printf("    X  %f  %f  %f  %f  %f\n", c->val[0], c->val[1], c->val[2], c->val[3], sum);
                    }
                }
                else if(sum >= 1.0) {  
                    // Changed 12 June 2005, from "sum > 1.0", to
                    // correct a bug where sum was exactly 1.0, which
                    // meant that c-val[dim-1] was not taken into
                    // account at all.

                    // we multiply PIVEC_MIN by sum, so that it will
                    // not go below PIVEC_MIN when it is divided by
                    // sum, immediately below.
                    c->val[mp->dim - 1] = (aTree->model->PIVEC_MIN[0] * sum) + (aTree->model->PIVEC_MIN[0] * ranDoubleUpToOne());
                    sum += c->val[mp->dim - 1];
                    for(i = 0; i < mp->dim; i++) {
                        c->val[i] = c->val[i] / sum;
                    }
                }
                if((1)) {
                    sum2 = 0.0;
                    for(i = 0; i < mp->dim; i++) {
                        sum2 += c->val[i];
                    }
                    if(fabs(1.0 - sum2) > 1.0e-14) {
                        printf("%%%% p4_unWindParameters()  part %i, comp %i, values do not sum to 1.0.  sum2=%g\n", 
                               pNum, mtNum, sum);
                        printf("%%%%  sum2 - 1.0 = %g\n", sum2 - 1.0);
                        printf("%%%% sum was %g  %f\n", sum, sum);
                        printf("%%%% fabs(sum - 1.0) = %g\n", fabs(sum - 1.0));
                        //exit(1);
                    }
                }
                if((1)) {
                    for(i = 0; i < mp->dim; i++) {
                        if(c->val[i] < aTree->model->PIVEC_MIN[0]) {
                            printf("%%%% p4_unWindParameters()  part %i, comp %i, val[%i] smaller than PIVEC_MIN %g\n", 
                               pNum, mtNum, i, aTree->model->PIVEC_MIN[0]);
                        }
                    }
                }
            }
#endif
        }

        // rMatrices
        for(mtNum = 0; mtNum < mp->nRMatrices; mtNum++) {
            r = mp->rMatrices[mtNum];
            if(r->free) {
                if(r->spec == RMATRIX_2P) {
                    if(parameters[pos] < aTree->model->KAPPA_MIN[0]) {
                        hit_limit = 1;
                        parameters[pos] = aTree->model->KAPPA_MIN[0];
                    }
                    else if(parameters[pos] > aTree->model->KAPPA_MAX[0]) {
                        hit_limit = 1;
                        parameters[pos] = aTree->model->KAPPA_MAX[0];
                    }
                    r->kappa[0] = parameters[pos];
                    pos++;
                }
                else {
                    if(aTree->model->rMatrixNormalizeTo1[0]) {
                        sum = 0.0;
                        for(i = 0; i < mp->dim - 2; i++) {
                            for(j = i + 1; j < mp->dim; j++) {
                                par = parameters[pos];
                                if(par < aTree->model->RATE_MIN[0]) {
                                    par = aTree->model->RATE_MIN[0] + (aTree->model->RATE_MIN[0] * ranDoubleUpToOne());
                                } else if(par > 0.999) {
                                    par = 0.999;
                                }
                                r->bigR[i][j] = par;
                                sum += par;
                                pos++;
                            }
                        }
                        // At this point, sum could be anywhere from a few times RATE_MIN to much greater than 1.0
                        if(sum < 1.0) {
                            //printf("sum < 1.0 ");
                            par = 1.0 - sum;
                            if(par < aTree->model->RATE_MIN[0]) {
                                hit_limit = 1;
                                par = aTree->model->RATE_MIN[0] + (aTree->model->RATE_MIN[0] * ranDoubleUpToOne());
                            } else if(par > 0.999) {
                                hit_limit = 1;
                                par = 0.999;
                            }
                            sum += par;
                            r->bigR[mp->dim - 2][mp->dim - 1] = par;
                            if(sum != 1.0) {
                                for(i = 0; i < mp->dim - 1; i++) {
                                    for(j = i + 1; j < mp->dim; j++) {
                                        r->bigR[i][j] /= sum;
                                    }
                                }
                            }
                            for(i = 0; i < mp->dim - 1; i++) {
                                for(j = i + 1; j < mp->dim; j++) {
                                    r->bigR[j][i] = r->bigR[i][j];
                                }
                            }
                        } else if(sum >= 1.0) {
                            //if(fabs(sum - 1.0) < 1.e-15) {
                            //	printf("equal  ");
                            //} else {
                            //	printf("more than  ");
                            //}
                            par = (aTree->model->RATE_MIN[0] * sum) + (aTree->model->RATE_MIN[0] * ranDoubleUpToOne());
                            sum += par;
                            r->bigR[mp->dim - 2][mp->dim - 1] = par;
                            for(i = 0; i < mp->dim - 1; i++) {
                                for(j = i + 1; j < mp->dim; j++) {
                                    r->bigR[i][j] /= sum;
                                }
                            }
                            for(i = 0; i < mp->dim - 1; i++) {
                                for(j = i + 1; j < mp->dim; j++) {
                                    r->bigR[j][i] = r->bigR[i][j];
                                }
                            }
                        }
#if 0
                        sum = 0.0;
                        for(i = 0; i < mp->dim - 1; i++) {
                            for(j = i + 1; j < mp->dim; j++) {
                                sum += r->bigR[i][j];
                            }
                        }
                        printf("rMatrix sum = %f\n", sum);
#endif
						
                    }
                    else {  // not normalized to 1
                        for(i = 0; i < mp->dim - 2; i++) {
                            for(j = i + 1; j < mp->dim; j++) {
                                if(parameters[pos] < aTree->model->RATE_MIN[0]) {
                                    hit_limit = 1;
                                    parameters[pos] = aTree->model->RATE_MIN[0];
                                }else if (parameters[pos] > aTree->model->RATE_MAX[0]) {
                                    hit_limit = 1;
                                    parameters[pos] = aTree->model->RATE_MAX[0];
                                }
                                r->bigR[i][j] = parameters[pos];
                                r->bigR[j][i] = parameters[pos];
                                pos++;
                            }
                        }
                    }
                }
            }
        }

        // gdasrvs
        for(mtNum = 0; mtNum < mp->nGdasrvs; mtNum++) {
            g = mp->gdasrvs[mtNum];
            if(g->free) {
                if(parameters[pos] < aTree->model->GAMMA_SHAPE_MIN[0]) {
                    hit_limit = 1;
                    //printf("    gamma min\n");
                    parameters[pos] = aTree->model->GAMMA_SHAPE_MIN[0];
                } else if(parameters[pos] > aTree->model->GAMMA_SHAPE_MAX[0]) {
                    hit_limit = 1;
                    //printf("    gamma max\n");
                    parameters[pos] = aTree->model->GAMMA_SHAPE_MAX[0];
                }
                g->val[0] = parameters[pos];
                pos++;
            }
        }

        // pInvar
        if(mp->pInvar->free) {
            if(parameters[pos] < aTree->model->PINVAR_MIN[0]) {
                hit_limit = 1;
                //printf("    pInvar min\n");
                parameters[pos] = aTree->model->PINVAR_MIN[0];
            } else if(parameters[pos] > aTree->model->PINVAR_MAX[0]) {
                hit_limit = 1;
                //printf("    pInvar max\n");
                parameters[pos] = aTree->model->PINVAR_MAX[0];
            }
            mp->pInvar->val[0] = parameters[pos];   
            pos++;
        }

    }



    // relRate
    if(aTree->model->relRatesAreFree) {
        for(pNum = 0; pNum < aTree->nParts - 1; pNum++) {
            mp = aTree->model->parts[pNum];
            if(parameters[pos] < aTree->model->RELRATE_MIN[0]) {
                hit_limit = 1;
                parameters[pos] = aTree->model->RELRATE_MIN[0];
            } else if(parameters[pos] > aTree->model->RELRATE_MAX[0]) {
                hit_limit = 1;
                parameters[pos] = aTree->model->RELRATE_MAX[0];
            }
            mp->relRate[0] = parameters[pos];
            pos++;
        }

        // Adjust so that the weighted average of the relative rates is 1.0
        // First find the total data length.
        totLen = 0;
        for(pNum = 0; pNum < aTree->nParts; pNum++) {
            totLen +=  aTree->data->parts[pNum]->nChar;
        }
        // What we want here is to take the sum of the (rates *
        // nChar) for all the parts except the last one, and then
        // subtract that from the total length, and divide the
        // difference by the nChar for the last part.
        sum = 0.0;
        for(pNum = 0; pNum < aTree->nParts - 1; pNum++) {
            sum = sum + (aTree->model->parts[pNum]->relRate[0] * aTree->data->parts[pNum]->nChar);
        }
        sum = (((double) totLen) -  sum) / aTree->data->parts[aTree->nParts - 1]->nChar;
        aTree->model->parts[aTree->nParts -1]->relRate[0] = sum;


#if 1
        sum = 0.0;
        //printf("        unwind: \n");
        for(pNum=0; pNum < aTree->nParts; pNum++) {
            //printf("   %i    %f    %f\n", aTree->data->parts[pNum]->nChar, 
            //	   ((double)aTree->data->parts[pNum]->nChar) / ((double) totLen), 
            //	   aTree->model->parts[pNum]->relRate[0]);
            sum += aTree->model->parts[pNum]->relRate[0] * (((double)aTree->data->parts[pNum]->nChar) / ((double) totLen));
        }
        //printf("sum = %f\n", sum);
        if(sum < 0.999999999999 || sum > 1.000000000001) {
            printf("unwind. relRates sum to %19.17f, should sum to 1.0.  Bad.", sum);
            exit(1);
        }
#endif
    }

    // brLens
    if(doBranchLengths) {
        for(i = 0; i < aTree->nNodes; i++) {
            n = aTree->nodes[i];
            if(n != aTree->root) {
                if(parameters[pos] < aTree->model->BRLEN_MIN[0]) {  // this week, its 1.0e-8
                    hit_limit = 1;
                    parameters[pos] = aTree->model->BRLEN_MIN[0] - (aTree->model->BRLEN_MIN[0] * 0.1) + ((aTree->model->BRLEN_MIN[0] * 0.1) * ranDoubleUpToOne());
                } else if(parameters[pos] > aTree->model->BRLEN_MAX[0]) {
                    hit_limit = 1;
                    parameters[pos] = aTree->model->BRLEN_MAX[0];
                }
                n->brLen[0] = parameters[pos];
                pos++;
            }
        }
    }
	
    if(hit_limit) {
        //printf("p4_unWindParameters() got hit_limit\n");
        p4_windUpParameters(aTree, parameters, NULL, NULL);
    }
}

//==============================================
// NLOPT_LN_BOBYQA
//==============================================


double p4_logLikeForNLOpt(unsigned nPrams, const double *parameters, double *grad, void *my_func_data)
{
    double logLike;

    p4_unWindParameters(thisTree, (double *)parameters);

#if 0
    {
        int i;
        FILE *dfile;

        if ((dfile = fopen("debugFile", "a+")) == NULL) {
            printf(" error: Couldn't open debugFile for appending \n");
            exit(1);
        }

        fprintf(dfile, "p4_treeOpt.c: p4_logLikeForNLOpt:\n");
        for(i = 0; i < nPrams; i++) {
            fprintf(dfile, "         parameters[%i] is %.12f\n", i, parameters[i]);
        }
        fprintf(dfile, "\n");
        for(i = 0; i < thisTree->nNodes; i++) {
            fprintf(dfile, "         brLen[%i] is %.12f\n", i, 
                    thisTree->nodes[i]->brLen[0]);
        }
        fclose(dfile);
    }
#endif
	
    //printf("  about to setPrams ... \n");
    p4_setPrams(thisTree);
    logLike = p4_treeLogLike(thisTree, 0);
    likelihoodEvaluations++;
    //printf("        logLikeForNLOpt finished.  logLike = %12.8f\n", logLike);

    return logLike;
}

void p4_allBOBYQAOptimize(p4_tree *aTree)
{


    int i;
    double  logLike = 0.0;
    //double  previousLogLike = 0.0;
    //double  diff;
    //int totalLikelihoodEvals = 0;
    double	*parameters = NULL;
    double	*lBounds = NULL;
    double	*uBounds = NULL;
    nlopt_result result;

    // You can't pass variables to minusLogLikeForBrent,
    // so I have to set file-wide variables.
    thisTree = aTree;
    doBranchLengths = 1;
    nPrams = aTree->model->nFreePrams + (aTree->nNodes - 1);

    if((0)) {
        printf("Starting p4_allBOBYQAOptimize: nPrams is %i\n", nPrams);
        printf("...about to windUpParameters.  PIVEC_MIN is %g\n", aTree->model->PIVEC_MIN[0]);
    }
    parameters = malloc(nPrams * sizeof(double));
    if(!parameters) {
        printf("Failed to allocate memory for opt parameters.\n");
        exit(1);
    }


    lBounds = malloc(nPrams * sizeof(double));
    if(!lBounds) {
        printf("Failed to allocate memory for opt lBounds.\n");
        exit(1);
    }

    uBounds = malloc(nPrams * sizeof(double));
    if(!uBounds) {
        printf("Failed to allocate memory for opt uBounds.\n");
        exit(1);
    }

    p4_windUpParameters(aTree, parameters, lBounds, uBounds);
    if((0)) {
        printf("p4_treeOptimize.c: allBOBYQAOptimize() starting with these params\n");
        printf("nPrams is %i\n", nPrams);
        for(i = 0; i < nPrams; i++) {
            printf("         parameters[%i] is %.12f\n", i, parameters[i]);
        }
        if((0)) {
            printf("\n");
            for(i = 0; i < aTree->nNodes; i++) {
                if(aTree->nodes[i] != aTree->root) {
                    printf("         branchLength[%i] is %.12f\n", i, aTree->nodes[i]->brLen[0]);
                }
            }
        }
    }
	
    // set up optimizer
    nlopt_opt opter=nlopt_create(NLOPT_LN_BOBYQA, nPrams);
    // lower and upper bound
    nlopt_set_lower_bounds(opter, lBounds);
    nlopt_set_upper_bounds(opter, uBounds);
    // objective function; two kinds, min and max.  Here I use max
    nlopt_set_max_objective(opter, p4_logLikeForNLOpt, NULL);
    // Stopping criterion, absolute, for like value.  Stop if it does not improve by this much.
    nlopt_set_ftol_abs(opter, 1.e-6);

    result = nlopt_optimize(opter, parameters, &logLike);
    if(result < 0 && (result != NLOPT_ROUNDOFF_LIMITED)) {
        printf("a result %i BOBYQA\n", (int)result);
    } 
    //else {
    //    printf("a result %i OK BOBYQA\n", (int)result);
    //}

    // repeat
    if(opter) {
        nlopt_destroy(opter);
        opter = NULL;
    }
    p4_unWindParameters(aTree, parameters);
    p4_setPrams(aTree);
    logLike = p4_treeLogLike(aTree, 0);
    likelihoodEvaluations++;
    p4_windUpParameters(aTree, parameters, lBounds, uBounds);

    opter=nlopt_create(NLOPT_LN_BOBYQA, nPrams);
    // lower and upper bound
    nlopt_set_lower_bounds(opter, lBounds);
    nlopt_set_upper_bounds(opter, uBounds);
    // objective function; two kinds, min and max.  Here I use max
    nlopt_set_max_objective(opter, p4_logLikeForNLOpt, NULL);
    // Stopping criterion, absolute, for like value.  Stop if it does not improve by this much.
    nlopt_set_ftol_abs(opter, 1.e-8);

    result = nlopt_optimize(opter, parameters, &logLike);
    if(result < 0 && (result != NLOPT_ROUNDOFF_LIMITED)) {
        printf("b result %i BOBYQA\n", (int)result);
    } 
    //else {
    //    printf("b result %i OK BOBYQA\n", (int)result);
    //}
    
     

    //printf("allBOBYQAOptimize(). %i full likelihood evaluations\n", totalLikelihoodEvals);
    //printf("allBOBYQAOptimize(). finished\n");

    if(opter) {
        nlopt_destroy(opter);
    }
    if(parameters) {
        free(parameters);
        parameters = NULL;
    }
    if(lBounds) {
        free(lBounds);
        lBounds = NULL;
    }
    if(uBounds) {
        free(uBounds);
        uBounds = NULL;
    }

}




void p4_newtAndBOBYQAOpt(p4_tree *aTree)
{
    double  logLike = 0.0;
    double  previousLogLike = 0.0;
    double  diff;
    int pass;
    int totalLikelihoodEvals = 0;
    double	*parameters = NULL;
    double	*lBounds = NULL;
    double	*uBounds = NULL;
    nlopt_result result;
    nlopt_opt    opter;

    // Rather than passing aTree to to logLikeForNLOpt,
    // I set file-wide variables; leftover from Brent
    thisTree = aTree;
    doBranchLengths = 0;
    nPrams = aTree->model->nFreePrams;
    
    printf("Starting p4_newtAndBOBYQAOpt().   nPrams is %i\n", nPrams);

    if(nPrams == 0) {
        p4_newtAround(aTree, 1.0, 10.0);
        p4_newtAround(aTree, 1.0e-1, 1.0);
        p4_newtAround(aTree, 1.0e-2, 0.1);
        p4_newtAround(aTree, 1.0e-5, 1.0e-7);
        return;
    }

    // If we are here, then nPrams > 0, so we use BOBYQA with newt.
    parameters = malloc(nPrams * sizeof(double));
    if(!parameters) {
        printf("Failed to allocate memory for opt parameters.\n");
        exit(1);
    }

    lBounds = malloc(nPrams * sizeof(double));
    if(!lBounds) {
        printf("Failed to allocate memory for opt lBounds.\n");
        exit(1);
    }

    uBounds = malloc(nPrams * sizeof(double));
    if(!uBounds) {
        printf("Failed to allocate memory for opt uBounds.\n");
        exit(1);
    }

    p4_windUpParameters(aTree, parameters, lBounds, uBounds);



    p4_newtAround(aTree, 1.0, 10.0);
    p4_newtAround(aTree, 1.0e-1, 1.0);
    p4_newtAround(aTree, 1.0e-5, 1.0e-7);

    likelihoodEvaluations = 0;
    previousLogLike = p4_treeLogLike(aTree, 0);
    printf("  x logLike = %f\n", previousLogLike);
    totalLikelihoodEvals += 1;	
    pass = 0;

    while(1) {
        //double prev;
        //double afterNewtLogLike;
        //double diff2;

        likelihoodEvaluations = 0;
        //printf("....About to newtAround ....\n");
        //prev = p4_treeLogLike(aTree, 0);
        p4_newtAround(aTree, 1.0e-5, 1.0e-7);
        //afterNewtLogLike = p4_treeLogLike(aTree, 0);
        //diff2 = afterNewtLogLike - prev;
        //printf("                          newt diff =     %.8f\n", diff2);
        //printf("....About to BOBYQA.  pass=%i\n", pass);
        

        // Initializing the optimizer inside the loop is in theory
        // inefficient; I should intialize it first, and then call
        // nlopt_optimize repeatedly in the loop.  And I used to do
        // that.  But I am doing it this way, inside the loop so that
        // I can "freshen" the paramters via p4_unWindParameters().

        // set up optimizer
        opter=nlopt_create(NLOPT_LN_BOBYQA, nPrams);
        // lower and upper bound
        nlopt_set_lower_bounds(opter, lBounds);
        nlopt_set_upper_bounds(opter, uBounds);
        // objective function; two kinds, min and max.  Here I use max
        nlopt_set_max_objective(opter, p4_logLikeForNLOpt, NULL);
        // Stopping criterion, absolute, for like value.  Stop if it does not improve by this much.
        nlopt_set_ftol_abs(opter, 1.e-8);

        result = nlopt_optimize(opter, parameters, &logLike);
        if((result < 0) && (result != NLOPT_ROUNDOFF_LIMITED)) {
            printf("result %i Failed newtAndBOBYQA\n", (int)result);
        } 
        
        //diff2 = logLike - afterNewtLogLike;
        //printf("                        BIBYQA diff =     %.8f\n", diff2);
        //printf("....finished BOBYQA.\n");

        // I suppose that when praxis exits, the best logLike is returned,
        // and the corresponding best parameters would be in the parameters
        // vector.  However, the models may not reflect that, and may have the parameters
        // of the last likelihood that was evaluated, which may not be the best one.
        // So make the model params reflect the parameters vector, by unwinding the 
        // parameters vector, and doing a setPrams.

        p4_unWindParameters(thisTree, parameters);
        p4_setPrams(thisTree);
        diff = logLike - previousLogLike;
        if((0)) {
            if(diff < -1.0e-6) {
                printf("    p4_newtAndBOBYQAOpt().  logLike=%f, diff=%f, got worse! (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            } else {
                printf("    p4_newtAndBOBYQAOpt().  logLike=%f, diff=%f   (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            }
        }
        totalLikelihoodEvals += likelihoodEvaluations;
        previousLogLike = logLike;
        if(fabs(diff) < 1.0e-6) {
            printf("finished loop because the improvement was %.2g\n", diff);
            break;
        }
        pass++;
        if(pass > 50) {
            printf("=============================================================\n");
            printf("p4_newtAndBOBYQAOpt().  Pass limit (50) exceeded without\n");
            printf("convergence.  Giving up.  This tree is not optimized!\n");
            break;
        }
    }
				

    logLike = p4_treeLogLike(aTree, 0);
    totalLikelihoodEvals++;
    if((0)) {
        printf("     p4_newtAndBOBYQAOpt().  logLike %f\n", logLike);
        printf("     %i total likelihood evaluations\n", totalLikelihoodEvals);
        //printf("Starting the opt again.  About to windUpParameters.\n");
    }

    //printf("p4_newtAndBOBYQAOpt().  %i full likelihood evaluations\n", totalLikelihoodEvals);

#if 0
    {
        int i;

        printf("p4_treeOptimize.c: newtAndBOBYQAOpt() after optimizing\n");
        printf("nPrams is %i\n", nPrams);
        for(i = 0; i < nPrams; i++) {
            printf("         parameters[%i] is %.12f\n", i, parameters[i]);
        }
        if(1) {
            printf("\n");
            for(i = 0; i < aTree->nNodes; i++) {
                if(aTree->nodes[i] != aTree->root) {
                    printf("         branchLength[%i] is %.12f\n", i, aTree->nodes[i]->brLen[0]);
                }
            }
        }
    }
#endif

    if(opter) {
        nlopt_destroy(opter);
    }
    if(parameters) {
        free(parameters);
        parameters = NULL;
    }
    if(lBounds) {
        free(lBounds);
        lBounds = NULL;
    }
    if(uBounds) {
        free(uBounds);
        uBounds = NULL;
    }
}




//==============================================
// Brent
//==============================================


double p4_minusLogLikeForBrent(double *parameters)
{
    double logLike;

    p4_unWindParameters(thisTree, parameters);

#if 0
    {
        int i;
        FILE *dfile;

        if ((dfile = fopen("debugFile", "a+")) == NULL) {
            printf(" error: Couldn't open debugFile for appending \n");
            exit(1);
        }

        fprintf(dfile, "p4_treeOpt.c: p4_minusLogLikeForBrent:\n");
        for(i = 0; i < nPrams; i++) {
            fprintf(dfile, "         parameters[%i] is %.12f\n", i, parameters[i]);
        }
        fprintf(dfile, "\n");
        for(i = 0; i < thisTree->nNodes; i++) {
            fprintf(dfile, "         brLen[%i] is %.12f\n", i, 
                    thisTree->nodes[i]->brLen[0]);
        }
        fclose(dfile);
    }
#endif
	
    //printf("  about to setPrams ... \n");
    p4_setPrams(thisTree);
    logLike = p4_treeLogLike(thisTree, 0);
    likelihoodEvaluations++;
    //printf("        minusLogLikeForBrent finished.  logLike = %12.8f\n", logLike);

    return -logLike;
}




//==============================================
// ALL BRENT POWELL OPTIMIZE
//==============================================

void p4_allBrentPowellOptimize(p4_tree *aTree)
{


    int i;
    double  logLike = 0.0;
    double  previousLogLike = 0.0;
    double  diff;
    int verbose = 0;
    int totalLikelihoodEvals = 0;
    double	*parameters = NULL;
    brent	*aBrent = NULL;

    // You can't pass variables to minusLogLikeForBrent,
    // so I have to set file-wide variables.
    thisTree = aTree;
    doBranchLengths = 1;
	
    // count up the number of parameters
    if(doBranchLengths) {
        nPrams = aTree->model->nFreePrams + (aTree->nNodes - 1);
    } else {
        nPrams = aTree->model->nFreePrams;
    }

    if((0)) {
        printf("Starting p4_allBrentPowellOptimize: nPrams is %i\n", nPrams);
        printf("...about to windUpParameters.\n");
    }
    parameters = malloc(nPrams * sizeof(double));
    if(!parameters) {
        printf("Failed to allocate memory for opt parameters.\n");
        exit(1);
    }
    p4_windUpParameters(aTree, parameters, NULL, NULL);
    if((0)) {
        printf("p4_treeOptimize.c: allBrentPowellOptimize() starting with these params\n");
        printf("nPrams is %i\n", nPrams);
        for(i = 0; i < nPrams; i++) {
            printf("         parameters[%i] is %.12f\n", i, parameters[i]);
        }
        if((0)) {
            printf("\n");
            for(i = 0; i < aTree->nNodes; i++) {
                if(aTree->nodes[i] != aTree->root) {
                    printf("         branchLength[%i] is %.12f\n", i, aTree->nodes[i]->brLen[0]);
                }
            }
        }
    }
    //p4_unWindParameters(aTree, parameters);
    //printf("parameters[0]=%f, parameters[1]=%f\n", parameters[0], parameters[1]);
    //exit(1);
	
    aBrent = newBrent(nPrams);
    //aTree->logLike = -praxis(aBrent, 1.0e-4, 1.0, nPrams, parameters, minusLogLikeForBrent);
    previousLogLike = p4_treeLogLike(aTree, 0);
    totalLikelihoodEvals = 1;
    //previousLogLike = treeLogLike(aTree, 0);
    diff = previousLogLike;
    if(verbose) printf("Starting praxis with loglike  %f...\n", previousLogLike);
    //if(verbose) dump_pdvector(parameters, nPrams);
    //windUpParameters(aTree, parameters, NULL, NULL, NULL);

    while(fabs(diff) > 1.0e-6) {
        likelihoodEvaluations = 0; // a static var.
        // In the praxis function, after aBrent, the next two numbers are tol and h.
        // These are transferred to internal praxis variables aBrent->toler and aBrent->htol.
        // Originally 1.0e-4, and 1.0, from Dave S.
        // I tried varying h, the second number.  Lowering to 0.5 or 0.1 made it faster and better, ie
        // it got to the true ml value without stopping.  With 0.1 it got to the ml in one go, but with higher
        // values it took two calls to get to the ml value-- Why does it get stuck in a non-optimum?
        // With h = 0.01, it got to ml in one go, but was twice as slow as h=0.1
        // 0.05 - 7104 (one go)
        // 0.1 -- 5216 likelihood evaluations (one go)
        // 0.2 -- 31000 !! (two goes)
        // 0.3 -- 10,000 (two goes)
        // 0.5 -- 9500 (two goes)
        // 1.0 -- 25000 (two goes)
		
        logLike = -praxis(aBrent, 1.0e-4, 0.1, nPrams, parameters, p4_minusLogLikeForBrent);

        // I suppose that when praxis exits, the best logLike is returned,
        // and the corresponding best parameters would be in the parameters
        // vector.  However, the models may not reflect that, and may have the parameters
        // of the last likelihood that was evaluated, which may not be the best one.
        // So make the model params reflect the parameters vector, by unwinding the 
        // parameters vector, and doing a setPrams.
        p4_unWindParameters(thisTree, parameters);
        p4_setPrams(thisTree);
        //diff = previousLogLike - logLike;
        diff = logLike - previousLogLike;
        if(verbose) {
            if(diff < 0.0) {
                printf("     logLike %f, diff %f, -> got worse!  (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            } else {
                printf("     logLike %f, diff %f  (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            }
        }
        //if(verbose) dump_pdvector(parameters, nPrams);
        totalLikelihoodEvals = totalLikelihoodEvals + likelihoodEvaluations;
        previousLogLike = logLike;
    }

    logLike = p4_treeLogLike(aTree, 0);

#if 1
    if((0)) {
        printf("     logLike %f\n", logLike);
        printf("     %i total likelihood evaluations\n", totalLikelihoodEvals);
        printf("Starting the opt again.  About to windUpParameters.\n");
    }

    // I have noticed that the occasional tree gets optimized to a
    // value that is not optimal.  Dang!  If I re-optimize, it finds
    // the better optimization.  I don't know what else to do about it
    // except repeat.

    p4_windUpParameters(aTree, parameters, NULL, NULL);
    diff = 1.0;
    while(fabs(diff) > 1.0e-6) {
        likelihoodEvaluations = 0;
        logLike = -praxis(aBrent, 1.0e-4, 0.05, 
                          nPrams, parameters, p4_minusLogLikeForBrent);
        p4_unWindParameters(thisTree, parameters); 
        p4_setPrams(thisTree);
        diff = previousLogLike - logLike;
        if(verbose) {
            if(diff > 0.0) {
                printf("     logLike %f, diff %f, -> got worse!  (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            } else {
                printf("     logLike %f, diff %f  (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            }
        }
        totalLikelihoodEvals = totalLikelihoodEvals + likelihoodEvaluations;
        previousLogLike = logLike;
    }

    logLike = p4_treeLogLike(aTree, 0);
    if(verbose) {
        printf("     logLike %f\n", logLike);
        printf("     %i total likelihood evaluations\n", totalLikelihoodEvals);
    }
#endif


    //printf("allBrentPowellOptimize(). %i full likelihood evaluations\n", totalLikelihoodEvals);
    //printf("allBrentPowellOptimize(). finished\n");

#if 0
    {
        printf("p4_treeOptimize.c: allBrentPowellOptimize() after optimizing\n");
        printf("nPrams is %i\n", nPrams);
        for(i = 0; i < nPrams; i++) {
            printf("         parameters[%i] is %.12f\n", i, parameters[i]);
        }
        if(1) {
            printf("\n");
            for(i = 0; i < aTree->nNodes; i++) {
                if(aTree->nodes[i] != aTree->root) {
                    printf("         branchLength[%i] is %.12f\n", i, aTree->nodes[i]->brLen[0]);
                }
            }
        }
    }
#endif
	
    if(aBrent) {
        freeBrent(aBrent);
        aBrent = NULL;
    }
    if(parameters) {
        free(parameters);
        parameters = NULL;
    }

}





//=====================================================================
// Fast optimize, alternating newt and either Brent-Powell or 1D Brent.
//=====================================================================

void p4_newtAndBrentPowellOpt(p4_tree *aTree)
{
    double  logLike = 0.0;
    double  previousLogLike = 0.0;
    double  diff;
    //double  diff2, prev, afterNewtLogLike;
    int pass;
    int totalLikelihoodEvals = 0;
    double	*parameters = NULL;
    brent	*aBrent = NULL;
    //int  nNewts;

    // You can't pass variables to minusLogLikeForBrent,
    // so I have to set file-wide variables.
    thisTree = aTree;
    doBranchLengths = 0;
	
    // count up the number of parameters
    if(doBranchLengths) {
        nPrams = aTree->model->nFreePrams + (aTree->nNodes - 1);
    } else {
        nPrams = aTree->model->nFreePrams;
    }
	
    //printf("Starting p4_newtAndBrentPowellOpt().   nPrams is %i\n", nPrams);
    //printf("aTree->model->rMatrixNormalizeTo1 = %i\n", aTree->model->rMatrixNormalizeTo1[0]);
    //printf("bigR\n");
    //dump_psdmatrix(aTree->model->parts[0]->rMatrices[0]->bigR, aTree->model->parts[0]->dim);
	
    //printf("    aTree->model->parts[0]->gdasrvs[0]
	
    if(nPrams == 0) {
        //p4_setPrams(aTree); // is this necessary?
        //logLikeBefore = p4_treeLogLike(aTree, 0);
        //printf("p4_newt() logLikeBefore=%f\n", logLikeBefore);

        p4_newtAround(aTree, 1.0, 10.0);
        p4_newtAround(aTree, 1.0e-1, 1.0);
        p4_newtAround(aTree, 1.0e-2, 0.1);
        p4_newtAround(aTree, 1.0e-5, 1.0e-7);
        //printf("newts finished. logLike=%f\n", aTree->logLike);
        //printf("newts finished. logLike=%f, nNewts=%i\n", aTree->logLike, nNewts);
        return;
    }
    else if(nPrams == 1) {
        //printf("p4_newtAndBrentPowellOpt() nPrams=1.\n");
        p4_newtAnd1DBrent(aTree);
        return;
    }

    // If we are here, then nPrams > 1, so we use Brent-Powell with newt.
    parameters = malloc(nPrams * sizeof(double));
    if(!parameters) {
        printf("Failed to allocate memory for opt parameters.\n");
        exit(1);
    }

    p4_windUpParameters(aTree, parameters, NULL, NULL);
    aBrent = newBrent(nPrams);
    p4_newtAround(aTree, 1.0, 10.0);
    p4_newtAround(aTree, 1.0e-1, 1.0);
    p4_newtAround(aTree, 1.0e-5, 1.0e-7);

    // Try to speed things up by doing a low-res praxis.  It did not speed things up by much.
    likelihoodEvaluations = 0;


    //printf("  w logLike = %f\n", p4_treeLogLike(aTree, 0));
    //printf("...About to praxis.\n");
    //praxis(aBrent, 1.0e-4, 10.0, nPrams, parameters, p4_minusLogLikeForBrent);
    //printf("...finished praxis.\n");
    //p4_unWindParameters(thisTree, parameters);
    //p4_setPrams(thisTree);
    //printf("first praxis %i like evaluations\n", likelihoodEvaluations);
    //exit(1);
    //totalLikelihoodEvals += likelihoodEvaluations;

    previousLogLike = p4_treeLogLike(aTree, 0);
    //printf("  x logLike = %f\n", previousLogLike);
    totalLikelihoodEvals += 1;	
    pass = 0;

    // Tune the tol and h for praxis.

    // Real protein, wag+gi
    // 1.0e-4, 1.0  -->  4 passes, 438 likes
    // 1.0e-4, 0.1  -->  4 passes, 352 likes

    // Grouped aa, from real protein, gtr+gi+free comps
    // 1.0e-4, 1.0  -->  8 passes,  8563 likes
    // 1.0e-4, 0.1  -->  8 passes,  9590 likes

    // Random DNA, gtr+ig+comp, but pInvar=0 and shape=infinity
    // 1.0e-4, 1.0  -->  23 passes,   42388 likes
    // 1.0e-4, 0.1  --> >30 passes,  >51956 likes, didn't converge

    // Random DNA, no c, pInvar=0 and shape=infinity, 3 partitions, complex model, lots of asrv
    // 1.0e-4, 1.0  -->  29 passes,  199802   likes
    // 1.0e-4, 0.1  --> >30 passes, >401570   likes, didn't converge


    // Sometimes the pass limit is exceeded without convergence.  It
    // appears that the last bit of optimization is slow.  So probly
    // there is potential for tuning this thing.  One possibility
    // might be to split it into 2 loops, where the second loop has
    // different optimizer tunings.  Maybe.
    while(1) {
        likelihoodEvaluations = 0;
        //printf("....About to newtAround ....\n");
        //prev = p4_treeLogLike(aTree, 0);
        p4_newtAround(aTree, 1.0e-5, 1.0e-7);
        //afterNewtLogLike = p4_treeLogLike(aTree, 0);
        //diff2 = afterNewtLogLike - prev;
        //printf("                          newt diff =     %.8f\n", diff2);
        //printf("....About to praxis.  pass=%i\n", pass);
        logLike = -praxis(aBrent, 1.0e-4, 1.0, nPrams, parameters, p4_minusLogLikeForBrent);
        //diff2 = logLike - afterNewtLogLike;
        //printf("                        praxis diff =     %.8f\n", diff2);
        //printf("....finished praxis.\n");

        // I suppose that when praxis exits, the best logLike is returned,
        // and the corresponding best parameters would be in the parameters
        // vector.  However, the models may not reflect that, and may have the parameters
        // of the last likelihood that was evaluated, which may not be the best one.
        // So make the model params reflect the parameters vector, by unwinding the 
        // parameters vector, and doing a setPrams.

        p4_unWindParameters(aTree, parameters);
        p4_setPrams(aTree);
        diff = logLike - previousLogLike;
        if((0)) {
            if(diff < -1.0e-6) {
                printf("    p4_newtAndBrentPowellOpt().  logLike=%f, diff=%f, got worse! (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            } else {
                printf("    p4_newtAndBrentPowellOpt().  logLike=%f, diff=%f   (%i likelihoodEvaluations)\n", 
                       logLike, diff, likelihoodEvaluations);
            }
        }
        totalLikelihoodEvals += likelihoodEvaluations;
        previousLogLike = logLike;
        if(fabs(diff) < 1.0e-6) {
            break;
        }
        pass++;
        if(pass > aTree->newtAndBrentPowellOptPassLimit[0]) {
            printf("=============================================================\n");
            printf("p4_newtAndBrentPowellOpt()\n"); 
            printf("Pass limit (%i) exceeded without convergence.\n", 
                   aTree->newtAndBrentPowellOptPassLimit[0]);
            printf("The result may still be useful.\n");
            printf("You may want to optimize again, and perhaps boost the var.newtAndBrentPowellOptPassLimit.\n");
            break;
        }
    }
				

    logLike = p4_treeLogLike(aTree, 0);
    totalLikelihoodEvals++;
    if((0)) {
        printf("     p4_newtAndBrentPowellOpt().  logLike %f\n", logLike);
        printf("     %i total likelihood evaluations\n", totalLikelihoodEvals);
        //printf("Starting the opt again.  About to windUpParameters.\n");
    }

    //printf("p4_newtAndBrentPowellOpt().  %i full likelihood evaluations\n", totalLikelihoodEvals);

#if 0
    {
        int i;

        printf("p4_treeOptimize.c: newtAndBrentPowellOpt() after optimizing\n");
        printf("nPrams is %i\n", nPrams);
        for(i = 0; i < nPrams; i++) {
            printf("         parameters[%i] is %.12f\n", i, parameters[i]);
        }
        if(1) {
            printf("\n");
            for(i = 0; i < aTree->nNodes; i++) {
                if(aTree->nodes[i] != aTree->root) {
                    printf("         branchLength[%i] is %.12f\n", i, aTree->nodes[i]->brLen[0]);
                }
            }
        }
    }
#endif

    if(aBrent) {
        freeBrent(aBrent);
        aBrent = NULL;
    }
    if(parameters) {
        free(parameters);
        parameters = NULL;
    }

#if 1
    // this is a bit of a hack, to fix up some pathological data opts.
    if(aTree->nNodes < 15) {
        p4_allBrentPowellOptimize(aTree);
    }
#endif	

}



// I tried using the gsl 1D brent minimizer, but it wouldn't work if
// the minimum was at a bracket border.  That would happen a lot, eg
// pInvar=0.0.  The gsl brent minimizer requires a bracket a < m < b,
// and a = m < b will not work.

void p4_newtAnd1DBrent(p4_tree *aTree)
{
    double   guess1, guess2, lBound, uBound, pram;
    double   diff, logLikeBefore;
    int      iter, max_iter;
    double   downFactor;

    //printf("p4_newtAnd1DBrent() here.\n");
    p4_newtAround(aTree, 1.0, 10.0);
    p4_newtAround(aTree, 1.0e-1, 1.0);
    p4_newtAround(aTree, 1.0e-2, 0.1);
    p4_bracketTheMinimum(aTree, &guess1, &guess2, &lBound, &uBound, &pram, 0.1);
    aTree->logLike = -LocalMin(guess1, guess2, 1.0e-5, 1.0e-10, 
                               p4_minusLogLikeForBrent, &pram);
    //printf("a p4_newtAnd1DBrent().  Got logLike %f\n", aTree->logLike);
    logLikeBefore = aTree->logLike;
    iter = 0;
    max_iter = 20;
    while(1) {
        p4_newtAround(aTree, 1.0e-5, 1.0e-7);
        downFactor = 0.1 / ((double)(iter + 1));
        //printf("iter=%i, downFactor = %f\n", iter, downFactor);
        p4_bracketTheMinimum(aTree, &guess1, &guess2, &lBound, &uBound, &pram, downFactor);
        aTree->logLike = -LocalMin(guess1, guess2, 1.0e-5, 1.0e-10, 
                                   p4_minusLogLikeForBrent, &pram);
        diff = aTree->logLike - logLikeBefore;
        //printf("a p4_newtAnd1DBrent().  diff=%.8f\n", diff);
        if(diff < -1.0e-6) {
            printf("a p4_newtAnd1DBrent().  Likelihood got worse.  Bad!\n");
            exit(1);
        }
        else if(fabs(diff) < 1.0e-6) {
            break;
        }
        else if(iter > max_iter) {
            printf("a p4_newtAnd1DBrent().  Max iterations reached.  Not converged!\n");
            break;
        }
        iter++;
        logLikeBefore = aTree->logLike;
    }
		
}


void p4_bracketTheMinimum(p4_tree *aTree, double *guess1, double *guess2, double *lBound, double *uBound, double *pram, double downFactor)
{
    p4_windUpParameters(aTree, pram, lBound, uBound);
    //printf("p4_bracketTheMinimum(). *pram = %f, *lBound = %f, *uBound = %f\n", *pram, *lBound, *uBound);
    //printf("p4_bracketTheMinimum(). *pram = %f, downFactor= %f\n", *pram, downFactor);
    *guess2 = *pram;
    *guess1 = *guess2 * (1.0 - downFactor);
    if(*guess1 <= *lBound || *guess1 >= *uBound) *guess1 = *lBound;
    if(*guess2 <= *lBound || *guess2 >= *uBound) *guess2 = *uBound;
    //printf("  Before BracketMinimum, guess1 = %f and guess2 = %f\n", *guess1, *guess2);
    BracketMinimum(guess1, guess2, p4_minusLogLikeForBrent, *lBound, *uBound);
    if(*guess1 <= *lBound || *guess1 >= *uBound) *guess1 = *lBound;
    if(*guess2 <= *lBound || *guess2 >= *uBound) *guess2 = *uBound;
    //printf("  After BracketMinimum, guess1 = %f and guess2 = %f\n", *guess1, *guess2);
	
}

// The next bit is from ...  (Thanks, Dave!)
/*	minfunc.c
|
|	Routines for function minimization based on those of Brent (1973; "Algorithms for
|	Minimization Without Derivatives"; Prentice-Hall).  Brent's routines were translated from
|	Algol to C and cleaned up (goto elimination, provision for dynamic memory allocation, etc.).
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
*/

 
/*--------------------------------------------------------------------------------------------------
|
|	LocalMin
|
|	Find a local minimum of a function of one variable using Brent's (1973) method.
|
|	This function is reentrant.
*/

double LocalMin(double a, double b, double eps, double t, double (*f)(double *), double *px)
	//double		a, b;	/* on input, (a,b) must bracket a local minimum */
	//double		eps;	/* t and eps define tol = eps|x| + t, f is never evaluated at two points */
	//double		t;		/*   closer together than tol;  eps should be > sqrt(DBL_EPSILON)        */
	//MinimizeFxn	f;		/* function to minimize */
	//double		*px;	/* value of x when f(x) is minimal */
{
    double		e, m, p, q, r, x, tol, t2, u, v, w, fu, fv, fw, fx,
        d = 0.0;	/* shuts up bogus lint warning */
    int			iter;

    v = w = x = a + CGOLD*(b - a);
    e = 0.0;
    fv = fw = fx = (*f)(&x);
    //RETURN_IF_ABORTED
	

    /* main loop */

    for (iter = 0; iter < MAX_ITER; iter++)
    {
        m = 0.5*(a + b);
        tol = eps*fabs(x) + t;
        t2 = 2.0*tol;

        /* check stopping criterion */
        if (fabs(x - m) <= t2 - 0.5*(b - a))
            break;

        p = q = r = 0.0;
        if (fabs(e) > tol)
        {
            /* fit parabola (trial) */
            r = (x - w)*(fx - fv);
            q = (x - v)*(fx - fw);
            p = (x - v)*q - (x - w)*r;
            q = 2.0*(q - r);
            if (q > 0.0)
                p = -p;
            else
                q = -q;
            r = e;
            e = d;	/* lint complains about possible use of 'd' before being set, but it's not
                           a problem because e=0 first time through, so "fabs(e) > tol" test fails */
        }
			
        /* Take parabolic-interpolation or golden-section step (note that Brent's Algol procedure
           had an error, p > q*(a-x) is correct) */
        if ((fabs(p) < fabs(0.5*q*r)) && (p > q*(a - x)) && (p < q*(b - x)))
        {
            /* parabolic interpolation step */
            d = p/q;
            u = x + d;
            /* don't evaluate f too close to a or b */
            if ((u - a < t2) || (b - u < t2))
                d = (x < m) ? tol : -tol;
        }
        else
        {
            /* "golden section" step */
            e = ((x < m) ? b : a) - x;
            d = CGOLD*e;
        }

        /* don't evaluate f too close to x */
        if (fabs(d) >= tol)
            u = x + d;
        else if (d > 0.0)
            u = x + tol;
        else 
            u = x - tol;

        fu = (*f)(&u);
        //RETURN_IF_ABORTED
		
        /* update a, b, v, w, and x */
        if (fu <= fx)
        {
            if (fu <= fx)
            {
                if (u < x)
                    b = x;
                else
                    a = x;
                v = w;
                fv = fw;
                w = x;
                fw = fx;
                x = u;
                fx = fu;
            }
        }
        else
        {
            if (u < x)
                a = u;
            else
                b = u;
            if ((fu <= fw) || (w == x))
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if ((fu <= fv) || (v == x) || (v == w))
            {
                v = u;
                fv = fu;
            }
        }
    }
    *px = x;
    return fx;
}



/*--------------------------------------------------------------------------------------------------
|
|	BracketMinimum
|
|	Bracket a function minimum using method of Press et al.
|
|	This function is reentrant.
*/

//void BracketMinimum(double *pA, double *pB, MinimizeFxn func, double minAllowed, double maxAllowed)
void BracketMinimum(double *pA, double *pB, double (*func)(double []), double minAllowed, double maxAllowed)
{
    double		a, b, c, x, fa, fb, fc, fx, xlim, r, q, temp;
    int applyMinMax = 1;
	
    fa = (*func)(pA);
    a = *pA;
    fb = (*func)(pB);
    b = *pB;
    //printf("beginning BracketMinimum. fa=%f, fb=%f\n", fa, fb);
    if (fb > fa) {
        // swap a and b so we can go downhill in direction of a to b
        temp = a;
        a = b;
        b = temp;
        temp = fa;
        fa = fb;
        fb = temp;
    }
    c = b + GOLD*(b - a);		// first guess for c
    if(applyMinMax) {
        if(c < minAllowed) {
            //printf("BracketMinimum: new guess c (now %f) forced to the minAllowed, %f\n", 
            //	   c, minAllowed);
            c = minAllowed;
        }
        else if(c > maxAllowed) {
            //printf("BracketMinimum: new guess c (now %f) forced to the maxAllowed, %f\n", 
            //	   c, maxAllowed);
            c = maxAllowed;
        }
    }
    //printf("c = %f\n", c);
    fc = (*func)(&c);

    //printf("\nBracketMinimum: initial try = (%g,%g,%g) => (%g,%g,%g)\n", a, b, c, fa, fb, fc);

    while (fb > fc) {
        // use inverse parabolic interpolation to compute a new trial point x (= the abcissa
        //   which is the minimum of a parabola through f(a), f(b), and f(c))  
        r = (b - a)*(fb - fc);
        q = (b - c)*(fb - fa);
        x = b - ((b - c)*q - (b - a)*r) / (2.0*SIGN(MAX(fabs(q-r), TINY), q-r));
        //new from peter.  Is this ok?
        if(applyMinMax) {
            if(x < minAllowed) {
                //printf("BracketMinimum:a: new guess x (now %f) forced to the minAllowed, %f\n", 
                //	   x, minAllowed);
                x = minAllowed;
            }
            else if(x > maxAllowed) {
                //printf("BracketMinimum:a: new guess x (now %f) forced to the maxAllowed, %f\n", 
                //	   x, maxAllowed);
                x = maxAllowed;
            }
        }
        xlim = b + GLIMIT*(c - b);
        if ((b - x)*(x - c) > 0.0) {
            // x is between b and c
            fx = (*func)(&x);
            //printf("  1: x => %g (fx = %g)\n", x, fx);
            if (fx < fc) {
                // there's a minimum between b and c
                a = b;
                // lines below are apparently do-nothing statements 
                //fa = fb;
                //b = x;
                //fb = fx;
                //printf("  A: (%g,%g,%g) => (%g,%g,%g)\n", a, b, c, fa, fb, fc);
                break;
            }
            else if (fx > fb) {
                // there's a minimum between a and x 
                c = x;
                // line below is apparently a do-nothing statement 
                //fc = fx;
                //printf("  B: (%g,%g,%g) => (%g,%g,%g)\n", a, b, c, fa, fb, fc);
                break;
            }
            // parabolic fit failed; "magnify" using default magnification 
            x = c + GOLD*(c - b);
            //peter was here
            if(applyMinMax) {
                if(x < minAllowed) {
                    //printf("BracketMinimum:b: new guess x (now %f) forced to the minAllowed, %f\n", 
                    //	   x, minAllowed);
                    x = minAllowed;
                }
                else if(x > maxAllowed) {
                    //printf("BracketMinimum:b: new guess x (now %f) forced to the maxAllowed, %f\n", 
                    //	   x, maxAllowed);
                    x = maxAllowed;
                }
            }
            fx = (*func)(&x);
            //printf("  2: x => %g (fx = %g)\n", x, fx);
        }
        else if ((c - x)*(x - xlim) > 0.0) {
            // x from parabolic fit is between c and its allowed limit 
            fx = (*func)(&x);
            //printf("  3: x => %g (fx = %g)\n", x, fx);
            if (fx < fc) {
                b = c; fb = fc;
                c = x; fc = fx;
                x = c + GOLD*(c - b);
                //peter was here
                if(applyMinMax) {
                    if(x < minAllowed) {
                        //printf("BracketMinimum:c: guess x (now %f) forced to the minAllowed, %f\n", 
                        //	   x, minAllowed);
                        x = minAllowed;
                    }
                    else if(x > maxAllowed) {
                        //printf("BracketMinimum:c: guess x (now %f) forced to the maxAllowed, %f\n", 
                        //	   x, maxAllowed);
                        x = maxAllowed;
                    }
                }
                fx = (*func)(&x);
                //printf("  4: x => %g (fx = %g)\n", x, fx);
            }
        }
        else if ((x - xlim)*(xlim - c) > 0.0) {
            // limit parabolic fit to its maximum allowed value 
            x = xlim;
            fx = (*func)(&x);
            //printf("  5: x => %g (fx = %g)\n", x, fx);
        }
        else {
            // reject parabolic x; use default magnification 
            x = c + GOLD*(c - b);
            //peter was here
            if(applyMinMax) {
                if(x < minAllowed) {
                    //printf("BracketMinimum:d: new guess x (now %f) forced to the minAllowed, %f\n", 
                    //	   x, minAllowed);
                    x = minAllowed;
                }
                else if(x > maxAllowed) {
                    //printf("BracketMinimum:d: new guess x (now %f) forced to the maxAllowed, %f\n", 
                    //	   x, maxAllowed);
                    x = maxAllowed;
                }
            }
            fx = (*func)(&x);
            //printf("  6: x => %g (fx = %g)\n", x, fx);
        }
		
        // new bracket is (b,c,x) 
        a = b; fa = fb;
        b = c; fb = fc;
        c = x; fc = fx;		
        //printf("  C: (%g,%g,%g) => (%g,%g,%g)\n", a, b, c, fa, fb, fc);
    }
    if (a < c) {
        *pA = a;
        *pB = c;
    }
    else {
        *pA = c;
        *pB = a;
    }
    //printf("  BracketMinimum returning (%g,%g)\n", *pA, *pB);
}


