void checkpoint(char *location);
void checkpointOneInt(char *location, int theInt);
double ranDoubleUpToOne(void);
void setBigQFromRMatrixDotCharFreq(double **theBigQ, double **theRMatrix, double *charFreq, int dim);
void normalizeBigQ(double **theBigQ, double *charFreq, int dim);
int indexOfIntInArray(int theInt, int *theIntArray, int arrayLength);
//void *pdbmalloc(size_t size);
//void pdbfree(void *ptr);

double logOfSum(double *theLogs, int len);
double newtonRaftery94_eqn16(double *logLikes, int len, double harmMean, double delta, int verbose);



// I swiped this from bambe.  Thanks, guys!
#ifdef ALPHA
#define SAFE_EXP(x) ((x)<-200.0 ? 0.0 : exp(x))
#else
#define SAFE_EXP(x) (exp(x))
#endif
