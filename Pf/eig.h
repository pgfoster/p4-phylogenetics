eig *allocEig(int theDim, double **theMat);
void allocComplexStuff(eig *anEig);
void freeEig(eig *anEig);
int eigensystem(eig *anEig);
void matrixExpTimesBranchLength(eig *anEig, double branchLength, double **result);
void matrixLog(eig *anEig, double **result);
void matrixPower(eig *anEig, double pow, double **result);
void firstDerivativeOfMatrixExpTimesBranchLength(eig *anEig, double branchLength, double **result, double rate);
void secondDerivativeOfMatrixExpTimesBranchLength(eig *anEig, double branchLength, double **result, double rate);
//void firstDerivativeOfMatrixExpTimesBranchLengthL(eig *anEig, double branchLength, double **result, double rate);
//void secondDerivativeOfMatrixExpTimesBranchLengthL(eig *anEig, double branchLength, double **result, double rate);



