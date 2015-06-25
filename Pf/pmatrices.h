double **psdmatrix(int dim);	// mallocs a square double matrix
void free_psdmatrix(double **m); // frees the square double matrix
void dump_psdmatrix(double **m, int dim); // prints a square double matrix
void copy_psdmatrix(double **from, double **to, int dim); // copies

double **pdmatrix(int rows, int cols); // mallocs a rectangular matrix
void free_pdmatrix(double **m); // frees the rectangular double matrix
int **psimatrix(int dim);	// mallocs a square int matrix
void free_psimatrix(int **m); // frees the square int matrix
void dump_psimatrix(int **m, int dim); // prints a square int matrix
int **pimatrix(int rows, int cols);  // mallocs a rectangular int matrix
void free_pimatrix(int **m);  
double *pdvector(int dim);	// mallocs a double vector
void dump_pdvector(double *vec, int dim);
int *pivector(int dim);		// mallocs an int vector

void transpose_psdmatrix(double **m, int dim);
double *dotMultMatrixByVector(double **mat, double *vec, int dim);
void dotMultMatrixByMatrix(double **leftmat, double **rightmat, int dim, double **result);
void multMatrixColumnsByVector(double **leftmat, double *rightvec, int dim, double **result);

void diagonalFillMatrix(double **sqMatrix, int dim);
void zeroFillSquareMatrix(double **sqMatrix, int dim);
