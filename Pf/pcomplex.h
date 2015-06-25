// Not GPL'd

pcomplex Complex(double a, double b);
pcomplex Cadd(pcomplex a, pcomplex b);
pcomplex Csub(pcomplex a, pcomplex b);
pcomplex Cmul(pcomplex a, pcomplex b);
pcomplex Conj(pcomplex a);
// Note that DLS has a CDiv function in linalg.  
// (With the slightly different spelling.)
pcomplex Cdiv(pcomplex a, pcomplex b);
double Cabs(pcomplex a);
pcomplex Csqrt(pcomplex a);
pcomplex RCmul(double a, pcomplex b);
pcomplex Cexp(pcomplex a);
pcomplex	Clog(pcomplex a);

pcomplex **pscmatrix(int dim);	// mallocs a square pcomplex matrix
void free_pscmatrix(pcomplex **m); // frees the square pcomplex matrix
void dump_pscmatrix(pcomplex **m, int dim); // prints a square pcomplex matrix
void copy_pscmatrix(pcomplex **from, pcomplex **to, int dim); // copies
void dump_complexVector(pcomplex *vec, int dim);

int ComplexInvertMatrix(pcomplex **a, int n, double *dwork, int *indx, 
										pcomplex **a_inv, pcomplex *col);
int ComplexLUDecompose(pcomplex **a, int n, double *vv, int *indx, double *pd);
void ComplexLUBackSubst(pcomplex **a, int n, int *indx, pcomplex *b);


