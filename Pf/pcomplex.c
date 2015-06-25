// Where did I steal this from?  Yang I think.  With lots of changes...
// And there is some code that I modified from Dave Swofford, below.
// Not GPL'd

#include "pftypes.h"
#include "pcomplex.h"
#include <math.h>
#include	<stdio.h>
#include <stdlib.h>	// for malloc and free
#include "pmatrices.h"
#include "util.h"


#define PIOVER2 1.57079632679489662

pcomplex Complex(double a, double b)
{
	pcomplex c;
	c.re = a;
	c.im = b;
	return c;
}

pcomplex Cadd(pcomplex a, pcomplex b)
{
	pcomplex c;
	c.re = a.re + b.re;
	c.im = a.im + b.im;
	return c;
}

pcomplex Csub(pcomplex a, pcomplex b)
{
	pcomplex c;
	c.re = a.re - b.re;
	c.im = a.im - b.im;
	return c;
}

pcomplex Cmul(pcomplex a, pcomplex b)
{
	pcomplex c;
	c.re = a.re * b.re - a.im * b.im;
	c.im = a.im * b.re + a.re * b.im;
	return c;
}

pcomplex Conj(pcomplex a)
{
	pcomplex c;
	c.re = a.re;
	c.im = -a.im;
	return c;
}

pcomplex Cdiv(pcomplex a, pcomplex b)
{
	pcomplex c;
	double r, den;
	if( fabs(b.re) >= fabs(b.im) ) {
		r = b.im / b.re;
		den = b.re + r * b.im;
		c.re = (a.re + r * a.im) / den;
		c.im = (a.im - r * a.re) / den;
	} else {
		r = b.re / b.im;
		den = b.im + r * b.re;
		c.re = (a.re * r + a.im) / den;
		c.im = (a.im * r - a.re) / den;
	}
	return c;
}

double Cabs(pcomplex a)
{
	double x, y, ans, temp;
	x = fabs(a.re);
	y = fabs(a.im);
	if(x == 0.0) 
		ans = y;
	else if(y == 0.0)
		ans = x;
	else if(x > y) {
		temp = y / x;
		ans = x * sqrt(1.0 + temp * temp);
	} else {
		temp = x / y;
		ans = y * sqrt(1.0 + temp * temp);
	}
	return ans;
}

pcomplex Csqrt(pcomplex a)
{
	pcomplex c;
	double x, y, w, r;
	if( (a.re == 0.0) && (a.im == 0.0) ) {
		c.re = 0.0;
		c.im = 0.0;
		return c;
	} else {
		x = fabs(a.re);
		y = fabs(a.im);
		if(x >= y) {
			r = y / x;
			w = sqrt(x) * sqrt(0.5 * (1.0 + sqrt(1.0 + r * r)));
		} else {
			r = x / y;
			w = sqrt(y) * sqrt(0.5 * (r + sqrt(1.0 + r * r)));
		}
		if(a.re >= 0.0) {
			c.re = w;
			c.im = a.im / (2.0 * w);
		} else {
			c.im = (a.im >= 0.0) ?  w : -w;
			c.re = a.im / (2.0 * c.im);
		}
		return c;
	}
}
			
pcomplex RCmul(double a, pcomplex b)
{
	pcomplex c;
	c.re = a * b.re;
	c.im = a * b.im;
	return c;
}

pcomplex Cexp(pcomplex a)
{
   pcomplex c;
   c.re = SAFE_EXP(a.re);
   if (fabs(a.im)==0) c.im = 0; 
   else  { c.im = c.re*sin(a.im); c.re*=cos(a.im); }
   return (c);
}



pcomplex	Clog(pcomplex a)
{
	pcomplex	c;
	c.re = log(Cabs(a));
	if(a.re == 0.0) {
		c.im = PIOVER2;
	} else {
		c.im = atan2(a.im, a.re);
	}
	return c;
}

pcomplex **pscmatrix(int dim)  // a square matrix
{
	int i;
	pcomplex **m;
	// allocate pointers to rows
	m=(pcomplex **) malloc((size_t)((dim)*sizeof(pcomplex*)));
	if (!m) {
		printf("allocation error in pscmatrix 1.\n");
		exit(1);
	}
	// allocate rows and set pointers to them 
	m[0]=(pcomplex *) malloc((size_t)((dim*dim)*sizeof(pcomplex)));
	if (!m[0]) {
		printf("allocation error in pscmatrix 2.\n");
		exit(1);
	}
	for(i=1;i<dim;i++) {
		m[i]=m[i-1]+dim;
	}
	// return pointer to array of pointers to rows
	return m;
}

void free_pscmatrix(pcomplex **m)
// free a pcomplex matrix allocated by pscmatrix()
// this works - I don't know why its a char *, maybe it should be a void *?

{
	free((char *) (m[0]));
	free((char *) (m));
}

void dump_pscmatrix(pcomplex **m, int dim)
{
	int	row, col;
	
	printf("{");
	for(row = 0; row < (dim - 1); row++) {
		printf("{");
		for(col = 0; col < (dim - 1); col++) {
			printf("%f + %f I, ", m[row][col].re, m[row][col].im);
			if(col == 1) printf("\n    ");
		}
		printf("%f + %f I},\n", 
			m[row][dim - 1].re, m[row][dim - 1].im);
	}
	printf("{");
	for(col = 0; col < (dim - 1); col++) {
		printf("%f + %f I, ", m[dim - 1][col].re, m[dim - 1][col].im);
		if(col == 1) printf("\n    ");
	}
	printf("%f + %f I}}", 
		m[dim - 1][dim - 1].re, m[dim - 1][dim - 1].im);
	printf("\n");
}

void copy_pscmatrix(pcomplex **from, pcomplex **to, int dim)
{
	int row, col;

	for(row = 0; row < dim; row++) {
		for(col = 0; col < dim; col++) {
			to[row][col] = from[row][col];
		}
	}

}

void dump_complexVector(pcomplex *vec, int dim)
{
	int	i;
	
	printf("{");
	for(i = 0; i < (dim - 1); i++) {
		printf("%f + %f I, ", vec[i].re, vec[i].im);
		if(i == 1) printf("\n    ");
	}
	printf("%f + %f I}\n", vec[dim - 1].re, vec[dim - 1].im);
}



// I used Dave Swoffords routines as a starting point, and complexified them.

/*----------------------------------------------------------------------------
|
|	ComplexInvertMatrix
|
|	Invert matrix 'a' using LU-decomposition technique, 
|   storing inverse in 'a_inv'.  Matrix 'a'
|	is destroyed.  Returns 1 if matrix is singular, 0 otherwise.
*/

int ComplexInvertMatrix(pcomplex **a, int n, double *dwork, int *indx,
                        pcomplex **a_inv, pcomplex *col)
//	pcomplex		**a;		matrix represented as vector of row pointers
//	int			n;			order of matrix
//	double		*dwork;		work vector of size n
//	int			*indx;		work vector of size n
//	pcomplex		**a_inv;	inverse of input matrix a (matrix a is destroyed)
//  pcomplex		*col;		for the second part, involving back subst

{
	int			rc, i, j;
	
	rc = ComplexLUDecompose(a, n, dwork, indx, (double *)NULL);
	//printf("ComplexLUDecompose() returned %i\n", rc);

	if (rc == 0) {
		for (j = 0; j < n; j++) {
			for (i = 0; i < n; i++)
				col[i] = Complex(0.0, 0.0);
			col[j] = Complex(1.0, 0.0);
			ComplexLUBackSubst(a, n, indx, col);
			for (i = 0; i < n; i++)
				a_inv[i][j] = col[i];
		}
	}
	return rc;
}

/*-------------------------------------------------------------------------
|
|	ComplexLUDecompose
|
|	Replace matrix 'a' with its LU-decomposition.  
|   Returns 1 if matrix is singular, 0
|	otherwise.
*/

int ComplexLUDecompose(pcomplex **a, int n, double *vv, int *indx, double *pd)
//	pcomplex		**a;	the matrix whose LU-decomposition is wanted
//	int			n;		order of a
//	double		*vv;	work vector of size n (stores implicit 
//							scaling of each row)
//	int			*indx;	=> row permutation according to partial 
//							pivoting sequence
//	double		*pd;	=> 1 if number of row interchanges was even, 
//							-1 if odd (NULL OK)
{
	int			i, imax, j, k;
	double		big, dum, temp, d;
	pcomplex		sum, cdum;

	d = 1.0;
	imax = 0; // only to shut the compiler up.

	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++) {
			if ((temp = Cabs(a[i][j])) > big)
				big = temp;
		}
		if (big == 0.0) {
                    printf("singular matrix in routine ComplexLUDecompose\n");
                    return 1;
		}
		vv[i] = 1.0 / big;
	}
	
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = a[i][j];
			for (k = 0; k < i; k++) 
				sum = Csub(sum, Cmul(a[i][k], a[k][j]));
			a[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum = Csub(sum, Cmul(a[i][k], a[k][j]));
			a[i][j] = sum;
			dum = vv[i] * Cabs(sum);
			if (dum >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				cdum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = cdum;
			}	
			d = -d;
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (a[j][j].re == 0.0 && a[j][j].im == 0.0)
			a[j][j] = Complex(1.0e-20, 1.0e-20);
		if (j != n - 1){
			cdum = Cdiv(Complex(1.0, 0.0), a[j][j]);
			for (i = j + 1; i < n; i++)
				a[i][j] = Cmul(a[i][j], cdum);
		}
	}

	if (pd != NULL)
		*pd = d;
	return 0;
}

/*----------------------------------------------------------------------
|
|	ComplexLUBackSubst
|
|	Perform back-substition into LU-decomposed matrix in order 
|	to obtain inverse.
*/

void ComplexLUBackSubst(pcomplex **a, int n, int *indx, pcomplex *b)
{
	int			i, ip, j,
				ii = -1;
	pcomplex		sum;

	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0) {
			for (j = ii; j <= i - 1; j++)
				sum = Csub(sum, Cmul(a[i][j], b[j]));
		} else if ((sum.re != 0.0) || (sum.im != 0.0))
			ii = i;
		b[i] = sum;
	}
	
	for (i = n - 1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j < n; j++)
			sum = Csub(sum, Cmul(a[i][j], b[j]));
		b[i] = Cdiv(sum, a[i][i]);
	}
}


