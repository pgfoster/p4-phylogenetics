#include <stdio.h>
#include <math.h>
#include <float.h>  // for DBL_EPSILON
#include <stdlib.h>
//#include "random.h"
//#include <libc.h>   // for srandom
#include "pmatrices.h"
#include "util.h"

//typedef double (*MinimizeFxn) (double *);

typedef struct {
	int	nl;
    double		dmin;
	double		ldt;
	double		qf1;
	double		qd0, qd1;
    double		m2, m4;	// constants
    double		eps2;	// a constant
    double		toler, htol;	// global versions of PrAxis arguments
    double		(*fxn)(double []);	// global versions of PrAxis arguments
    int		n;		// global versions of PrAxis arguments
    double		*gx;		// global versions of PrAxis arguments
    double		gfx;
    double		**vv;
    double		*d, *q0, *q1;
    double		*xnew;
    double		vsmall;
    double		large;
    double		vlarge;
    double		*y;
    double		*z;

} brent;


brent *newBrent(int dim);
void freeBrent(brent *aBrent);
double praxis(brent *aBrent, double tol, double h, int dim, double *xx, double (*theFunction)(double []));
void minFit(brent *aBrent, double eps, double tol, double **ab, double *q, double *e);
void sortDV(brent *aBrent);
void lineMin(brent *aBrent, int j, int nits, double *pd2, double *px1, double f1, int fk);
double fLin(brent *aBrent, int j, double lambda);
void qUad(brent *aBrent);

#ifndef FALSE
#	define FALSE	0
#	define TRUE		1
#endif

#define MAX_ITER	50
#define CGOLD		0.3819660113	// (3 - sqrt(5))/2 
#define GOLD		1.618034
#if 0
#	define GLIMIT		100.0	/* DLS 07mar98: this was the original setting (see below) */
#else
	/* DLS 07mar98: the original GLIMIT setting of 100.0 was sometimes allowing too great an
       extrapolation, causing problems with gamma shape parameter estimation in PAUP*.  Setting
	   to a smaller value shouldn't hurt anything; in the worst case it will just take a little
	   longer to find the bracket */
#	define GLIMIT		10.0
#endif

#define TINY		1.0e-20
#define SIGN(a, b)	((b) > 0.0 ? fabs(a) : -fabs(a))
#define FOREVER ;;
#if !defined(MAX)
#	define MAX(x, y)	((x) >= (y)	? (x) : (y))
#endif
