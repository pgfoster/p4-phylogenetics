// PGF got this from Dave Swofford.  Thanks Dave!  Not GPL'd.

/*	linalg.h
|
|	Prototypes for matrix-inversion and eigensystem functions
|
|	Copyright (c) 1995 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
|
|	NOTE: if ANSI function prototypes are not supported, define NO_PROTOTYPES
|		  before including this file.
*/

#ifdef __cplusplus
extern "C" {
#endif

#if !defined(PROTO_)
#	ifdef NO_PROTOTYPES
#		define PROTO_(a) ()
#	else
#		define PROTO_(a) a
#	endif
#endif

#define RC_COMPLEX_EVAL 2	/* code that complex eigenvalue obtained */

extern int	InvertMatrix PROTO_((double **, int, double *, int *, double **));
extern int	LUDecompose PROTO_((double **, int, double *, int *, double *));
extern int	EigenRealGeneral PROTO_((int, double **, double *, double *, double **, int *, double *));

#ifdef __cplusplus
}
#endif
