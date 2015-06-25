// PGF got this from Dave Swofford.  Thanks, Dave!  It is not GPL'd, obviously.

/*	gamma.h, now paupgamma.h
|
|	Prototypes for functions relating to gamma and other statistical distributions.
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
|	All rights reserved.
*/

#ifdef __cplusplus
extern "C" {
#endif

extern double	BinomP (int, int);
extern double	Chi2P (double, int);
extern double	GammaD (double, double, int *);
extern double	PaupLnGamma (double, int *);
extern double	NormP (double, double *);
extern double	PPChi2 (double, double, double, int *);
extern double	StudentTP (double, int);

#ifdef __cplusplus
}
#endif
