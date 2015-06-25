// PGF got this from Dave Swofford.  Thanks, Dave!  It is not GPL'd,
// but it is gopyright by Dave.  It is ok with Dave if it is used, as
// long as acknowledgement to Dave is given.

/*	gamma.c
|
|	This file contains miscellaneous routines related to statistical distributions including gamma
|	and beta functions.
|
|	Copyright (c) 1998 by David L. Swofford, Smithsonian Institution.
*/

//#define DOES_MATH
//#include "dls.h"
//#include "paup.h"
#define ENTRY
#define LOCAL static
#include "dlsgamma.h"
//#include "like.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "util.h"


static double	BetaCF (double, double, double);
static double	BetaI (double, double, double);
static double	PPND16 (double, int *);

#if !defined(MIN)
#	define MIN(x, y)	((x) <= (y)	? (x) : (y))
#endif

/*--------------------------------------------------------------------------------------------------
|
|	GammaD
|
|	Compute incomplete gamma integral.
|
|	Algorithm AS239  Appl. Stat. vol. 37, no. 3 (1988).
|
| 	Auxiliary functions required: ALOGAM = logarithm of the gamma function, and ALNORM = algorithm
|	AS66 (DLS: I instead use "NORMP" from Alan Miller of CSIRO, which is supposedly more accurate).
*/

#define C_B2	0.33333333333333331
#define PLIMIT	1000.0

ENTRY double GammaD(x, p, pErrorCode)
	double		x, p;
	int			*pErrorCode;
{
	double		a, b, c, an, rn, pn1, pn2, pn3, pn4, pn5, pn6, arg, tmp, result;

	/* validate X and P */
    if ((p <= 0.0) || (x < 0.0))
    	{
    	if (pErrorCode != NULL)
			*pErrorCode = 1;
		return 0.0;
	    }

	if (pErrorCode != NULL)
		*pErrorCode = 0;

    if (x == 0.0)
		return 0.0;

	/* use a normal approximation if P > PLIMIT */
    if (p > PLIMIT)
    	{
		pn1 = sqrt(p)*3.0*(pow(x/p, C_B2) + 1.0/(p*9.0) - 1.0);
		return NormP(pn1, (double *)NULL);
		}

	/* if X is extremely large compared to P then set GAMMAD = 1 */
    if (x > 1e8)
		return 1.0;

    if ((x <= 1.0) || (x < p))
    	{
		/* Use Pearson's series expansion. (Note that P is not large enough to force overflow in
		   PaupLnGamma.)  No need to test error return from PaupLnGamma since P > 0. */
		arg = p*log(x) - x - PaupLnGamma(p + 1.0, (int *)NULL);
		c = 1.0;
		result = 1.0;
		a = p;
		do	{
			a += 1.0;
			c = c * x/a;
			result += c;
			}
			while (c > 1e-14);
		arg += log(result);
		result = 0.0;
		if (arg >= -88.0)
		    result = SAFE_EXP(arg);
    	}
    else
	    {
		/* use a continued fraction expansion */
		arg = p*log(x) - x - PaupLnGamma(p, (int *)NULL);
		a = 1.0 - p;
		b = a + x + 1.0;
		c = 0.0;
		pn1 = 1.0;
		pn2 = x;
		pn3 = x + 1.0;
		pn4 = x*b;
		result = pn3/pn4;
		for (;;)
			{
			a += 1.0;
			b += 2.0;
			c += 1.0;
			an = a*c;
			pn5 = b*pn3 - an*pn1;
			pn6 = b*pn4 - an*pn2;
			if (fabs(pn6) > 0.0)
				{
			    rn = pn5/pn6;
			    tmp = rn*1e-14;
			    if (fabs(result - rn) <= MIN(1e-14, tmp))
			    	break;
			    result = rn;
				}
		
			pn1 = pn3;
			pn2 = pn4;
			pn3 = pn5;
			pn4 = pn6;
			if (fabs(pn5) >= 1e37)
				{
				/* re-scale terms in continued fraction if terms are large */
			    pn1 /= 1e37;
			    pn2 /= 1e37;
			    pn3 /= 1e37;
			    pn4 /= 1e37;
				}
			}
		arg += log(result);
		result = 1.0;
		if (arg >= -88.0)
		    result = 1.0 - SAFE_EXP(arg);
	    }

    return result;
}

/*--------------------------------------------------------------------------------------------------
|
|	PPChi2
|
|	Get percentage points of the chi-squared distribution ("PPCHI2").
|
|	Algorithm used is:
|
|	Best, D. J. and D. E. Roberts.  1975.  The percentage points of the Chi2 distribution. 
|	Applied Statistics 24:385-388.  (Algorithm AS91)
|
|	To evaluate the percentage points of the chi-squared probability distribution function:
|
|		p must lie in the range 0.000002 to 0.999998,
|		v must be positive,
|		g must be supplied and should be equal to ln(gamma(v/2.0))
|
|	This routine incorporates the suggested changes in AS R85 (vol. 40(1), pp. 233-235, 1991) which
|	should eliminate the need for the limited range for p above, though these limits have not been
|	removed from the routine.  DLS: I removed them.
|
|	If IFAULT = 4 is returned, the result is probably as accurate as the machine will allow.
|
|	Auxiliary routines required: PPND = AS 111 (or AS 241) and GAMMAD = AS 239.
*/

#define AA	0.6931471806
#define E 	5e-7
#define C1	0.01
#define C2	0.222222
#define C3	0.32
#define C4	0.4
#define C5	1.24
#define C6	2.2
#define C7	4.67
#define C8	6.66
#define C9	6.73
#define C10	13.32
#define C11	60.0
#define C12	70.0
#define C13	84.0
#define C14	105.0
#define C15	120.0
#define C16	127.0
#define C17 140.0
#define C18	1175.0
#define C19	210.0
#define C20	252.0
#define C21	2264.0
#define C22	294.0
#define C23	346.0
#define C24	420.0
#define C25	462.0
#define C26	606.0
#define C27	672.0
#define C28	707.0
#define C29	735.0
#define C30	889.0
#define C31	932.0
#define C32	966.0
#define C33	1141.0
#define C34	1182.0
#define C35	1278.0
#define C36	1740.0
#define C37	2520.0
#define C38	5040.0

ENTRY double PPChi2(p, v, g, ifault)
	double		p, v, g;
	int			*ifault;
{
    int			i, if1;
    double		a, b, c, q, t, x, p1, p2, s1, s2, s3, s4, s5, s6, ch, xx, tmp;

	/* test arguments and initialize */
    *ifault = 1;
    if ((p < 0.0) || (p >= 1.0)) /* DLS: I allow any p in [0,1) rather than [0.000002,0.999998] */
		return -1.0;			 /*      as suggested by comment above                          */
    *ifault = 2;
    if (v <= 0.0)
		return -1.0;

    *ifault = 0;
    xx = 0.5*v;
    c = xx - 1.0;

    if (v < -C5*log(p))
		{
		/* starting approximation for small chi-squared */
	    ch = pow(p*xx*SAFE_EXP(g + xx*AA), 1.0/xx);
	    if (ch < E)
			return ch;
	    }
	else
    	{
	    if (v > C3)
	    	{
		    x = PPND16(p, &if1);	/* DLS: I use AS 241 here rather than AS 111 */
		
			/* starting approximation using Wilson and Hilferty estimate */
		    p1 = C2/v;
		    tmp = x*sqrt(p1) + 1.0 - p1;
		    ch = v*tmp*tmp*tmp;
		
			/* starting approximation for p tending to 1 */
		    if (ch > C6*v + 6.0)
				ch = -2.0*(log(1.0 - p) - c*log(0.5*ch) + g);
			}
		else
			{
			/* starting approximation for v less than or equal to 0.32 */
		    ch = C4;
		    a = log(1.0 - p);
			do	{
			    q = ch;
			    p1 = 1.0 + ch*(C7 + ch);
			    p2 = ch*(C9 + ch*(C8 + ch));
			    t = -0.5 + (C7 + 2.0*ch)/p1 - (C9 + ch*(C10 + 3.0*ch))/p2;
			    ch -= (1.0 - SAFE_EXP(a + g + 0.5*ch + c*AA)*p2/p1)/t;
			    }
			    while (fabs(q/ch - 1.0) > C1);
		    }
		}

	/* call to algorithm AS 239 and calculation of seven-term Taylor series */
    for (i = 1; i <= 20; ++i)
    	{
		q = ch;
		p1 = 0.5*ch;
		p2 = p - GammaD(p1, xx, &if1);
		if (if1 == 0)
			{
			t = p2*SAFE_EXP(xx*AA + g + p1 - c*log(ch));
			b = t/ch;
			a = 0.5*t - b*c;
			s1 = (C19 + a*(C17 + a*(C14 + a*(C13 + a*(C12 + C11*a)))))/C24;
			s2 = (C24 + a*(C29 + a*(C32 + a*(C33 + C35*a))))/C37;
			s3 = (C19 + a*(C25 + a*(C28 + C31*a)))/C37;
			s4 = (C20 + a*(C27 + C34*a) + c*(C22 + a*(C30 + C36*a)))/C38;
			s5 = (C13 + C21*a + c*(C18 + C26*a))/C37;
			s6 = (C15 + c*(C23 + C16*c))/C38;
			ch += t*(1.0 + 0.5*t*s1 - b*c*(s1 - b*(s2 - b*(s3 - b*(s4 - b*(s5 - b*s6))))));
			if (fabs(q/ch - 1.0) > E)
    			return ch;
			}
		else
			{
			*ifault = 3;
			return -1.0;
			}
	    }

	*ifault = 4;
    return ch;
}

/*--------------------------------------------------------------------------------------------------
|
|	PaupLnGamma
|
|	Uses Lanczos-type approximation to ln(gamma) for z > 0.
|
|	Reference:
|		Lanczos, C. 1964. A precision approximation of the gamma function. J. SIAM Numer. Anal. B
|		1:86-96. 
|
|	Accuracy: About 14 significant digits except for small regions in the vicinity of 1 and 2.
|
|	Originally coded in Fortran by Alan Miller (CSIRO Division of Mathematics & Statistics)
|	("latest revision - 17 April 1988").  Translated to C (and modified slightly) by David Swofford
|	(Smithsonian Institution).
*/

#define LN_SQRT_PI	0.9189385332046727
#define A0			0.9999999999995183

ENTRY double PaupLnGamma(z, pErrCode)
	double		z;				/* must be positive */
	int			*pErrCode;		/* can be NULL if caller knows that z is OK on input */
{
	int			j, ier;
	double		lngamma, tmp;

    static double a[8]
      = {  676.5203681218835,
	     -1259.139216722289,
	       771.3234287757674,
	      -176.6150291498386,
	        12.50734324009056,
	        -0.1385710331296526,
	         9.934937113930748e-6,
	         1.659470187408462e-7
	    };

	ier = 0;
	lngamma = 0.0;
    if (z <= 0.0)
		ier = 1;
	else
		{
	    tmp = z + 7.0;
	    for (j = 7; j >= 0; --j)
	    	{
			lngamma += a[j]/tmp;
			tmp -= 1.0;
			}
	    lngamma += A0;
		lngamma = log(lngamma) + LN_SQRT_PI - (z + 6.5) + (z - 0.5)*log(z + 6.5);
		}
	
	if (pErrCode != NULL)
		*pErrCode = ier;
	return lngamma;
}

/*--------------------------------------------------------------------------------------------------
|
|	PPND16
|
|
|	Produces the normal deviate Z corresponding to a given lower tail area of P; Z is accurate to
|	about 1 part in 10**16.
|
|	Algorithm AS241:  Appl. Stat. vol. 37, no. 3 (1988).
*/

LOCAL double PPND16(p, ifault)
	double		p;
	int			*ifault;
{
	double		result, q, r;

    *ifault = 0;
    q = p - 0.5;
    if (fabs(q) <= 0.425)
    	{
		r = 0.180625 - q*q;
		result = q*(((((((r*2509.0809287301226727  + 33430.575583588128105)*r
		                + 67265.770927008700853)*r + 45921.953931549871457)*r
		                + 13731.693765509461125)*r + 1971.5909503065514427)*r
		                + 133.14166789178437745)*r + 3.387132872796366608)
		         / (((((((r*5226.495278852854561   + 28729.085735721942674)*r
		                + 39307.89580009271061)*r  + 21213.794301586595867)*r
		                + 5394.1960214247511077)*r + 687.1870074920579083)*r
		                + 42.313330701600911252)*r + 1.0);
    	}
    else
    	{
		if (q < 0.0)
		    r = p;
		else
		    r = 1.0 - p;

		if (r <= 0.0)
			{
		    *ifault = 1;
		    return 0.0;
			}
		r = sqrt(-log(r));
		if (r <= 5.0)
			{
		    r += -1.6;
		    result = (((((((r*7.7454501427834140764e-4 + 0.0227238449892691845833)*r
		                   + 0.24178072517745061177)*r + 1.27045825245236838258)*r
		                   + 3.64784832476320460504)*r + 5.7694972214606914055)*r 
		                   + 4.6303378461565452959)*r  + 1.42343711074968357734)
		          / (((((((r*1.05075007164441684324e-9 + 5.475938084995344946e-4)*r
		                 + 0.0151986665636164571966)*r + 0.14810397642748007459)*r
		                 + 0.68976733498510000455)*r   + 1.6763848301838038494)*r
		                 + 2.05319162663775882187)*r   + 1.0);
			}
		else
			{
		    r += -5.0;
		    result = (((((((r*2.01033439929228813265e-7 + 2.71155556874348757815e-5)*r
		                  + 0.0012426609473880784386)*r + 0.026532189526576123093)*r
		                  + 0.29656057182850489123)*r   + 1.7848265399172913358)*r
		                  + 5.4637849111641143699) *r   + 6.6579046435011037772)
		        / (((((((r*2.04426310338993978564e-15   + 1.4215117583164458887e-7)*r
		                + 1.8463183175100546818e-5)*r   + 7.868691311456132591e-4)*r
		                + 0.0148753612908506148525)*r   + 0.13692988092273580531)*r
		                + 0.59983220655588793769)* r    + 1.0);
			}
		if (q < 0.0)
		    result = -result;
	    }
    return result;
}

/*--------------------------------------------------------------------------------------------------
|
|	NormP
|
|	Calculate the tail area under the normal curve (accurate to 1.e-15).
|
|	Algorithm: based upon algorithm 5666 for the error function, from:
|		Hart, J.F. et al, 'Computer Approximations', Wiley 1968
|
|	Input:
|	 	Z = no. of standard deviations from the mean
|
|	Output:
|	 	P, Q = probabilities to the left & right of Z.   P + Q = 1.
|		PDF = the probability density.
|		(DLS change: I only return P.  Caller can calculate Q = 1 - P if it wants it.  My version
|		returns P as function result.)
|
|	Originally coded in Fortran by Alan Miller (CSIRO Division of Mathematics & Statistics)
|	("latest revision - 30 March 1986").  Translated to C (and modified slightly) by David Swofford
|	(Smithsonian Institution).
*/

#define P0		220.2068679123761
#define P1		221.2135961699311
#define P2		112.0792914978709
#define P3		33.912866078383
#define P4		6.37396220353165
#define P5		0.7003830644436881
#define P6		0.03526249659989109
#define Q0		440.4137358247522
#define Q1		793.8265125199484
#define Q2		637.3336333788311
#define Q3		296.5642487796737
#define Q4		86.78073220294608
#define Q5		16.06417757920695
#define Q6		1.755667163182642
#define Q7		0.08838834764831844
#define CUTOFF	7.071				/* = 10/sqrt(2) */
#define SQRT2PI	2.506628274631001

ENTRY double NormP(z, pPdf)
	double		z;
	double		*pPdf;		/* NULL OK if caller doesn't want it */
{
	double		p, zabs, expntl, pdf;

    zabs = fabs(z);
    if (zabs > 37.0)
    	{
		pdf = 0.0;
		p = (z > 0.0) ? 1.0 : 0.0;
    	}
    else
    	{
	    expntl = SAFE_EXP(-0.5*zabs*zabs);
	    pdf = expntl/SQRT2PI;
	
	    if (zabs < CUTOFF)
			p = expntl*((((((P6*zabs + P5)*zabs + P4)*zabs + P3)*zabs  + P2)*zabs + P1)*zabs + P0) /
			 (((((((Q7*zabs + Q6)*zabs + Q5)*zabs + Q4)*zabs + Q3)*zabs + Q2)*zabs + Q1)*zabs + Q0);
	    else
			p = pdf/(zabs + 1.0/(zabs + 2.0/(zabs + 3.0/(zabs + 4.0/(zabs + 0.65)))));
	
	    if (z > 0.0)
			p = 1.0 - p;
	    }
	
	if (pPdf != NULL)
		*pPdf = pdf;
	return p;
}

/*--------------------------------------------------------------------------------------------------
|
|	Chi2P
|
|	Return P-value for specified chi-squared value and degrees-of-freedom.
*/

ENTRY double Chi2P(chiSq, df)
	double		chiSq;
	int			df;
{
	if (chiSq <= 0.0)
		return 1.0;
	else
		return 1.0 - GammaD(0.5*chiSq, 0.5*df, (int *)NULL);
}

/*--------------------------------------------------------------------------------------------------
|
|	BinomP
|
|	Perform a binomial test of the two-tailed hypothesis H0: p = 0.5 and return P-value.
*/

ENTRY double BinomP(n, n1)
	int			n,		/* total number of counts */
				n1;		/* number of counts in first category */
{
	int			n2;
	double		prob;

	n2 = n - n1;
	if (n1 == n2)
		return 1.0;
	else
		{
		if (n1 > n2)
			prob = BetaI((double)n1, (double)(n - n1 + 1), 0.5);
		else
			prob = BetaI((double)n2, (double)(n - n2 + 1), 0.5);
		
		return prob * 2.0;		/* two-tailed test; symmetric when p = 0.5 */
		}
}

/*--------------------------------------------------------------------------------------------------
|
|	StudentTP
|
|	Return P-value (two-tailed) for Student's t distribution.
*/

ENTRY double StudentTP(t, df)
	double		t;
	int			df;
{
	return BetaI(0.5*df, 0.5, df/(t*t + df));
}

/*--------------------------------------------------------------------------------------------------
|
|	BetaI
|
|	Evaluate incomplete beta function.
|
|	This code is based upon the "betai" routine from Press et al. (1988), "Numerical Recipes in C",
|	Cambridge University Press.
*/

LOCAL double BetaI(alpha, beta, x)
	double		alpha, beta, x;
{
	double		bt;

	if ((x == 0.0) || (x == 1.0))
		bt = 0.0;
	else
		bt = SAFE_EXP(PaupLnGamma(alpha + beta, (int *)NULL) - PaupLnGamma(alpha, (int *)NULL)
		     - PaupLnGamma(beta, (int *)NULL) + (alpha * log(x)) + (beta * log(1.0 - x)));
	if (x < ((alpha + 1.0) / (alpha + beta + 2.0)))
		return bt * BetaCF(alpha, beta, x) / alpha;
	else
		return 1.0 - (bt * BetaCF(beta, alpha, 1.0 - x) / beta);
}

/*--------------------------------------------------------------------------------------------------
|
|	BetaCF
|
|	Evaluate continued fraction for incomplete beta function (used by BetaI).
|
|	This code is based upon the "betacf" routine from Press et al. (1988), "Numerical Recipes in C",
|	Cambridge University Press.
*/

#define MAX_ITER	100
#define EPSIL		1.0e-7

LOCAL double BetaCF(alpha, beta, x)
	double		alpha, beta, x;
{
	int			m;
	double		qap, qam, qab, em, tem, d, bz, bp, bpp, ap, app, aold,
				bm = 1.0,
				az = 1.0,
				am = 1.0;

	qab = alpha + beta;
	qap = alpha + 1.0;
	qam = alpha - 1.0;
	bz = 1.0 - (qab * x / qap);
	for (m = 1; m <= MAX_ITER; m++)
		{
		em = (double)m;
		tem = em + em;
		d = em * (beta - em) * x / ((qam + tem) * (alpha + tem));
		ap = az + d*am;
		bp = bz + d*bm;
		d = -(alpha + em) * (qab + em) * x / ((qap + tem) * (alpha + tem));
		app = ap + d*az;
		bpp = bp + d*bz;
		aold = az;
		am = ap/bpp;
		bm = bp / bpp;
		az = app/bpp;
		bz = 1.0;
		if (fabs(az - aold) < (EPSIL * fabs(az)))
			return az;
		}
	//ErrorBell();
	printf("Internal error in 'BetaCF': MAX_ITER exceeded\n");
	return 0.0;
}
