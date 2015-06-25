/* PGF got this code from Dave L Swofford.  Thanks Dave!  He got it
   from a book, and re-coded the fortran into c.  This is not GPL'd; I
   assume it is copyrighted by DLS.  PGF added the newBrent() and
   freeBrent() functions.  And a little more tweaking.

   It is OK with Dave that this code be used as long as acknowlegement
   is given to him.
   -PGF

*/


#include "brent.h"





brent *newBrent(int dim)
{
	brent	*aBrent;

	aBrent = malloc(sizeof(brent));

    aBrent->vv = psdmatrix(dim);

    if ((aBrent->d = (double *)malloc(dim * sizeof(double))) == NULL) {
        printf("praxis: unable to malloc aBrent->d\n");
        exit(1);
    }
    if ((aBrent->q0 = (double *)malloc(dim * sizeof(double))) == NULL) {
        printf("praxis: unable to malloc aBrent->q0\n");
        exit(1);
    }
    if ((aBrent->q1 = (double *)malloc(dim * sizeof(double))) == NULL) {
        printf("praxis: unable to malloc aBrent->q1\n");
        exit(1);
    }
    if ((aBrent->xnew = (double *)malloc(dim * sizeof(double))) == NULL) {
        printf("praxis: unable to malloc aBrent->xnew\n");
        exit(1);
    }
    if ((aBrent->y = (double *)malloc(dim * sizeof(double))) == NULL) {
        printf("praxis: unable to malloc aBrent->y\n");
        exit(1);
    }
    if ((aBrent->z = (double *)malloc(dim * sizeof(double))) == NULL) {
        printf("unable to malloc aBrent->z\n");
        exit(1);
    }



	return aBrent;
}


void freeBrent(brent *aBrent)
{
    if(aBrent->vv) free_psdmatrix(aBrent->vv);
    if(aBrent->d) free(aBrent->d);
    if(aBrent->q0) free(aBrent->q0);
    if(aBrent->q1) free(aBrent->q1);
    if(aBrent->xnew) free(aBrent->xnew);
    if(aBrent->y) free(aBrent->y);
    if(aBrent->z) free(aBrent->z);

	free(aBrent);

}



double praxis(brent *aBrent, double tol, double h, int dim, double *xx, double (*theFunction)(double []))
{
	int			i, j, k, k2, illc, kl, kt, ktm;
	double			scbd, ldfac; 
	double			sf, df, f1, lds, t2, sl, dn, s, sz;
	


	aBrent->eps2 = DBL_EPSILON*DBL_EPSILON;
	aBrent->m2 = sqrt(DBL_EPSILON);
	aBrent->m4 = sqrt(aBrent->m2);
	
	// machine-dependent initializations 
	aBrent->vsmall = aBrent->eps2*aBrent->eps2;
	aBrent->large = 1.0/aBrent->eps2;
	aBrent->vlarge = 1.0/aBrent->vsmall;

	aBrent->n = dim;
	aBrent->fxn = theFunction;
	aBrent->gx = xx;
	aBrent->toler = tol;
	aBrent->htol = h;
	
	
	// heuristic numbers:
	//   - if axes may be badly scaled (which should be avoided if possible), 
	//   		set scbd=10, otherwise 1
	//   - if the problem is known to be ill-conditioned, set illc=TRUE, 
	//   		otherwise FALSE
	//   - ktm+1 is the number of iterations without improvement before 
	//   		the algorithm terminates
	//     (see section 7.6 of Brent's book).  ktm=4 is very cautious; 
	//	 	usually ktm=1 is satisfactory.
	

	scbd = 1.0;
	illc = FALSE;
	ktm = 1;


	ldfac = illc ? 0.1 : 0.01;
	kt = aBrent->nl = 0;
	aBrent->qf1 = aBrent->gfx = (*aBrent->fxn)(aBrent->gx);
	aBrent->toler = t2 = aBrent->eps2 + fabs(aBrent->toler);
	aBrent->dmin = aBrent->eps2;


	if (aBrent->htol < 100.0*aBrent->toler)
		aBrent->htol = 100.0*aBrent->toler;
	aBrent->ldt = aBrent->htol;

	for(i = 0; i < aBrent->n; i++) {
		for(j = 0; j < aBrent->n; j++) {
			aBrent->vv[i][j] = 0.0;
			if(i == j) aBrent->vv[i][j] = 1.0;
		}
	}

	aBrent->d[0] = aBrent->qd0 = 0.0;
	for (i = 0; i < aBrent->n; i++) {
		aBrent->q0[i] = 0.0;	// DLS: this wasn't included in Brent's code, 
						// but it's important
		aBrent->q1[i] = aBrent->gx[i];
	}


	// ------ main loop ------ //

	for(FOREVER)
		{
		//printf("a aBrent->nl = %i, aBrent->gfx = %f, aBrent->gx[0] = %f, aBrent->gx[1] = %f\n", aBrent->nl, aBrent->gfx, aBrent->gx[0], aBrent->gx[1]);
		//printf("a aBrent->nl = %i, aBrent->gx[0] = %f\n", aBrent->nl, aBrent->gx[0]);
		sf = aBrent->d[0];
		aBrent->d[0] = s = 0.0;
		
		// minimize along first direction 
	
		lineMin(aBrent, 0, 2, &aBrent->d[0], &s, aBrent->gfx, FALSE);

		if (s < 0.0)
			{
			for (i = 0; i < aBrent->n; i++)
				aBrent->vv[i][0] = -aBrent->vv[i][0];
			}
		if ((sf <= 0.9*aBrent->d[0]) || (0.9*sf >= aBrent->d[0]))
			{
			for (i = 1; i < aBrent->n; i++)
				aBrent->d[i] = 0.0;
			}
	
		for (k = 1; k < aBrent->n; k++)
			{
			for (i = 0; i < aBrent->n; i++)
				aBrent->y[i] = aBrent->gx[i];
			sf = aBrent->gfx;
			illc = illc || (kt > 0);
	
			for (FOREVER)
				{
				kl = k;
				df = 0.0;
				if (illc)
					{
					// random step to get off resolution valley 
					for (i = 0; i < aBrent->n; i++)
						{
						s = aBrent->z[i] = (0.1*aBrent->ldt + t2*pow(10.0, kt))*(ranDoubleUpToOne() - 0.5);
						for (j = 0; j < aBrent->n; j++)
							aBrent->gx[j] += s*aBrent->vv[j][i];
						// peter was here, Aug 24_00, these next few lines, try to avoid negs
						//for (j = 0; j < aBrent->n; j++) {
						//	if(aBrent->gx[j] < 0.0) {
						//		aBrent->gx[j] = 1.0e-12 * ranDoubleUpToOne();
						//		printf("Brent deneggify 1 here.\n");
						//	}
						//}
						
						}
					aBrent->gfx = (*aBrent->fxn)(aBrent->gx);

					}
				for (k2 = k; k2 < aBrent->n; k2++)
					{
					sl = aBrent->gfx;
					s = 0.0;
					// minimize along "non-conjugate" directions
                    lineMin(aBrent, k2, 2, &aBrent->d[k2], &s, aBrent->gfx, FALSE);
					
					if (illc)
						{
						sz = s + aBrent->z[k2];
						s = aBrent->d[k2]*sz*sz;
						}
					else
						s = sl - aBrent->gfx;
					if (df < s)
						{
						df = s;
						kl = k2;
						}
					}
				if (!illc && (df < fabs(100.0*DBL_EPSILON*aBrent->gfx)))
					illc = TRUE;	// no success with illc=FALSE so try once with illc=TRUE
				else
					break;
				}

			
			for (k2 = 0; k2 < k; k2++)
				{
				// minimize along "conjugate" directions
				s = 0.0;
                lineMin(aBrent, k2, 2, &aBrent->d[k2], &s, aBrent->gfx, FALSE);

				}
	
			f1 = aBrent->gfx;
			aBrent->gfx = sf;
			lds = 0.0;
			for (i = 0; i < aBrent->n; i++)
				{
				sl = aBrent->gx[i];
				aBrent->gx[i] = aBrent->y[i];
				aBrent->y[i] = (sl -= aBrent->y[i]);
				lds += sl*sl;
				}
			lds = sqrt(lds);
			if (lds > aBrent->eps2)
				{
				// throw away direction kl 
				for (i = kl - 1; i >= k; i--)
					{
					for (j = 0; j < aBrent->n; j++)
						aBrent->vv[j][i + 1] = aBrent->vv[j][i];
					aBrent->d[i + 1] = aBrent->d[i];
					}
					
				// set new "conjugate" aBrent->direction ... 
				aBrent->d[k] = 0.0;
				for (i = 0; i < aBrent->n; i++)
					aBrent->vv[i][k] = aBrent->y[i]/lds;
				
				// ... and minimize along it
                lineMin(aBrent, k, 4, &aBrent->d[k], &lds, f1, TRUE);

				if (lds <= 0.0)
					{
					lds = -lds;
					for (i = 0; i < aBrent->n; i++)
						aBrent->vv[i][k] = -aBrent->vv[i][k];
					}
				}
			aBrent->ldt *= ldfac;
			if (aBrent->ldt < lds)
				aBrent->ldt = lds;

			t2 = 0.0;
			for (i = 0; i < aBrent->n; i++)
				t2 += aBrent->gx[i]*aBrent->gx[i];
			t2 = aBrent->m2*sqrt(t2) + aBrent->toler;
			
			// see if step length exceeds half the tolerance 
			kt = (aBrent->ldt > 0.5*t2) ? 0 : kt + 1;
			if (kt > ktm)
				{
				return aBrent->gfx;
				}
			}
		
		// try quadratic extrapolation in case we are stuck in a curved valley 		
		qUad(aBrent);
	
		// calculate V = U.(D^(-1/2))  (note: 'v' currently contains U) 
	
		dn = 0.0;
		for (i = 0; i < aBrent->n; i++)
			{
			aBrent->d[i] = 1.0/sqrt(aBrent->d[i]);
			if (aBrent->d[i] > dn)
				dn = aBrent->d[i];
			}
		for (j = 0; j < aBrent->n; j++)
			{
			s = aBrent->d[j]/dn;
			for (i = 0; i < aBrent->n; i++)
				aBrent->vv[i][j] *= s;
			}
	
		if (scbd > 1.0)
			{
			// scale axes in attempt to reduce condition number 
			s = aBrent->vlarge;
			for (i = 0; i < aBrent->n; i++)
				{
				sl = 0.0;
				for (j = 0; j < aBrent->n; j++)
					sl += aBrent->vv[i][j]*aBrent->vv[i][j];
				aBrent->z[i] = sqrt(sl);
				if (aBrent->z[i] < aBrent->m4)
					aBrent->z[i] = aBrent->m4;
				if (s > aBrent->z[i])
					s = aBrent->z[i];
				}
			for (i = 0; i < aBrent->n; i++)
				{
				sl = s/aBrent->z[i];
				aBrent->z[i] = 1.0/sl;
				if (aBrent->z[i] > scbd)
					{
#					if 0	// DLS/POL 1/6/97: Borland C says this is a do-nothing statement,
                           			//                   and it looks like they're right
					sl = 1.0/scbd;
#					endif
					aBrent->z[i] = scbd;
					}
				}
			}
			
		// Find the singular value decomposition of v.  This gives the eigenvalues and principal axes
		//   of the approximating quadratic form without squaring the condition number. 
		// transpose v for MinFit 
		for (i = 1; i < aBrent->n; i++)
			{
			for (j = 0; j < i; j++)
				{
				s = aBrent->vv[i][j];
				aBrent->vv[i][j] = aBrent->vv[j][i];
				aBrent->vv[j][i] = s;
				}
			}
		minFit(aBrent, DBL_EPSILON, aBrent->vsmall, aBrent->vv, aBrent->d, aBrent->y);	// ('y' is just a scratch vector) 
		
		if (scbd > 1.0)
			{
			// unscaling 
			for (i = 0; i < aBrent->n; i++)
				{
				s = aBrent->z[i];
				for (j = 0; j < aBrent->n; j++)
					aBrent->vv[i][j] *= s;
				}
			for (i = 0; i < aBrent->n; i++)
				{
				s = 0.0;
				for (j = 0; j < aBrent->n; j++)
					s += aBrent->vv[j][i]*aBrent->vv[j][i];
				s = sqrt(s);
				aBrent->d[i] *= s;
				s = 1.0/s;
				for (j = 0; j < aBrent->n; j++)
					aBrent->vv[j][i] *= s;
				}
			}
	
		for (i = 0; i < aBrent->n; i++)
			{
			s = dn*aBrent->d[i];
			if (s > aBrent->large)
				aBrent->d[i] = aBrent->vsmall;
			else if (s < aBrent->eps2)
				aBrent->d[i] = aBrent->vlarge;
			else
				aBrent->d[i] = 1.0/(s*s);
			}
			
		// sort new eigenvalues and eigenvectors
        sortDV(aBrent);
		
		aBrent->dmin = aBrent->d[aBrent->n - 1];
		if (aBrent->dmin < aBrent->eps2)
			aBrent->dmin = aBrent->eps2;
		
		illc = (aBrent->m2*aBrent->d[0] > aBrent->dmin);
		}

}

void minFit(brent *aBrent, double eps, double tol, double **ab, double *q, double *e)
{
	int			i, j, k, l, l2, kt;
	double		c, f, g, h, s, x, y2, z2;

	l = 0; // to shut the compiler up.
	/* Householder's reduction to bidiagonal form */
	g = x = 0.0;
	for (i = 0; i < aBrent->n; i++)
		{
		e[i] = g;
		s = 0.0;
		l = i + 1;
		for (j = i; j < aBrent->n; j++)
			s += ab[j][i]*ab[j][i];
		if (s < tol)
			g = 0.0;
		else
			{
			f = ab[i][i];
			g = (f < 0.0) ? sqrt(s) : -sqrt(s);
			h = f*g - s;
			ab[i][i] = f - g;
			for (j = l; j < aBrent->n; j++)
				{
				f = 0.0;
				for (k = i; k < aBrent->n; k++)
					f += ab[k][i]*ab[k][j];
				f /= h;
				for (k = i; k < aBrent->n; k++)
					ab[k][j] += f*ab[k][i];
				}
			}
		q[i] = g;
		s = 0.0;
		for (j = l; j < aBrent->n; j++)
			s += ab[i][j]*ab[i][j];
		
		if (s < tol)
			g = 0.0;
		else
			{
			f = ab[i][i+1];
			g = (f < 0.0) ? sqrt(s) : -sqrt(s);
			h = f*g - s;
			ab[i][i+1] = f - g;
			for (j = l; j < aBrent->n; j++)
				e[j] = ab[i][j]/h;
			for (j = l; j < aBrent->n; j++)
				{
				s = 0.0;
				for (k = l; k < aBrent->n; k++)
					s += ab[j][k]*ab[i][k];
				for (k = l; k < aBrent->n; k++)
					ab[j][k] += s*e[k];
				}
			}
		y2 = fabs(q[i]) + fabs(e[i]);
		if (y2 > x)
			x = y2;
		}

	/* accumulation of right-hand transformations */
	for (i = aBrent->n-1; i >= 0; i--)
		{
		if (g != 0.0)
			{
			h = ab[i][i+1]*g;
			for (j = l; j < aBrent->n; j++)
				ab[j][i] = ab[i][j]/h;
			for (j = l; j < aBrent->n; j++)
				{
				s = 0.0;
				for (k = l; k < aBrent->n; k++)
					s += ab[i][k]*ab[k][j];
				for (k = l; k < aBrent->n; k++)
					ab[k][j] += s*ab[k][i];
				}
			}
		for (j = l; j < aBrent->n; j++)
			ab[i][j] = ab[j][i] = 0.0;
		ab[i][i] = 1.0;
		g = e[i];
		l = i;
		}
		
	/* diagonalization of the bidiagonal form */
	eps *= x;
	for (k = aBrent->n-1; k >= 0; k--)
		{
		kt = 0;
		
		test_splitting:

		if (++kt > 30)
			{
			e[k] = 0.0;
			}
		for (l2 = k; l2 >= 0; l2--)
			{
			l = l2;
			if (fabs(e[l]) <= eps)
				goto test_convergence;
			if (fabs(q[l-1]) <= eps)
				break;
			}

		/* cancellation of e[l] if l > 1 */
		c = 0.0;
		s = 1.0;
		for (i = l; i <= k; i++)
			{
			f = s*e[i];
			e[i] *= c;
			if (fabs(f) <= eps)
				break;
			g = q[i];
			if (fabs(f) < fabs(g))
				h = fabs(g)*sqrt(1.0 + (f/g)*(f/g));
			else if (f != 0.0)
				h = fabs(f)*sqrt(1.0 + (g/f)*(g/f));
			else
				h = 0.0;
			q[i] = h;
			if (h == 0.0)
				g = h = 1.0;	/* note: this replaces q[i]=h=sqrt(g*g+f*f) which may give
				                   incorrect results if the squares underflow or if f=g=0 */
			c = g/h;
			s = -f/h;
			}
			
		test_convergence:
		
		z2 = q[k];
		if (l != k)
			{			
			/* shift from bottom 2*2 minor */
			
			x = q[l];
			y2 = q[k-1];
			g = e[k-1];
			h = e[k];
			f = ((y2-z2)*(y2+z2) + (g-h)*(g+h))/(2.0*h*y2);
			g = sqrt(f*f + 1.0);
			s = (f < 0.0) ? f - g : f + g;
			f = ((x-z2)*(x+z2) + h*(y2/s - h))/x;
				
			/* next QR transformation */
			
			c = s = 1.0;
			for (i = l + 1; i <= k; i++)
				{
				g = e[i];
				y2 = q[i];
				h = s*g;
				g *= c;
				if (fabs(f) < fabs(h))
					z2 = fabs(h)*sqrt(1.0 + (f/h)*(f/h));
				else if (f != 0.0)
					z2 = fabs(f)*sqrt(1.0 + (h/f)*(h/f));
				else
					z2 = 0.0;
				e[i-1] = z2;
				if (z2 == 0.0)
					z2 = f = 1.0;
				c = f/z2;
				s = h/z2;
				f = x*c + g*s;
				g = -x*s + g*c;
				h = y2*s;
				y2 *= c;
				for (j = 0; j < aBrent->n; j++)
					{
					x = ab[j][i-1];
					z2 = ab[j][i];
					ab[j][i-1] = x*c + z2*s;
					ab[j][i] = -x*s + z2*c;
					}
				if (fabs(f) < fabs(h))
					z2 = fabs(h)*sqrt(1.0 + (f/h)*(f/h));
				else if (f != 0.0)
					z2 = fabs(f)*sqrt(1.0 + (h/f)*(h/f));
				else
					z2 = 0.0;
				q[i-1] = z2;
				if (z2 == 0.0)
					z2 = f = 1.0;
				c = f/z2;
				s = h/z2;
				f = c*g + s*y2;
				x = -s*g + c*y2;
				}
			e[l] = 0.0;
			e[k] = f;
			q[k] = x;
			goto test_splitting;
			}
		
		if (z2 < 0.0)
			{
			/* q[k] is made non-negative */
			q[k] = -z2;
			for (j = 0; j < aBrent->n; j++)
				ab[j][k] = -ab[j][k];
			}
		}
}

void sortDV(brent *aBrent)
{
	int		i, j, k;
	double	s;

	for (i = 0; i < aBrent->n - 1; i++)
		{
		k = i;
		s = aBrent->d[i];
		for (j = i + 1; j < aBrent->n; j++)
			{
			if (aBrent->d[j] > s)
				{
				k = j;
				s = aBrent->d[j];
				}
			}
		if (k > i)
			{
			aBrent->d[k] = aBrent->d[i];
			aBrent->d[i] = s;
			for (j = 0; j < aBrent->n; j++)
				{
				s = aBrent->vv[j][i];
				aBrent->vv[j][i] = aBrent->vv[j][k];
				aBrent->vv[j][k] = s;
				}
			}
		}
}

void lineMin(brent *aBrent, int j, int nits, double *pd2, double *px1, double f1, int fk)
{
	int		i, k, need_d2z, success;
	double	x1, x2, xm, f0, f2, fm, d1, d2, t2, s, sf1, sx1;

	/* copy args passed by reference to locals (will pass back at end) */
	d2 = *pd2;
	x1 = *px1;
	
	sf1 = f1;
	sx1 = x1;
	k = 0;
	xm = 0.0;
	f0 = fm = aBrent->gfx;
	need_d2z = (d2 < DBL_EPSILON);	/* if TRUE, we need f''(0) */
	
	/* find step size */
	s = 0.0;
	for (i = 0; i < aBrent->n; i++)
		s += aBrent->gx[i]*aBrent->gx[i];
	s = sqrt(s);
	t2 = aBrent->m4*sqrt(fabs(aBrent->gfx)/(need_d2z ? aBrent->dmin : d2) + s*aBrent->ldt) + aBrent->m2*aBrent->ldt;
	s = aBrent->m4*s + aBrent->toler;
	if (need_d2z && (t2 > s))
		t2 = s;
	if (t2 < aBrent->eps2)
		t2 = aBrent->eps2;
	if (t2 > 0.01*aBrent->htol)
		t2 = 0.01*aBrent->htol;
	if (fk && (f1 <= fm))
		{
		xm = x1;
		fm = f1;
		}
	if (!fk || (fabs(x1) < t2))
		{
		x1 = (x1 >= 0.0) ? t2 : -t2;
        f1 = fLin(aBrent, j, x1);
		}
	if (f1 <= fm)
		{
		xm = x1;
		fm = f1;
		}
	
	/* find a distance x2 ("lambda*") that approximately minimizes f in the chosen direction */
	
	do	{
		if (need_d2z)
			{
			/* evaluate FLin at another point and estimate the second derivative */
			x2 = (f0 < f1) ? -x1 : 2.0*x1;
            f2 = fLin(aBrent, j, x2);
			//f2 = [self flinWithJ: j lambda: x2];
			if (f2 <= fm)
				{
				xm = x2;
				fm = f2;
				}
			d2 = (x2*(f1 - f0) - x1*(f2 - f0))/(x1*x2*(x1 - x2));
			}
	
		/* estimate first derivative at 0 */	
		d1 = (f1 - f0)/x1 - x1*d2;
		need_d2z = TRUE;			/* reset flag in case we don't exit loop */
		
		/* predict minimum */
		if (d2 <= aBrent->eps2)
			x2  = (d1 < 0.0) ? aBrent->htol : -aBrent->htol;
		else
			x2 = -0.5*d1/d2;
		if (fabs(x2) > aBrent->htol)
			x2 = (x2 > 0.0) ? aBrent->htol : -aBrent->htol;
		
		/* evaluate f at predicted minimum */
		do	{
            f2 = fLin(aBrent, j, x2);
			//f2 = [self flinWithJ: j lambda: x2];

			success = TRUE;
			if ((k < nits) && (f2 > f0))
				{
				/* no success so halve interval and try again */
				success = FALSE;
				k++;
				if ((f0 < f1) && (x1*x2 > 0.0))
					break;
				x2 *= 0.5;
				}
			}
			while (!success);
		}		
		while (!success);

	aBrent->nl++;	/* increment one-dimensional search counter */
	if (f2 > fm)
		x2 = xm;
	else
		fm = f2;
		
	/* get new estimate of second derivative */
	if (fabs(x2*(x2 - x1)) > aBrent->eps2)
		d2 = (x2*(f1 - f0) - x1*(fm - f0))/(x1*x2*(x1 - x2));
	else if (k > 0)
		d2 = 0.0;
	if (d2 < aBrent->eps2)
		d2 = aBrent->eps2;
	x1 = x2;
	aBrent->gfx = fm;
	if (sf1 < aBrent->gfx)
		{
		aBrent->gfx = sf1;
		x1 = sx1;
		}
		
	/* update x for linear search but not for parabolic search */
	if (j >= 0)
		{
		for (i = 0; i < aBrent->n; i++)
			aBrent->gx[i] += x1*aBrent->vv[i][j];
		}

	*px1 = x1;
	*pd2 = d2;
}

double fLin(brent *aBrent, int j, double lambda)
{
	int		i;
	double	qa, qb, qc;

	if (j >= 0)
		{
		/* linear search */
		for (i = 0; i < aBrent->n; i++)
			aBrent->xnew[i] = aBrent->gx[i] + lambda*aBrent->vv[i][j];
		}
	else
		{
		/* search along a parabolic space curve */
		qa = lambda*(lambda - aBrent->qd1)/(aBrent->qd0*(aBrent->qd0 + aBrent->qd1));
		qb = (lambda + aBrent->qd0)*(aBrent->qd1 - lambda)/(aBrent->qd0*aBrent->qd1);
		qc = lambda*(lambda + aBrent->qd0)/(aBrent->qd1*(aBrent->qd0 + aBrent->qd1));
		
		/* previous three points were stored as follows: x' in q0, x'' in gx, and x''' in q1;
		   see comments in 'Quad' */
		for (i = 0; i < aBrent->n; i++)
			aBrent->xnew[i] = qa*aBrent->q0[i] + qb*aBrent->gx[i] + qc*aBrent->q1[i];
		}

		return (*aBrent->fxn)(aBrent->xnew);
}

void qUad(brent *aBrent)
{
	int		i;
	double	lambda, s;
	double	qa, qb, qc;
	
	/* q0 and q1 contain previous two points */

	s = aBrent->gfx;
	aBrent->gfx = aBrent->qf1;
	aBrent->qf1 = s;
	aBrent->qd1 = 0.0;
	for (i = 0; i < aBrent->n; i++)
		{
		/* copy x to q1 for use in next cycle (but save current q1 in x so we can calculate the
		   norm below) */
		s = aBrent->gx[i];
		aBrent->gx[i] = aBrent->q1[i];
		aBrent->q1[i] = s;			/* copy original x to q1 */
		aBrent->qd1 += (aBrent->q1[i] - aBrent->gx[i])*(aBrent->q1[i] - aBrent->gx[i]);
		}
	aBrent->qd1 = sqrt(aBrent->qd1);
	if ((aBrent->qd0 > 0.0) && (aBrent->qd1 > 0.0) && (aBrent->nl >= 3*aBrent->n*aBrent->n))
		{
		s = 0.0;
		lambda = aBrent->qd1;
        lineMin(aBrent, -1, 2, &s, &lambda, aBrent->qf1, TRUE);
		qa = lambda*(lambda - aBrent->qd1)/(aBrent->qd0*(aBrent->qd0 + aBrent->qd1));
		qb = (lambda + aBrent->qd0)*(aBrent->qd1 - lambda)/(aBrent->qd0*aBrent->qd1);
		qc = lambda*(lambda + aBrent->qd0)/(aBrent->qd1*(aBrent->qd0 + aBrent->qd1));
		}
	else
		{
		aBrent->gfx = aBrent->qf1;
		qa = qb = 0.0;
		qc = 1.0;
		}
	aBrent->qd0 = aBrent->qd1;
	for (i = 0; i < aBrent->n; i++)
		{
		s = aBrent->q0[i];							/* save current q0 for calculation below */
		aBrent->q0[i] = aBrent->gx[i];						/* copy current q1 (now in gx) to q0 for next cycle */
		aBrent->gx[i] = qa*s + qb*aBrent->gx[i] + qc*aBrent->q1[i];	/* gx now contains q1, and q1 contains gx */
		}
}




