typedef struct {
	int      dim;
	double   *prams;
	double   *lowerBounds;
	double   *upperBounds;
	double   yDiffRequested;
	double   (*funk)(double []);
	double   **p;                     // dim + 1 rows, dim cols
	double   *y;
	int      funkEvaluationCount;
	int      nMaxEvaluationsAllowed;
	double   *pSum;
	double   *pTry;
	int      *compStarts;
	int      compLen;
} simplex;



simplex *newSimplex(
					int       dim,
					double   *prams,
					double   *lowerBounds,
					double   *upperBounds,
					double    yDiffRequested,
					int       nMaxEvaluationsAllowed,
					double   (*funk)(double []),
					int      *compStarts,
					double    startFactor);
void freeSimplex(simplex *aSimp);
double amoeba(simplex *aSimp);
double amotry(simplex *aSimp, int ihi, double fac);
