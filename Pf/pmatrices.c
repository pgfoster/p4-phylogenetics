#include <stdio.h>	// for printf
#include "pmatrices.h"
#include <stdlib.h>	// for malloc and free


double **psdmatrix(int dim) // square double
{
	int i;
	double **m;
	/* allocate pointers to rows */
	m=(double **) malloc(dim * sizeof(double*));
	if (!m) {
		printf("allocation error in psdmatrix 1.\n");
		exit(1);
	}
	// allocate rows and set pointers to them 
	m[0]=(double *) malloc(dim * dim * sizeof(double));
	if (!m[0]) {
		printf("allocation error in psdmatrix 2.\n");
		exit(1);
	}
	for(i=1;i<dim;i++) {
		m[i]=m[i-1]+dim;
	}
	// return pointer to array of pointers to rows
	return m;
}

void free_psdmatrix(double **m)
// free a double matrix allocated by psdmatrix()
// this works - I don't know why its a char *, maybe it should be a void *?

{
	free((char *) (m[0]));
	free((char *) (m));
}

void dump_psdmatrix(double **m, int dim)
{
	int	row, col;

#if 0
	printf("{");
	for(row = 0; row < (dim - 1); row++) {
		printf("{");
		for(col = 0; col < (dim - 1); col++) {
			printf("%f,", m[row][col]);
		}
		printf("%f},\n", m[row][dim - 1]);
	}
	printf("{");
	for(col = 0; col < (dim - 1); col++) {
		printf("%f,", m[dim - 1][col]);
	}
	printf("%f}}", m[dim - 1][dim - 1]);
	printf("\n");
#endif

	for(row = 0; row < dim; row++) {
		for(col = 0; col < dim; col++) {
			printf(" & %8.6f", m[row][col]);
		}
		printf("\n");
	}
}

void copy_psdmatrix(double **from, double **to, int dim)
{
	int row, col;

	for(row = 0; row < dim; row++) {
		for(col = 0; col < dim; col++) {
			to[row][col] = from[row][col];
		}
	}

}

double **pdmatrix(int rows, int cols)
// allocate a double rectangular matrix with subscript 
// range m[0..rows-1][0..cols-1]
{
	int i;
	double **m;
	// allocate pointers to rows 
	m=(double **) malloc(rows * sizeof(double*));
	if (!m) {
		printf("allocation error in pdmatrix 1.\n");
		exit(1);
	}
	// allocate rows and set pointers to them 
	m[0]=(double *) malloc(rows * cols * sizeof(double));
	if (!m[0]) {
		printf("allocation error in pdmatrix 2.\n");
		exit(1);
	}
	for(i=1;i<rows;i++) {
		m[i]=m[i-1]+cols;
	}
	// return pointer to array of pointers to rows
	return m;
}

void free_pdmatrix(double **m)
{
	free((char *) (m[0]));
	free((char *) (m));
}



int **psimatrix(int dim)
// allocate an int square matrix with subscript 
// range m[0..dim-1][0..dim-1]
{
	int i;
	int **m;
	// allocate pointers to rows 
	m=(int **) malloc(dim * sizeof(int*));
	if (!m) {
		printf("allocation error in psimatrix 1.\n");
		exit(1);
	}
	// allocate rows and set pointers to them 
	m[0]=(int *) malloc(dim * dim * sizeof(int));
	if (!m[0]) {
		printf("allocation error in psimatrix 2.\n");
		exit(1);
	}
	for(i=1;i<dim;i++) {
		m[i]=m[i-1]+dim;
	}
	// return pointer to array of pointers to rows
	return m;
}

void free_psimatrix(int **m)
{
	free((char *) (m[0]));
	free((char *) (m));
}

void dump_psimatrix(int **m, int dim)
{
	int	row, col;
	
	printf("{");
	for(row = 0; row < (dim - 1); row++) {
		printf("{");
		for(col = 0; col < (dim - 1); col++) {
			printf("%i,", m[row][col]);
		}
		printf("%i},\n", m[row][dim - 1]);
	}
	printf("{");
	for(col = 0; col < (dim - 1); col++) {
		printf("%i,", m[dim - 1][col]);
	}
	printf("%i}}", m[dim - 1][dim - 1]);
	printf("\n");
}

int **pimatrix(int rows, int cols)
// allocate an int rectangular matrix with subscript
// range m[0..rows-1][0..cols-1]
{
    int i;
    int **m;
    /* allocate pointers to rows */
    m=(int **) malloc(rows * sizeof(int*));
    if (!m) {
        printf("allocation error in pimatrix 1.\n");
        exit(1);
    }
    // allocate rows and set pointers to them
    m[0]=(int *) malloc(rows * cols * sizeof(int));
    if (!m[0]) {
        printf("allocation error in pimatrix 2.\n");
        exit(1);
    }
    for(i=1;i<rows;i++) {
        m[i]=m[i-1]+cols;
    }
    // return pointer to array of pointers to rows
    return m;
}

void free_pimatrix(int **m)
{
    free((char *) (m[0]));
    free((char *) (m));
}




double *pdvector(int dim)
{
	double *vec;
	
	vec = (double *)malloc(dim * sizeof(double));
	if(!vec) printf("memory allocation error, pdvector\n");
	return vec;
}

void dump_pdvector(double *vec, int dim)
{
	int	k;
	
	//printf("{");
	for(k = 0; k < dim; k++) {
		printf(" & %8.6f", vec[k]);
	}
	//printf("%f}\n", vec[dim - 1]);
	printf("\n");
}
	


int *pivector(int dim)
{
	int *vec;
	
	vec = (int *)malloc(dim * sizeof(int));
	if(!vec) printf("memory allocation error, pivector\n");
	return vec;
}

void transpose_psdmatrix(double **m, int dim)
{
	double temp;
	int row, col;
	
	for(row = 0; row < dim; row++) {
		for(col = 0; col < row; col++) {
			temp = m[row][col];
			m[row][col] = m[col][row];
			m[col][row] = temp;
		}
	}
	
}

double *dotMultMatrixByVector(double **mat, double *vec, int dim)
{
	double *result;
	int	row, col;
	
	result = pdvector(dim);
	for(row = 0; row < dim; row++) {
		result[row] = 0.0;
		for(col = 0; col < dim; col++) {
			result[row] = result[row] + (mat[row][col] * vec[col]);
		}
	}
	return result;

}

void dotMultMatrixByMatrix(double **leftmat, double **rightmat, int dim, double **result)
{
	int	row, col, k;
	
	for(row = 0; row < dim; row++) {
		for(col = 0; col < dim; col++) {
			result[row][col] = 0.0;
			for(k = 0; k < dim; k++) {
				result[row][col] = result[row][col] + 
						(leftmat[row][k] * rightmat[k][col]);
			}
		}
	}
}

void	multMatrixColumnsByVector(double **leftmat, double *rightvec, int dim, double **result)
{
	int	row, col;
	
	for(row = 0; row < dim; row++) {
		for(col = 0; col < dim; col++) {
			result[row][col] = leftmat[row][col] * rightvec[col];
		}
	}
}

void diagonalFillMatrix(double **sqMatrix, int dim)
{
	int i, j;
	
	for(i = 0; i < dim; i ++){
		for(j = 0; j < dim; j++){
			if(i == j) {
				sqMatrix[i][j] = 1.0;
			} else {
				sqMatrix[i][j] = 0.0;
			}
		}
	}
}

void zeroFillSquareMatrix(double **sqMatrix, int dim)
{
    int i, j;

    for(i = 0; i < dim; i ++){
        for(j = 0; j < dim; j++){
			sqMatrix[i][j] = 0.0;
        }
    }
}

