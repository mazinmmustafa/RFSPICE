#ifndef MATRIX_H
#define MATRIX_H

// Definitions
#include "myLib.h"

// Matrix
typedef struct Matrix{
	int rows, cols;
	complex double **data;
}Matrix;
#define DefaultMatrix {0, 0, NULL}

// Functions
void printMatrix(const Matrix A);
void zerosMatrix(Matrix *A);
void onesMatrix(Matrix *A);
void eyeMatrix(Matrix *A);
void allocateMatrix(Matrix *A, int rows, int cols);
void deallocateMatrix(Matrix *A);
void saveMatrix(const Matrix A, char *fileName);
void addMatrix(const Matrix A, const Matrix B, Matrix *C);
void subMatrix(const Matrix A, const Matrix B, Matrix *C);
void multMatrix(const Matrix A, const Matrix B, Matrix *C);
void scaleMatrix(Matrix *A, complex double a);
void copyMatrix(const Matrix A, Matrix *B);
void solveLinearMatrix(const Matrix A, const Matrix b, Matrix *x);
void inverseMatrix(Matrix A, Matrix *B);

#endif
