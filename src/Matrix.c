//
#include "../include/Matrix.h"

void printMatrix(const Matrix A){
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			printf("(%11.4E, %11.4E) ",
			// printf("(%0.3f, %0.3f) ",
			creal(A.data[i][j]), cimag(A.data[i][j]));
		}
		printf("\n");
	}
}

void zerosMatrix(Matrix *A){
	#define A (*A)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			A.data[i][j] = 0.0;
		}
	}
	#undef A
}

void onesMatrix(Matrix *A){
	#define A (*A)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			A.data[i][j] = 1.0;
		}
	}
	#undef A
}

void eyeMatrix(Matrix *A){
	#define A (*A)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			if (i==j){
				A.data[i][j] = 1.0;
			}else{
				A.data[i][j] = 0.0;
			}
		}
	}
	#undef A
}

void allocateMatrix(Matrix *A, int rows, int cols){
	#define A (*A)
	assert(A.data==NULL);
	assert(rows>0&&cols>0);
	//
	A.rows = rows;
	A.cols = cols;
	A.data = malloc(rows*sizeof(complex double*));
	for (int i=0; i<rows; i++){
		A.data[i] = malloc(cols*sizeof(complex double));
	}
	zerosMatrix(&A);
	#undef A
}

void deallocateMatrix(Matrix *A){
	#define A (*A)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	//
	for (int i=0; i<A.rows; i++){
		free(A.data[i]);
	}
	free(A.data);
	A.rows = 0;
	A.cols = 0;
	A.data = NULL;
	#undef A
}

void saveMatrix(const Matrix A, char *fileName){
	FILE *file=fopen(fileName, "w");
	assert(file!=NULL);
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	assert(fileName!=NULL);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			fprintf(file, "(%11.4E, %11.4E) ",
			creal(A.data[i][j]), cimag(A.data[i][j]));
		}
		fprintf(file, "\n");
	}
	fclose(file);
}

void addMatrix(const Matrix A, const Matrix B, Matrix *C){
	// C = A+B
	#define C (*C)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	assert(B.data!=NULL);
	assert(B.rows>0&&B.cols>0);
	assert(C.data!=NULL);
	assert(C.rows>0&&C.cols>0);
	//
	assert(A.rows==B.rows);
	assert(A.rows==C.rows);
	assert(A.cols==B.cols);
	assert(A.cols==C.cols);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			C.data[i][j] = A.data[i][j]+B.data[i][j];
		}
	}
	#undef C
}

void subMatrix(const Matrix A, const Matrix B, Matrix *C){
	// C = A-B
	#define C (*C)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	assert(B.data!=NULL);
	assert(B.rows>0&&B.cols>0);
	assert(C.data!=NULL);
	assert(C.rows>0&&C.cols>0);
	//
	assert(A.rows==B.rows);
	assert(A.rows==C.rows);
	assert(A.cols==B.cols);
	assert(A.cols==C.cols);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			C.data[i][j] = A.data[i][j]-B.data[i][j];
		}
	}
	#undef C
}

void multMatrix(const Matrix A, const Matrix B, Matrix *C){
	// C = A*B
	#define C (*C)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	assert(B.data!=NULL);
	assert(B.rows>0&&B.cols>0);
	assert(C.data!=NULL);
	assert(C.rows>0&&C.cols>0);
	//
	assert(A.cols==B.rows);
	assert(A.rows==C.rows);
	assert(B.cols==C.cols);
	//
	complex double sum;
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<B.cols; j++){
			sum = 0.0;
			for (int k=0; k<A.cols; k++){
				sum+=A.data[i][k]*B.data[k][j];
			}
			C.data[i][j] = sum;
		}
	}
	#undef C
}

void scaleMatrix(Matrix *A, complex double a){
	// A => A*a
	#define A (*A)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			A.data[i][j] = a*A.data[i][j];
		}
	}
	#undef A
}

void copyMatrix(const Matrix A, Matrix *B){
	// A => B
	#define B (*B)
	assert(A.data!=NULL);
	assert(A.rows>0&&A.cols>0);
	assert(B.data!=NULL);
	assert(B.rows>0&&B.cols>0);
	//
	assert(A.rows==B.rows);
	assert(A.cols==B.cols);
	//
	for (int i=0; i<A.rows; i++){
		for (int j=0; j<A.cols; j++){
			B.data[i][j] = A.data[i][j];
		}
	}
	#undef B
}

void PivotSwap(int N, Matrix *A, Matrix *b, int k){
	#define A (*A)
	#define b (*b)
	complex double Amax=A.data[k][k];
	int p=k;
	for (int i=k+1; i<N; i++){
		if (cabs(A.data[i][k])>cabs(Amax)){
			Amax = A.data[i][k];
			p = i;
		}
	}
	if (p!=k){
		complex double D;
		D = b.data[p][0];
		b.data[p][0] = b.data[k][0];
		b.data[k][0] = D;
		for (int j=0; j<N; j++){
			D = A.data[p][j];
			A.data[p][j] = A.data[k][j];
			A.data[k][j] = D;
		}
	}
	#undef A
	#undef b
}

void solveLinearMatrix(const Matrix AIn, const Matrix bIn, Matrix *x){
	// Solve A*x=b
	#define x (*x)
	assert(AIn.data!=NULL);
	assert(AIn.rows>0&&AIn.cols>0);
	assert(bIn.data!=NULL);
	assert(bIn.rows>0&&bIn.cols>0);
	assert(x.data!=NULL);
	assert(x.rows>0&&x.cols>0);
	//
	assert(AIn.rows==AIn.cols);
	int N=AIn.rows;
	assert(bIn.rows==N);
	assert(bIn.cols==1);
	assert(x.rows==N);
	assert(x.cols==1);
	//
	Matrix A=DefaultMatrix;
	Matrix b=DefaultMatrix;
	allocateMatrix(&A, N, N);
	allocateMatrix(&b, N, 1);
	copyMatrix(AIn, &A);
	copyMatrix(bIn, &b);
	// Forward reduction
	for (int k=0; k<N; k++){
		PivotSwap(N, &A, &b, k);
		if (cabs(A.data[k][k])==0.0){
			printf("Matrix is not invertible!\n");
			exit(1);
		}
		for (int i=k+1; i<N; i++){
			b.data[i][0]-=A.data[i][k]*b.data[k][0]/A.data[k][k];
			for (int j=k+1; j<N; j++){
				A.data[i][j]-=A.data[i][k]*A.data[k][j]/A.data[k][k];
			}
			A.data[i][k] = 0.0;
		}
	}
	// Backward substitution
	complex double sum;
	for (int i=N-1; i>-1; i--){
		sum = 0.0;
		for (int j=i+1; j<N; j++){
			sum+=A.data[i][j]*x.data[j][0];
		}
		x.data[i][0] = (b.data[i][0]-sum)/A.data[i][i];
	}
	//
	deallocateMatrix(&A);
	deallocateMatrix(&b);
	#undef x
}

void PivotSwapMatrix(int N, Matrix *A, Matrix *B, int k){
	#define A (*A)
	#define B (*B)
	complex double Amax=A.data[k][k];
	int p=k;
	for (int i=k+1; i<N; i++){
		if (cabs(A.data[i][k])>cabs(Amax)){
			Amax = A.data[i][k];
			p = i;
		}
	}
	if (p!=k){
		complex double D;
		for (int j=0; j<N; j++){
			D = A.data[p][j];
			A.data[p][j] = A.data[k][j];
			A.data[k][j] = D;
			D = B.data[p][j];
			B.data[p][j] = B.data[k][j];
			B.data[k][j] = D;
		}
	}
	#undef A
	#undef B
}

void inverseMatrix(Matrix AIn, Matrix *G){
	#define G (*G)
	assert(AIn.rows==AIn.cols);
	int N=AIn.rows;
	assert((G.rows==G.cols)&&(G.rows==N));
	Matrix A=DefaultMatrix;
	allocateMatrix(&A, N, N);
	//
	copyMatrix(AIn, &A);
	eyeMatrix(&G);
	complex double factor;
	// Forward
	for (int j=0; j<N-1; j++){
		PivotSwapMatrix(N, &A, &G, j);
		for (int i=j+1; i<N; i++){
			factor=A.data[i][j]/A.data[j][j];
			if (cabs(A.data[j][j])==0.0){
				printf("Singular matrix!\n");
				exit(1);
			}
			for (int k=0; k<N; k++){
				A.data[i][k]-=A.data[j][k]*factor;
				G.data[i][k]-=G.data[j][k]*factor;
			}
		}
	}
	// Backward
	for (int j=N-1; j>0; j--){
		PivotSwapMatrix(N, &A, &G, j);
		for (int i=j-1; i>-1; i--){
			factor=A.data[i][j]/A.data[j][j];
			if (cabs(A.data[j][j])==0.0){
				printf("Singular matrix!\n");
				exit(1);
			}
			for (int k=0; k<N; k++){
				A.data[i][k]-=A.data[j][k]*factor;
				G.data[i][k]-=G.data[j][k]*factor;
			}
		}
	}
	// Diagonal
	for (int i=0; i<N; i++){
		factor=A.data[i][i];
		if (cabs(A.data[i][i])==0.0){
			printf("Singular matrix!\n");
			exit(1);
		}
		for (int k=0; k<N; k++){
			A.data[i][k]/=factor;
			G.data[i][k]/=factor;
		}
	}
	deallocateMatrix(&A);
	#undef G
}