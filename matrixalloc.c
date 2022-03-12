/*
 * Program: matrixalloc.c
 * Purpose:
 * Goals:
 *  1-)rearrange loops to increase spatial locality of a program
 *  2-)exploit cache to improve the performance of a program
 *  3-)write cache-friendly code for matrix operations
 * Date:   12/03/2021
 */

#include <stdio.h>
#include <stdlib.h>

void initmat(double *, int);
void madd(double *A, double *B, double *C, int N);
void mmul_ijk(double *A, double *B, double *C, int N);
void mmul_ikj(double *A, double *B, double *C, int N);
void mmul_jik(double *A, double *B, double *C, int N);
void mmul_jki(double *A, double *B, double *C, int N);
void mmul_kij(double *A, double *B, double *C, int N);
void mmul_kji(double *A, double *B, double *C, int N);

void main(int argc, char *argv[])
{
   double *A;
   double *B;
   double *C;
   int N = atoi(argv[1]);

   // Allocate N-by-N matrix in the heap.
   A = (double *) malloc(N * N * sizeof(double));
   B = (double *) malloc(N * N * sizeof(double));
   C = (double *) malloc(N * N * sizeof(double));


   // Initialize the matrix.
   initmat(A, N);
   initmat(B, N);

   // Call function
   //mmul_ijk(A,B,C,N);
   //mmul_ikj(A,B,C,N);
   //mmul_jik(A,B,C,N);
   //mmul_jki(A,B,C,N);
   //mmul_kij(A,B,C,N);
   mmul_kji(A,B,C,N);


   // Free the allocated space
   free(A);
   free(B);
   free(C);
}

/*
 * Function:  initmat
 * Purpose:   initilize the matrix
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void initmat(double *mat, int N)
{
   int i, j;

   for (i = 0; i < N; i++)
      for (j = 0; j < N; j++)
   	    mat[i*N+j] = i + j;
}

/*
 * Function:  madd
 * Purpose:   adds two matrixes
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void madd(double *A, double *B, double *C, int N){
  int i, j;
  for (i = 0; i < N; i++){
    for (j = 0; j < N; j++){
      C[i*N+j] = A[i*N+j] + B[i*N+j];
    }
  }
}

/*
 * Function:  mmil_ijk
 * Purpose:   multiplies two matrixes, algorithm(A=row-wise, B=column-wise,C=fixed)
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void mmul_ijk(double *A, double *B, double *C, int N){
  int i, j, k, sum;
  for (i=0; i<N; i++)  {
    for (j=0; j<N; j++) {
      sum = 0.0;
      for (k=0; k<N; k++){
        sum += A[i*N+k] * B[k*N+j];
      }
      C[i*N+j] = sum;
    }
  }
}

/*
 * Function:  mmil_ikj
 * Purpose:   multiplies two matrixes, algorithm(A=fixed, B=row-wise,C=row-wise)
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void mmul_ikj(double *A, double *B, double *C, int N){
  int i, j, k, r;
  for (i=0; i<N; i++)  {
    for (k=0; k<N; k++) {
      r = A[i*N+k];
      for (j=0; j<N; j++){
        C[i*N+j] += r * B[k*N+j];
      }
    }
  }
}

/*
 * Function:  mmil_jik
 * Purpose:   multiplies two matrixes, algorithm(A=row-wise, B=column-wise,C=fixed)
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void mmul_jik(double *A, double *B, double *C, int N){
  int i, j, k, sum;
  for (j = 0; j < N; j++){
    for (i = 0; i < N; i++) {
	     sum = 0.0;
	     for (k = 0; k < N; k++){
	       sum += A[i*N+k]*A[k*N+j];
       }
	     C[i*N+j] += sum;
    }
  }
}

/*
 * Function:  mmil_jki
 * Purpose:   multiplies two matrixes, algorithm(A=column-wise, B=fixe,C=column-wise)
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void mmul_jki(double *A, double *B, double *C, int N){
  int i, j, k, r;
  for (j = 0; j < N; j++){
    for (k = 0; k < N; k++) {
	     r = A[k*N+j];
	     for (i = 0; i < N; i++){
         C[i*N+j] += A[i*N+k]*r;
       }
    }
  }
}

/*
 * Function:  mmil_kij
 * Purpose:   multiplies two matrixes, algorithm(A=fixed, B=row-wise,C=row-wise)
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void mmul_kij(double *A, double *B, double *C, int N){
  int i, j, k, r;
  for (k = 0; k < N; k++){
    for (i = 0; i < N; i++) {
	     r = A[i*N+k];
	     for (j = 0; j < N; j++){
	       C[i*N+j] += r*A[k*N+j];
       }
    }
  }
}

/*
 * Function:  mmil_kji
 * Purpose:   multiplies two matrixes, algorithm(A=column-wise, B=fixed,C=column-wise)
 *
 * Parameters:
 *    mat: input matrix
 *    N: 1 dimensional size of the matrix
 *
 * Returns: NA
 *
 */
void mmul_kji(double *A, double *B, double *C, int N){
  int i, j, k, r;
  for (k = 0; k < N; k++){
    for (j = 0; j < N; j++) {
	     r = A[k*N+j];
	     for (i = 0; i < N; i++){
	       C[i*N+j] += A[i*N+k]*r;
       }
    }
  }
}
