#ifndef MATRIX
#define MATRIX

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* 
 * Allocates int matrix of size m x n on heap. Organized in row-major order.
 * Returns pointer to array of int* (pointers to each row).
 */
int** get_int_matrix(const int m, const int n) {
	int** x = malloc(m * sizeof(int*));
	if (x == NULL) {
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	for (int i=0; i< m; i++) {
		x[i] = calloc(n, sizeof(int));
		if (x[i] == NULL) {
			fprintf(stderr, "Error: Out of memory!\n");
			exit(1);
		}
	}

	return x;
}

/* 
 * Allocates int array of length n on heap. Returns pointer to array.
 */
int* get_int_array(const int n) {
	int* x = calloc(n, sizeof(int));
	if (x == NULL){
		fprintf(stderr, "Error: Out of memory.\n");
		exit(1);
	}
	return x;
}


/* 
 * Allocates int matrix of size m x n on heap. Organized in row-major order.
 * Returns pointer to array of int* (pointers to each row).
 */
double** get_double_matrix(const int m, const int n) {
	double** x = malloc(m * sizeof(double*));
	if (x == NULL) {
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}

	for (int i=0; i< m; i++) {
		x[i] = calloc(n, sizeof(double));
		if (x[i] == NULL) {
			fprintf(stderr, "Error: Out of memory!\n");
			exit(1);
		}
	}

	return x;
}

/* 
 * Allocates double array of length n on heap. Returns pointer to array.
 */
double* get_double_array(const int n) {
	double* x = calloc(n, sizeof(double));
	if (x == NULL) {
		fprintf(stderr, "Error: Out of memory!\n");
		exit(1);
	}
	return x;
}


// Get a double cube, of dimensions m x n x o. 
double*** get_double_cube(const int m, const int n, const int o) {
	double*** x = malloc(m * sizeof(double**));
	if (x == NULL) {
		fprintf(stderr, "Error: Out of memory.\n");
		exit(1);
	}
	for (int i = 0; i < m; i++) {
		x[i] = get_double_matrix(n, o);
	}
	return x;
}


/* 
 * Prints int matrix, primarily for debugging purposes.
 */
void print_int_matrix(int** mat, int m, int n) {
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			printf(" %d", mat[i][j]);
		}
		printf("\n");
	}
}

void print_double_matrix(double** mat, int m, int n) {
	for (int i = 0; i< m; i++) {
		for (int j = 0; j< n; j++) {
			printf(" %f", mat[i][j]);
		}
		printf("\n");
	}
}

/* 
 * Prints int array, primarily for debugging purposes.
 */
void print_int_array(int* arr, int n) {
	for (int i = 0; i < n; i++)
		printf(" %d", arr[i]);
	printf("\n");
}

void print_double_array(double* arr, int n) {
	for (int i = 0; i < n; i++)
		printf(" %f", arr[i]);
	printf("\n");
}


// m = number of rows, since stored in row major order
void free_int_matrix(int** matrix, const int m) {
	for (int i = 0; i < m; i++) {
		free(matrix[i]);
	}
	free(matrix);
}

void free_double_matrix(double** matrix, const int m) {
	for (int i = 0; i < m; i++) {
		free(matrix[i]);
	}
	free(matrix);
}

void free_double_cube(double*** cube, const int m, const int n) {
	for(int i = 0; i < m; i++) {
		free_double_matrix(cube[i], n);
	}
	free(cube);
}


void free_int_array(int* arr) {
	free(arr);
}

void free_double_array(double* arr) {
	free(arr);
}


bool equal_int_matrix(int** S, int** T, int n) {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (S[i][j] != T[i][j])
				return false;
		}
	}
	return true;
}

// Boolean matrix sum add1 + add2 = sum. 0 if both 0, otherwise 1. 
void boolean_matrix_sum(int** add1, int** add2, int** sum, int n) {
	for (int i = 0; i< n; i++) {
		for (int j = 0; j < n; j++) {
			if (add1[i][j] || add2[i][j]) {
				sum[i][j] = 1;
			} else {
				sum[i][j] = 0;
			}
		}
	}
}

// Boolean matrix product -- factor1 * factor2 = prod. prod cell = 1 if product is nonzero, otherwise 0.
void boolean_matrix_product(int** factor1, int** factor2, int** prod, int n) {
	for (int i = 0; i< n; i++) {
		for (int j = 0; j< n; j++) {
			bool changed = false;
			for (int k = 0; k< n; k++) {
				if (factor1[i][k] && factor2[k][j]) { // if both are nonzero, then product will be nonzero
					prod[i][j] = 1;
					changed = true;
					break;
				}
			}
			if (!changed) prod[i][j] = 0;
		}
	}
}

// for optimize, why not use libraries?
// T is the transitive closure of relation A. Does not modify matrix A. 
void transitive_closure(int** A, int** T, int n) {
	int** R = get_int_matrix(n, n);
	int** S = get_int_matrix(n, n);
	for (int i = 0; i< n; i++) {
		for (int j = 0; j< n; j++) {
			T[i][j] = A[i][j];
		}
	} 
	
	bool cont = true;
	while(cont) { // run at least once
		
		// Set S= T
		for (int i = 0; i< n; i++) {
			memmove(S[i], T[i], sizeof(int) * n);
		}
		// T = A*S + A
		boolean_matrix_product(A, S, R, n);
		boolean_matrix_sum(R, A, T, n);
		cont = !equal_int_matrix(S, T, n);
	}

	free_int_matrix(R, n);
	free_int_matrix(S, n);
}


int abs_matrix_diff(int** A, int** B, int n) {
	int diff = 0; 
	for(int i = 0; i< n; i++) {
		for(int j = 0; j< n; j++) {
			if (A[i][j] != B[i][j]) diff++;
		}
	}
	return diff;
}

void fprint_int_matrix(FILE* log, int** mat, const int m, const int n) {
	for (int i = 0; i < m; i++) {
		fprintf(log, "%d", mat[i][0]);
		for (int j = 1; j < n; j++) {
			fprintf(log, ",%d", mat[i][j]);
		}
		fprintf(log, "\n");
	}
}

void fprint_int_array(FILE* log, int* arr, const int n) {
	fprintf(log, "%d", arr[0]);
	for (int i = 1; i < n; i++) {
		fprintf(log, ",%d", arr[i]);
	}
	fprintf(log, "\n");
}

int** copy_int_matrix(int** mat, const int m, const int n) {
	int** copy = get_int_matrix(m, n);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j< n; j++) {
			copy[i][j] = mat[i][j];
		}
	}
	return copy;
}

void fprint_double_array(FILE* log, double* arr, const int n) {
	for (int i = 0; i < n; i++) {
		fprintf(log, "%g ", arr[i]);
	}
	fprintf(log, "\n");
}

#endif
