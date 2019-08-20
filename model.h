/**
 * Model header file for model struct and helper functions. 
 * Based on implementation by Thomas Sakoparnig, based on the ct_cbn package.
 * 
 * Kevin Chen, July 6th, 2016.
 */

#ifndef MODEL_H_
#define MODEL_H_
#include "h-esbcn.h"

typedef struct {
	int** P;  // cover relations of the event poset

	int n;  // number of events

	/* EVENTS (NODES) -- Length n */
	int* num_pa;  // number of parents for each node (event). has length num_events.
	int** pa;  // pa[i] is the parent set of event i, where 0 <= i < num_events

	/* GENOTYPES -- Length <= 2^(n) - 1*/
	int* num_ch; // number of children for each genotype
	int** ch_diff; // contains lists of differing indexes of children for ith genotype in valids.
	int* g_num_pa; // number of parents for each genotype
	int** g_pa; // contains lists of integer representations of parents. 
	int** g_pa_diff;

	int m; // number of valid genotypes
	int* valids; // list of all valid genotypes

} model;


void free_model(model M, bool em_fields) {
	free_int_matrix(M.P, M.n); 
	free_int_matrix(M.pa, M.n);
	free(M.num_pa);

	// need to check if they're initialized or not. These are throwing issues in the MCMC since they're not allocated. 
	if (em_fields) {
		free(M.num_ch);
		free_int_matrix(M.ch_diff, M.m);
		free(M.g_num_pa);
		free_int_matrix(M.g_pa, M.m);
		free_int_matrix(M.g_pa_diff, M.m);
		free(M.valids);
	}

}

// Functions like strcpy - takes in src and destination model pointers, copies information from one model to the other
// Allocates new arrays for heap elements, be sure to free duplicate as well. 
void model_copy(model* src, model* dest) {
	int n = src->n;
	dest->n = n;

	//duplicate cover relations and linear extension
	int** P = get_int_matrix(n, n);
	for (int i= 0; i< n; i++) {
		for (int j = 0; j< n; j++) {
			P[i][j] = src->P[i][j];
		}
	}
	dest->P = P;

	//duplicate parents array
	int* num_pa = get_int_array(n);
	int** pa = malloc(n * sizeof(int*));
	for (int i = 0; i< n; i++) {
		int num_parents = src->num_pa[i];
		num_pa[i] = num_parents;
		pa[i] = get_int_array(num_parents);
		for (int j = 0; j< num_parents; j++) {
			pa[i][j] = src->pa[i][j];
		}
	}
	dest->num_pa = num_pa;
	dest->pa = pa;
}

void free_parents(model* M) {
	free(M->num_pa);
	free_int_matrix(M->pa, M->n);
}

void print_model(model M, int* theta_types) {
	printf("\n\n-- Model instance: --\n\n");
	printf("Number of events = %d\n", M.n);
	printf("Cover Relations: \n");
	print_int_matrix(M.P, M.n, M.n);

	printf("\n\nnum_parents: parent sets\n\n");
	for (int i = 0; i< M.n; i++) {
		printf("event %d - %d: ", i, M.num_pa[i]);
		if (M.num_pa[i] > 1) {
			if (theta_types[i] == 0) printf("AND ");
			if (theta_types[i] == 1) printf("OR  ");
			if (theta_types[i] == 2) printf("XOR ");
		} else {
			printf("--- ");
		}
		for(int j = 0; j < M.num_pa[i]; j++) {
			printf(" %d", M.pa[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

bool equal_models(model* A, model* B) {
	if (A->n != B->n) {
		printf("n off.\n");
		return false;
	}

	int n = A->n;
	for (int i = 0; i< n; i++) {
		for (int j = 0; j< n; j++) {
			if (A->P[i][j] != B->P[i][j]) {
				printf("posets off.\n");
				return false;
			}
		}
	}

	return true;
}


#endif
