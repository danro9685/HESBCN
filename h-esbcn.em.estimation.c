#include <stdio.h>
#include <stdlib.h>
#include "h-esbcn.h"
#include "model.h"
#include "queue.h"
#include "markov.h"

/** 
 * Free memory used in program. 
 */
void cleanup(model M, data* D, int N_u, int* theta_types_plus, double* lambdas) {
	free_int_matrix(GENOTYPE, pow2(M.n));
	free_data_array(D, N_u, M.n);
	free_model(M, true);
	free_int_array(theta_types_plus);
	free_double_array(lambdas);
}


int main(int argc, char* argv[]) {
	// get commandline args
	char poset[10];
	char thetas[10];
	char output[100];
	char data_stem[100];
	get_em_inputs(argc, argv, poset, thetas, output, data_stem);

	model M;
	int N_u; // number of unique genotypes (patterns) observed
	data* D; // Stores read-in data set information
	process_patterns(data_stem, &M, &N_u, &D);	//read and process patterns, generate data set D. 
	get_poset(&M, poset);
	int* theta_types = get_theta_types(thetas, M.n);

	//populate the genotype global so that each genotype is only computed once. 
	precompute_binary(M.n);

	// dummy initializations to work with previous code. 
	M.num_pa = malloc(sizeof(int));
	M.pa = get_int_matrix(M.n, sizeof(int));

	set_theta_types(&M, theta_types, theta_types);


	// initialize EM, run.
	int* theta_types_plus = get_int_array(M.n+ 1);
	for (int i = 0; i < M.n; i++) {
		theta_types_plus[i] = theta_types[i];
	}
	free_int_array(theta_types);
	theta_types_plus[M.n] = -1;

	add_one_event(&M, D, N_u, theta_types_plus);

	// EM algorithm to guess lambdas. 
	double* lambdas = get_double_array(M.n);
	em_initialization(&M, theta_types_plus, D, N_u, lambdas); //  

	expectation_maximization(&M, theta_types_plus, D, N_u, lambdas, output);

	cleanup(M, D, N_u, theta_types_plus, lambdas);
}
