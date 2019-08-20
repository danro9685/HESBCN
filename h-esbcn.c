#include "h-esbcn.h"
#include "mcmc.h"
#include "markov.h"
#include <math.h>

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

int main(int argc, char** argv) {
	char filestem[512];
	char output[512];
	char reg_scheme[10]; // default value

	// initialize some parameters 
	unsigned int seed;
	int number_samples, record_ith, burn_in;
	initialize_MH_params(&seed, &number_samples, &record_ith, &burn_in, reg_scheme);

	check_file_inputs(argc, argv, filestem, output, &number_samples, &seed, reg_scheme);	//process args, get input filename.
	model M;

	int N_u; // number of unique genotypes (patterns) observed
	data* D; // Stores read-in data set information
	process_patterns(filestem, &M, &N_u, &D);	//read and process patterns, generate data set D. 


	//populate the genotype global so that each genotype is only computed once. 
	precompute_binary(M.n);

	// initialize these two fields for future use
	M.P = get_int_matrix(M.n, M.n);

	// initialize these fields for set_theta_types
	M.num_pa = get_int_array(M.n); // create an array that stores the sizes of parents. 
	M.pa = (int**) malloc(M.n * sizeof(int*)); // create an array to hold arrays of integers
	for (int i = 0; i< M.n; i++) {
		M.pa[i] = (int*) malloc(4); // 4 is machine friendly - in reality, this initialization is for use with the following few function calls. 
	}
	int* theta_types = get_int_array(M.n);
	for (int i = 0; i< M.n; i++) { // with the empty poset, all events guaranteed to have 0 parents. 
		theta_types[i] = -1;
	}


	set_theta_types(&M, theta_types, theta_types);
	new_total_ordering(&M, theta_types, D, N_u); 

	metropolis_hastings(&M, D, N_u, burn_in, number_samples, record_ith, seed, output, theta_types, reg_scheme);

	set_theta_types(&M, theta_types, theta_types);

	/* Notes: 
	 * Diagnosis of cancer will now be treated as an event, stored in the last index.
	 * Add one event for diagnosis
	 */

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
	return 0;
}
