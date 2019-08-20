#ifndef MARKOV_H
#define MARKOV_H

#include <stdbool.h>
#include <math.h>
#include <float.h>

#include "model.h"
#include "matrix.h"
#include "h-esbcn.h"
#include "a-queue.h"
#include "likelihood_computation.c"

#define EM_MAX_ITER 10000
#define EM_ACCURACY 1e-8

int* VALIDS; // list of valid genotypes. Right now, repurposed to be 2^n size and 1 if index valid, 0 otherwise. Cuts out need for computation to check compatibility
int* ABS_TO_VALID; // array where index i stores the position of genotype i in the smaller valid list in the model.  
int* LIFT; // m lengthed array where LIFT[i] gives the compressed integer representation of the parent set lift. 

int NUM_END_STATES; // Number of end states, which are averaged to determine expected times for each event and likelihood.
int* END_STATES; // array of length NUM_END_STATES containing the indexes of all end states. We define an end state as a genotype w/ no children.


void get_em_inputs(int argc, char** argv, char* poset, char* thetas, char* output, char* data_stem) {
	int c = 0;
	bool poset_flag = false;
	bool theta_flag = false;
	bool output_flag = false;
	bool input_flag = false;
	bool help_flag = false;
	static struct option long_options[] =
        {
          {"poset",   	  		required_argument, 0, 'p'},
          {"thetas",   	  		required_argument, 0, 't'},
          {"output",  				required_argument, 0, 'o'},
          {"data",  					required_argument, 0, 'd'},
          {"help", 						no_argument, 			 0, 'h'},
          {NULL, 0, NULL,0,}
        };
	while ((c = getopt_long(argc, argv, "d:p:t:o:h", long_options, NULL)) != -1 ) {
		switch(c) {			
		case 'd': // get input filename, satisfy boolean check for provided input. 
				strcpy(data_stem, optarg);
				input_flag = true;
				break;
			case 'p':
				strcpy(poset, optarg);
				poset_flag = true;
				break;
			case 't': // get input filename, satisfy boolean check for provided input. 
				strcpy(thetas, optarg);
				theta_flag = true;
				break;
			case 'o':
				strcpy(output, optarg);
				output_flag = true;
				break;
			case 'h':
				help_flag = true;
				break;
		}
	}

	if (help_flag) {
		printf("Usage: ./h-esbcn.em.estimation (-(option) (argument) )*\n");
		printf("where input option (-i or --input), output option(-o or --output), and type (-t or --type) are required.\n");
		printf("Options include:\n");
		printf("-p | --poset \tPath to input file with inferred poset, comma separated. Required argument.\n");
		printf("-t | --thetas \tPath to input file with theta types. Required argument.\n");
		printf("-o | --output \tPath to output file to write results to.\n");
		printf("-d | --dataset \tPath to input file with sample data. Please include file extension as well. ex: test1.txt\n");
		printf("-h | --help \tPrint help information for program.\n");
		exit(1);
	}
	if (!input_flag) {
		fprintf(stderr, "Error: No input file specified. Please provide the filename (without filetype) of a file containing formatted sample data using the \"-d\" or \"--dataset\" options.\n");
		exit(1);
	}
	if (!poset_flag) {
		fprintf(stderr, "Error: No input file specified. Please provide the filename of a file containing formatted conditional probability table using the \"-i\" or \"--input\" options.\n");
		exit(1);
	}
	if (!theta_flag) {
		fprintf(stderr, "Error: No theta-type file specified. Please provide the filename of a file containing comma separated theta types with the -t option.\n");
		exit(1);
	}
	if (!output_flag) {
		fprintf(stderr, "Error: No output path specified. Please provide the output destination using the \"-o\" or \"--output\" options.\n");
		exit(1);
	}
}


int* get_theta_types(char* thetas, int n) {
	FILE* f = open_file(thetas, "r");
	char* line = NULL;
	size_t temp;
	getline(&line, &temp, f);
	char* token = strtok(line, ",");
	int* theta_types = get_int_array(n);
	for (int i = 0; i < n; i++) {
		theta_types[i] = atoi(token);
		token = strtok(NULL, ",");
	}
	fclose(f);
	return theta_types;
}


void get_poset(model* M, char* poset) {
	FILE* f = open_file(poset, "r");
	int n = M->n;
	size_t temp = 0;
	char* line = NULL; 
	M->P = get_int_matrix(n, n);
	for (int i = 0; i < n; i++) {
		getline(&line, &temp, f);
		for (int j = 0; j < n; j++) {
			M->P[i][j] = 0;
			if (line[2 * j] == '1') {
				M->P[i][j] = 1;
			}
		}
	}
	fclose(f);
}


// Add 1 event, indicating whether genotype has been diagnosed as cancer. This creates 2 forms for each genotype - diagnosed and undiganosed. 
// Simplified slightly by placing the diagnosis bit in the least significant bit, such that all evens are undiagnosed and odds are diagnosed. 
void add_one_event(model* M, data* D, int num_unique, int* theta_types) {
	int n = M->n;
	int new_n = M->n + 1;

 	// new matrix to include binary event. Is this actually needed?
	int** new_poset = get_int_matrix(new_n, new_n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			new_poset[i][j] = M->P[i][j];
		}
	}
	free_int_matrix(M->P, n);
	M->P = new_poset;

	free_parents(M);
	M->n = new_n;

	// initialize these fields for set_theta_types
	M->num_pa = get_int_array(new_n); // create an array that stores the sizes of parents. 
	M->pa = (int**) malloc(new_n * sizeof(int*)); // create an array to hold arrays of integers
	for (int i = 0; i< new_n; i++) {
		M->pa[i] = (int*) malloc(4); // 4 is machine friendly - in reality, this initialization is for use with the following few function calls. 
	}
	set_theta_types(M, theta_types, theta_types);

	// Extend GENOTYPE to accomodate the new event. 
	free_int_matrix(GENOTYPE, pow2(n));
	precompute_binary(new_n); // set up new set of all binaries


	for (int i = 0; i < num_unique; i++) {
		int* observed = D[i].g;
		int* augmented = get_int_array(new_n);
		for (int j = 0; j < n; j++) {
			augmented[j] = observed[j];
		}
		augmented[n] = 1;
		D[i].g = augmented;
		free_int_array(observed);
	}
}


// BFS method of building the list of valid genotypes. Gives a particular random ordering, for book-keeping purposes. 
void build_valids(model* M, int* theta_types) {
	int num_events = M->n;
	int lattice_size = pow2(num_events);
	vector* v = new_vector(); // fill with correctly sized and bfs ordered genes.
	VALIDS = get_int_array(lattice_size);

	queue* q = new_queue();
	enqueue(q, 0);
	while (!is_empty(q)) {
		int curr = dequeue(q);
		if (!VALIDS[curr]) {
			VALIDS[curr] = 1;
			push_back(v, curr);

			int* g = GENOTYPE[curr];
			for (int i = 0; i < num_events; i++) {
				if (!g[i]) {
					g[i] = 1;
					int compressed = binary_compression(g, num_events);
					if (compatible_with_poset(compressed, M, theta_types)) {
						enqueue(q, compressed);
					}
					g[i] = 0;
				}
			}
		}
	}

	free_queue(q);
	M->m = v->size;
	M->valids = free_vector(v, false);

	ABS_TO_VALID = get_int_array(pow2(num_events));
	for (int i = 0; i < M->m; i++) {
		int genotype = M->valids[i];
		ABS_TO_VALID[genotype] = i;
	}
}

// Find the children of each genotype, making note of the indices where they differ from the parents. 
// Note: Repeats work - if possible, try to merge this method with build_valids.
void find_children(model* M, int* theta_types) {
	int lattice_size = M->m; 
	int num_events = M->n;

	// find children
	M->num_ch = get_int_array(lattice_size); // heap allocated and pre-initialized to 0. Record number of children for each genotype in lattice. 
	M->ch_diff = malloc(lattice_size * sizeof(int*));
	if (M->ch_diff == NULL) {
		fprintf(stderr, "Error: out of memory.\n");
		exit(1);
	} 

	vector* v;
	// find valid children for each genotype
	for (int i = 0; i< lattice_size; i++) {
		v = new_vector();
		int parent = M->valids[i];
		int* parent_g = GENOTYPE[parent];

		// compile a list of valid and possible children, given the poset.
		for (int j = 0; j < num_events; j++) {
			if (!parent_g[j]) {
				parent_g[j] = 1;
				int compressed = binary_compression(parent_g, num_events);
				if (VALIDS[compressed]) {
					push_back(v, j);
				}

				parent_g[j] = 0;
			}
		}
		M->num_ch[i] = v->size;
		M->ch_diff[i] = free_vector(v, false);
	}
}


// Quick way to generate a list of valid parent genotypes for each valid genotype. Suboptimal, but works for now. 
void find_parents(model* M) {
	int n = M->n;
	int m = M->m;
	vector** parents = malloc(sizeof(vector*) * m);
	vector** indices = malloc(sizeof(vector*) * m);
	if (parents == NULL || indices == NULL) {
		fprintf(stderr, "Out of memory: exiting.\n");
		exit(1);		
	}

	for (int i = 0; i < m; i++) {
		parents[i] = new_vector();	
		indices[i] = new_vector();
		if (parents[i] == NULL || indices[i] == NULL) {
			fprintf(stderr, "Out of memory: exiting.\n");
			exit(1);
		}		
	}

	for (int i = 0; i < m; i++) { // loop through lattice
		int parent = M->valids[i];
		int* g = GENOTYPE[parent];
		for (int j = 0; j < M->num_ch[i]; j++) {
			// push the parent into all of the children's parent vectors
			int diff_index = M->ch_diff[i][j];

			g[diff_index] = 1; // get child genotype
			int child = binary_compression(g, n);
			int child_index = ABS_TO_VALID[child];
			push_back(parents[child_index], parent);
			push_back(indices[child_index], diff_index);
			g[diff_index] = 0; // back to parent genotype
		}
	}

	M->g_num_pa = get_int_array(m);
	M->g_pa = malloc(sizeof(int*) * m);
	M->g_pa_diff = malloc(sizeof(int*) * m);

	for (int i = 0; i < m; i++) {
		M->g_num_pa[i] = parents[i]->size;
		M->g_pa[i] = free_vector(parents[i], false);
		M->g_pa_diff[i] = free_vector(indices[i], false);
	}

	free(parents);
	free(indices);
}

// Find genotypic end states (genotypes that cannot acquire more mutations). 
// We use genotypic completion over functional completion, simply to get as much information as possible about all events. 
void define_end_state(model* M) {
	int m = M->m;
	int* num_ch = M->num_ch;
	vector* end_states = new_vector();
	for (int i = 0; i < m; i++) {
		if (num_ch[i] == 0) { // if there are no genotypic children, then genotype is an end state.
			push_back(end_states, i);
		}
	}

	NUM_END_STATES = end_states->size;
	END_STATES = free_vector(end_states, false);
}


// Initialize fields needed in Data struct for EM. 
void data_preprocess(model* M, int* theta_types, data* D, int num_unique) {
	int n = M->n;
	int m = M->m;
	int* valids = M->valids;
	int** original_poset = M->P;
	for (int  i = 0; i < num_unique; i++) {
		int* observed = D[i].g;
		int compressed = binary_compression(observed, n);
		D[i].is_compatible = (VALIDS[compressed] == 1); // is_compatible set to true if compressed is 
		D[i].subset = get_int_array(m);

		D[i].P = get_int_matrix(n, n);

		// For each data point, layer diagnostic event into a copy of the poset, assuming that diagnosis happened as a result of the genotype. 
		for (int j = 0; j < n - 1; j++) {
			for (int k = 0; k < n - 1; k++) {
				D[i].P[j][k] = M->P[j][k];
			}
			D[i].P[j][n - 1] = observed[j];
			D[i].P[n - 1][j] = (observed[j] + 1) % 2; 
		}

		// Modification over original induced refinements - invalid over induced refinement only if observed[i] = 1 and g[j] = 0
		M->P = D[i].P;
		int num_compatible = 1;
		D[i].subset[0] = 1;
		for (int j = 1; j < m; j++) {
			D[i].subset[j] = 0;
			if (compatible_with_poset(valids[j], M, theta_types)) {		
				D[i].subset[j] = 1;
				num_compatible += 1;
			}
		}
		D[i].t = get_double_array(n);
	}

	M->P = original_poset;
}

// Output frequencies of events relative to their parent sets in the data. Useful as a sanity check.
void freqs(model* M, int* theta_types, data* D, int num_unique) {
	int n = M->n;
	for (int i = 0; i < n; i++) {
		printf("event %d\t", i);
		int event_count = 0;
		int parent_count = 0;
		if (M->num_pa[i] == 0) {
			for (int j = 0; j < num_unique; j++) {
				if (D[j].g[i]) event_count += D[j].count;
			}
			printf("Observed %d of %d\n", event_count, NUM_SAMPLES_READ);
		} else { // conjunctive
			for (int j = 0; j < num_unique; j++) {
				if (D[j].g[i]) {
					event_count += D[j].count;
				}
				
				if (check_parent_set(M, theta_types, D[j].g, i)) {
					parent_count += D[j].count;
				}
			}
			printf("Observed %d of %d\n", event_count, parent_count);
		}
	}
	printf("\n\n");
}

// Generate parent "lift genotypes", indicating whether each event in a given genotype has its parent set is fulfilled.
// To save space, these are also stored as ints, and array representations are available via the GENOTYPE global.
void generate_parent_lift(model* M, int* theta_types) {
	int m = M->m;
	int n = M->n;
	LIFT = get_int_array(m);
	int g[n];

	for (int i = 0; i < m; i++) {
		int* genotype = GENOTYPE[M->valids[i]];
		for (int j = 0; j < n; j++) { // generate binary array indicating whether event j's parent set is fulfilled
			if (check_parent_set(M, theta_types, genotype, j)) {
				g[j] = 1;
			} else {
				g[j] = 0;
			}
		}

		LIFT[i] = binary_compression(g, n);
	}
}


// Compute initial lambda values based on dataset frequencies, which are used for probability computations. 
// We employ add-one smoothing in num_primed to ensure that num_primed > num_mutated (p != 1) and we don't get a lambdas[i] = 1/0 situation.
void compute_lambdas(model* M, int* theta_types, double* lambdas, data* D, int num_unique) {
	int n = M->n - 1; // last event is diagnostic event
	for (int i = 0; i < n; i++) { // tally up counts for each genotype where parent set is fulfilled.
		int num_mutated = 0; // num with parent set fulfilled and mutated. 
		int num_primed = 1; // num with parent set fulfilled
		for (int j = 0; j < num_unique; j++) {
			if (check_parent_set(M, theta_types, D[j].g, i)) {
				num_mutated += D[j].count * D[j].g[i];
				num_primed += D[j].count;
			}
		}
		double p = (double) num_mutated / num_primed;
		lambdas[i] = p / (1.0 - p);
	}
	lambdas[n] = 1.0;
}

// Precompute marginal lambdas for each genotype. lambda_exit used later to estimate probability that a genotype obtains a given child node.
// Intuitively, lambda_exit values are the marginal lambdas for each genotype (combined lambdas for the possible mutations, given a state). 
void compute_lambda_exit(model* M, double* lambdas, double* lambda_exit) {
	int m = M->m;
	for (int i = 0; i < m; i++) {
		lambda_exit[i] = 0.0;
		for (int j = 0; j < M->num_ch[i]; j++) {
			int index = M->ch_diff[i][j];
			lambda_exit[i] += lambdas[index];
		}
	}
}

// Compute probabilities of each genotype, based on existing lambda and lambda_exit values. 
// lambdas[diff_index] / lambda_exit[index] = probability of obtaining event, given genotype index.
void compute_prob(model* M, double* lambdas, double* lambda_exit, double* prob, int* subset) {
	int m = M->m;

	prob[0] = 1.0; // 0 contains wildtype
	for (int i = 1; i < m; i++) {
		prob[i] = 0.0;
		if (subset[i]) {
			for (int j = 0; j < M->g_num_pa[i]; j++) {
				int index = ABS_TO_VALID[M->g_pa[i][j]];
				
				if (subset[index]) { 
					int diff_index = M->g_pa_diff[i][j];
					prob[i] += (lambdas[diff_index] / lambda_exit[index]) * prob[index];	
				}
			}
		}
	}
}


// event is locked out if it's in an XOR parent relationship but is not the observed parent while child is present. 
// @param index: refers to index of genotype in M->valids. Use ABS_TO_VALID to get the conversion from integer representation.
bool locked_out(model* M, int* theta_types, int index, int event) {
	int* g = GENOTYPE[M->valids[index]];
	if (g[event]) return false; // if event is already on, it's not locked out.
	
	// event must be off at this point. 
	int n = M->n;
	for (int i = 0; i < n; i++) {
		if (M->P[event][i] && theta_types[i] == 2 && g[i]) { // if event is a parent of i, i's parents are XOR, and i present -> event can't manifest. 
			return true;
		}
	}
	return false;
}


// Compute expected waiting times for each event (time from parent set fulfillment -> mutation), taking into account diagnostic time event in an observed data point. 
// Computations based on derivations in Beerenwinkel 2007.
void compute_exp(model* M, int* theta_types, double* lambdas, double* lambda_exit, double* prob, int* subset, double** exp) {
	int n = M->n;
	int m = M->m;

	for (int pos = 0; pos < n; pos++) {
		for (int i = 0; i < m; i++) {
			exp[pos][i] = 0.0;
			if (subset[i]) {
				for (int c=0; c<M->g_num_pa[i]; c++) {
					int parent_index = ABS_TO_VALID[M->g_pa[i][c]];  // index of c-th parent (in J_P) of i
					if (subset[parent_index]) {
						int diff_index = M->g_pa_diff[i][c];  // index of differing event with parent
						int* g = GENOTYPE[M->valids[parent_index]]; // parent genotype
						int* l_g = GENOTYPE[LIFT[parent_index]]; // get the parent lift array

						exp[pos][i] += (lambdas[diff_index] / lambda_exit[parent_index]) * exp[pos][parent_index];
						if (l_g[pos] && !g[pos]) { 
							exp[pos][i] += (lambdas[diff_index] * prob[parent_index]) /(lambda_exit[parent_index] * lambda_exit[parent_index]);
						}
					}
				}
			}
		}
	}
}

// Assign lambda MLE based on expected times and observed datapoints. 
void lambdas_MLE(const int n, data* D, const int num_unique, double* lambdas, const int num_compatible) {
	for (int i = 0; i < n; i++) {
		double sum = 0.0;
		for (int j = 0; j < num_unique; j++) {
			if (D[j].is_compatible) {
				sum += D[j].count * D[j].t[i];
			}
		}
		lambdas[i] = (double) num_compatible / sum;
	}
}

// Take norm of normalized difference between x and y, w.r.t. x. 
double norm_diff(double* x, double* y, const int n) {
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		double summand = (x[i] - y[i]) / x[i];
		sum += summand * summand;
	}
	return sum;
}

// Gets number of datapoints compatible with the current poset. 
int get_num_compatible(data* D, int num_unique) {
	int num_compatible = 0;
	for (int i = 0; i < num_unique; i++) {
		if (D[i].is_compatible) num_compatible += D[i].count;
	}
	return num_compatible;
}


// Initial EM for a warm start on the nested EM. Assumes zero noise and no censoring for simplicity. 
// Needed to run the EM_epsilon in the first iteration of the main nested EM loop.
double EM_lambdas(model* M, int* theta_types, data* D, int num_unique, double* lambdas, double* percent_compatible, double eps) {
	double loglik = 0.0;
	double loglik_new = 0.0;
	double log_percent = log(*percent_compatible);

	int m = M->m;
	int n = M->n;
	double S = lambdas[n - 1];

	double* lambdas_new = get_double_array(n); // one per event
	double* lambda_exit = get_double_array(m); // one per genotype
	double* prob = get_double_array(m); // one per genotype
	double** exp = get_double_matrix(n, m);

	double* lambdas_best = get_double_array(n); // save the best ones

	// precompute this value, so that we don't re-compute it at each iteration. 
	int num_compatible = get_num_compatible(D, num_unique);
	int iter = 0;

	while ((norm_diff(lambdas, lambdas_new, n) > EM_ACCURACY || (iter < 10)) && iter < 10) { 
		if (iter >  10 && loglik_new < loglik) { // shouldn't happen, but it happens. could have to do with our non-conjunctive model. 
			fprintf(stderr, "Error: likelihood decreasing in EM_lambdas. loglik_new = %f, loglik = %f. Exiting.\n", loglik_new, loglik);
			print_double_array(lambdas, n);
			print_double_array(lambdas_new, n);
			break;
		}
		// update parameters
		loglik = loglik_new;
		loglik_new = 0.0;
		for (int i = 0; i < n; i++) {
			lambdas_new[i] = lambdas[i];
		}

		// == E step ==
		compute_lambda_exit(M, lambdas, lambda_exit); // update lambda_exit values

		for (int i = 0; i < num_unique; i++) {
			if (D[i].is_compatible) {
				compute_prob(M, lambdas, lambda_exit, prob, D[i].subset);

				compute_exp(M, theta_types, lambdas, lambda_exit, prob, D[i].subset, exp);

				for (int j = 0; j < NUM_END_STATES; j++) {
					int index = END_STATES[j];
					if (D[i].subset[index]) { // if end state is permissible 
						loglik_new += D[i].count * (log_percent + log(prob[index]));
					}
				}

				for (int j = 0; j < n; j++) {
					double sum = 0.0;
					for (int k = 0; k < NUM_END_STATES; k++) {
						int index = END_STATES[k];
						sum += exp[j][index] / prob[index];
					}
					D[i].t[j] = sum / NUM_END_STATES;
				}
			}
		}

		// == M step ==
		lambdas_MLE(n, D, num_unique, lambdas, num_compatible);

		double fac = S / lambdas[n - 1];
		for (int i = 0; i < n; i++) {
			lambdas[i] *= fac;
		}
		iter++;
	}

	free_double_array(lambdas_new);
	free_double_array(lambda_exit);
	free_double_array(prob);
	free_double_matrix(exp, n);
	free_double_array(lambdas_best);

	return loglik_new;
}

// Estimates lambda values for a warm start to nested EM. 
double estimate_parameters(model* M, int* theta_types, data* D, int num_unique, double* lambdas, double* percent_compatible) {
	int n = M->n;
	int num_compatible = 0;
	for (int i = 0; i < num_unique; i++) {
		if (D[i].is_compatible) num_compatible += D[i].count;
	}

	double fraction_compatible = (double) num_compatible / NUM_SAMPLES_READ;
	*percent_compatible = fraction_compatible;

	double q_eps = 1.0 / (pow2(M->n) - M->m); // the +1 includes the diagnosis time. Do we need it?
	
	double loglik_incompatible = 0.0;
	if (NUM_SAMPLES_READ > num_compatible) {
		loglik_incompatible = (NUM_SAMPLES_READ - num_compatible) * (log(1.0 - fraction_compatible) + log(q_eps));
	}

	double* lambdas_em = get_double_array(n);
	compute_lambdas(M, theta_types, lambdas_em, D, num_unique);
	double loglik = 0;
	double loglik_best = -DBL_MAX;

	// EM for lambdas
	loglik = EM_lambdas(M, theta_types, D, num_unique, lambdas_em, percent_compatible, q_eps) + loglik_incompatible;
	
	// if more likely than current lambdas, set to new. 
	// Previously ran a few iterations of EM here, but decided to just do 1 for a warm start and push iterations to main loop. 
	if (loglik > loglik_best) {
		loglik_best = loglik;
		memcpy(lambdas, lambdas_em, sizeof(double) * n);
	}

	free_double_array(lambdas_em);
	return loglik_best;
}


// Initializes structures in Model, Data, lambdas, and some globals. 
void em_initialization(model* M, int* theta_types, data* D, int num_unique, double* lambdas) {
	// generate list of all valid genotypes.
	build_valids(M, theta_types);
	find_children(M, theta_types);

	find_parents(M);

	// initialize data time fields and construct sublattices
	data_preprocess(M, theta_types, D, num_unique);

	generate_parent_lift(M, theta_types);
	
	define_end_state(M);

	// get lambda values
	double percent_compatible;
	estimate_parameters(M, theta_types, D, num_unique, lambdas, &percent_compatible);

	// output parent/child frequencies
	freqs(M, theta_types, D, num_unique);
}


// Uses BFS ordering for a DP approach to computing probabilities for ALL valid genotypes (no subset restriction).
// Note: we assume in this model that there are no events that influence diagnosis time (last event), but that's not true is it?
// Most obviously, diagnosing a wildtype as cancer should have a probability of near 0, while genotypes with other mutations
// should have higher diagnostic chances. 
void compute_all_prob(model* M, double* lambdas, double* lambda_exit, double* prob) {
	int m = M->m;
	int n = M->n;

	prob[0] = 1.0;
	double prob_diag = lambdas[n-1];

	int index_1 = ABS_TO_VALID[1]; // 1 is the genotype for no events but diagnosed. 
	prob[index_1] = prob[0] * prob_diag / lambda_exit[0];


	for (int i = 1; i < m; i++) {
		int genotype = M->valids[i];
		if (genotype % 2 == 0) { // undiagnosed event. Since diagnosis is put in the least significant bit, all undiagnosed genotypes are even. 
			prob[i] = 0.0;
			for (int j = 0; j < M->g_num_pa[i]; j++) {
				int parent_index = ABS_TO_VALID[M->g_pa[i][j]];
				int diff_index = M->g_pa_diff[i][j];
				prob[i] += (lambdas[diff_index] / lambda_exit[parent_index]) * prob[parent_index];
			}
			int diagnosed_index = ABS_TO_VALID[genotype + 1];
			prob[diagnosed_index] = prob[i] * prob_diag / lambda_exit[i];

		}
	}
}

// Compute conditional probabilities of the observed dataset, factoring in the provided epsilon. 
void compute_condprob(model* M, data* D, int num_unique, double* prob, double epsilon, double** condprob, double* marginals) {
	int n = M->n;
	int m = M->m;

	for (int i = 0; i < num_unique; i++) {
		double marginal = 0.0;
		int compressed = binary_compression(D[i].g, n);
		for (int j = 1; j < m; j++) {
			condprob[j][i] = 0.0;
			if (M->valids[j] % 2) { // only observed patterns, since this function computes conditional probabilities of the observed dataset. 
				int d = hamming_distance(M->valids[j], compressed, n, NULL);
				condprob[j][i] = prob[j] * power(epsilon, d) * power(1 - epsilon, n - d);
				marginal += condprob[j][i];
			}
		}

		for (int j = 1; j < m; j++) {
			if (M->valids[j] % 2) {
				condprob[j][i] /= marginal;
			}
		}

		marginals[i] = marginal;
	}
}

// EM to find maximum likelihood error rate (epsilon), given some set of lambda values. 
double EM_epsilon(model* M, int* theta_types, data* D, int num_unique, double* prob, double** condprob, double* epsilon, double* marginals) {
	double loglik = 0;
	double loglik_new = 0;
	double delta_epsilon = 0.0;
	int iter = 0;
	int n = M->n;
	int m = M->m;
	int em_iter = 10;

	while (iter < 10 && (iter < em_iter / 5 || fabs(delta_epsilon) / *epsilon > 1e-4)) {
		if ((iter > em_iter / 2) && (loglik_new < loglik)) {
			fprintf(stderr, "Warning in EM_epsilon: Likelihood is decreasing.\n");
			printf("loglik_new = %g, loglik = %g\n", loglik_new, loglik);
			break;
		}

		loglik = loglik_new;
		loglik_new = 0;

		compute_condprob(M, D, num_unique, prob, *epsilon, condprob, marginals);

		double H = 0.0; // entropy 
		for (int i = 0; i < num_unique; i++) {
			int compressed = binary_compression(D[i].g, n);
			for (int j = 1; j < m; j++) {
				if (M->valids[j] % 2) {
					int d = hamming_distance(M->valids[j], compressed, n, NULL);
					H += d * condprob[j][i] * D[i].count;
				}
			}
			loglik_new += log(marginals[i]) * D[i].count;
		}
		H /= (NUM_SAMPLES_READ * n);

		delta_epsilon = H - *epsilon;
		*epsilon = H;

		iter++;
	}
	return loglik_new;
}


// Return true if genotype at x comes after y. Strict toggles '='
// Tried to remove assumption of ordering, considering only whether genotypes come before one another in time ( if possible ). 
bool is_after(model* M, int* theta_types, int x, int y, bool strict) {
	// return true if x has more parent sets fulfilled, along the *same lineage path*. 
	if (x == y) return !strict;
	int n = M->n;

	int* x_lift = GENOTYPE[LIFT[x]];
	int* y_lift = GENOTYPE[LIFT[y]];
	int* g_x = GENOTYPE[M->valids[x]];
	int* g_y = GENOTYPE[M->valids[y]];

	// if y -> x, then x has all non-xor parent sets in y. 
	for (int i = 0; i < n; i++) {
		// y has a mutation, then x must have it too.
		if (g_y[i] && !g_x[i]) return false;

		// at this point, either both have mutation or mutation not in y.
		if (!g_y[i] && y_lift[i]) { // if parent set is AND or OR, then must be in x also.
			if (theta_types[i] <= 0) {
				if (!x_lift[i]) return false;
			} else if (theta_types[i] == 1) { // check that all parents in y are also in x
				if (!x_lift[i]) return false;
				for (int j = 0; j < M->num_pa[i]; j++) {
					if (g_y[j] && !g_x[j]) return false;
				}		
			} else { // if parent set is XOR, then is only ok if x has the mutation in y AND another parent.
				// if here, then i is primed but not mutated in y.
				// if i is primed in x, it must also be with the same parent.
				int* parents = M->pa[i];
				if (x_lift[i]) {
					for (int j = 0; j < M->num_pa[i]; j++) {
						int parent = parents[j];
						if (g_y[parent] != g_x[parent]) return false;
					}
				} else { // if i is NOT primed in x, it must have too many parents.
					int num_x_parents = 0;
					for (int j = 0; j < M->num_pa[i]; j++) {
						if (g_x[parents[j]]) num_x_parents++;
					}
					if (num_x_parents <= 1) return false;
				}
			}

		}
	}
	return true;
}


// Compute the prob[G|X]. Probability of each genotype, given the noisy observation X. 
// censprob is a m x m matrix, where censprob[i][j] is the probability of genotype i given genotype j.
// Should be unchanged by OR/XOR.
void compute_censored_prob(model* M, int* theta_types, double* lambdas, double* lambda_exit, double* prob, double** censprob) {
	int m = M->m;

	// clear out the censprob array.
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			censprob[i][j] = 0;
		}
	}

	for (int i = 0; i < m; i++) {
		if (M->valids[i] % 2) { // observed (diagnosed) patterns only
			censprob[i][i] = prob[i];
			for (int j = 0; j < m; j++) {
				if (is_after(M, theta_types, j, i, true)) { // if j "comes after" i
					for (int k = 0; k < M->g_num_pa[j]; k++) {
						int parent = M->g_pa[j][k];
						int parent_index = ABS_TO_VALID[parent]; // kth parent index. 

						if (is_after(M, theta_types, parent_index, i, false)) { // if j's parent is before i
							int diff_index = M->g_pa_diff[j][k];
							censprob[i][j] += lambdas[diff_index] / lambda_exit[parent_index] * censprob[i][parent_index];
						}
					}
				}
			}
		}
	}
}

// compute expected waiting times for each event, factoring in censored probabities. 
// Assumes that conditional expectations shouldn't change from the conjunctive model, since we use parent set fulfillment as the condition
// to check whether a genotype contributes to each event's Expected time. 
void compute_all_exp(model* M, int* theta_types, double* lambdas, double* lambda_exit, double* prob, double** exp, double** censprob) {
	int n = M->n;
	int m = M->m;
	double*** censexp = get_double_cube(n, m, m);

	// clear out previous values for safetys
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			exp[i][j] = 0.0;
		}
	}

	// calculate uncensored expectation.
	for (int pos = 0; pos < n; pos++) {
		for (int i=0; i<m; i++) {// First calculate events BEFORE censoring. i.e. unobserved
			exp[pos][i] = 0.0;
			for (int c = 0; c < M->g_num_pa[i]; c++) {
				int parent = M->g_pa[i][c];  // index of c-th parent (in J_P) of i
				int parent_index = ABS_TO_VALID[parent];
				if (parent % 2 == 0) { // BEFORE
					int diff_index = M->g_pa_diff[i][c];  // index of differing event
					int* g = GENOTYPE[parent];

					exp[pos][i] += (lambdas[diff_index] / lambda_exit[parent_index]) * exp[pos][parent_index];
					int* l_g = GENOTYPE[LIFT[parent_index]];
					if (l_g[pos] && !g[pos]) {
						exp[pos][i] += (lambdas[diff_index] * prob[parent_index]) / (lambda_exit[parent_index] * lambda_exit[parent_index]);
					}
				}
			}
		}

		// ------ Second half --------
		// Now calculate the events AFTER censoring. i.e. observation
		for (int i = 0; i < m; i++) {
			if (M->valids[i] % 2) { 
				censexp[pos][i][i] = exp[pos][i];
				for (int l = 0; l < m; l++) {
					if (is_after(M, theta_types, l, i, true)) {
						for (int c = 0; c < M->g_num_pa[l]; c++) {
							int parent = M->g_pa[l][c];  // index of c-th parent (in J_P) of i
							int parent_index = ABS_TO_VALID[parent];
							if (is_after(M, theta_types, parent_index, i, false)) {
								int j = M->g_pa_diff[l][c];  // index of differing event
								int* g = GENOTYPE[parent];

								censexp[pos][i][l] += (lambdas[j] / lambda_exit[parent_index]) * censexp[pos][i][parent_index];
								int* l_g = GENOTYPE[LIFT[parent_index]];
								if (l_g[pos] && !g[pos] && !locked_out(M, theta_types, parent_index, pos)) {
									censexp[pos][i][l] += (lambdas[j] * censprob[i][parent_index]) / (lambda_exit[parent_index] * lambda_exit[parent_index]);
								}
							}
						}
					}
				}

				// Set the expected time for an event to the average of all end states. 
				double sum = 0;
				for (int j = 0; j < NUM_END_STATES; j++) {
					int index = END_STATES[j];
					sum += censexp[pos][i][index];
				}
				exp[pos][i] = sum / NUM_END_STATES;
			}
		}
	}
	free_double_cube(censexp, n, m);
}


// check if two genotypes (index in model) have same lineage.
// We are interested in making comparisons between two genotypes ONLY if they have the same lineage (are comparable). 
bool same_lineage(model* M, int* theta_types, int index1, int index2) {
	return is_after(M, theta_types, index1, index2, false) || is_after(M, theta_types, index2, index1, false);
}


// compute conditional expectation given the observed genotype
// Modified to only use expectations from genotypes in the same lineage.
void compute_condexp(model* M, int* theta_types, data* D, int num_unique, double* prob, double** condprob, double** exp) {
	int n = M->n;
	int m = M->m;

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < num_unique; j++) { // for all datapoints
			int compressed = binary_compression(D[j].g, n);
			int index = ABS_TO_VALID[compressed];
			double time = 0;
			for (int k = 0; k < m; k++) {  // for all valid genotypes k
				int geno = M->valids[k];
				if (geno % 2 && same_lineage(M, theta_types, index, k)) { // if parent genotype is observed
					time += exp[i][k] / prob[k] * condprob[k][j]; // time is the sum of E[T_i] * prob[k] / condprob[k][j]
				}
			}
			D[j].t[i] = time;
		}
	}
}



// Update lambda values in the outer EM loop. 
void em_mle(model* M, int* theta_types, const int n, data* D, int num_unique, double* lambdas, double S) {
	for (int i = 0; i < n; i++) {
		double sum = 0.0;
		int count = 0;
		for (int j = 0; j < num_unique; j++) {
			int index = ABS_TO_VALID[binary_compression(D[j].g, n)];
			if (!locked_out(M, theta_types, index, i)) {
				count += D[j].count;
				sum += D[j].count * D[j].t[i];
			}
		}
		lambdas[i] = (double) count / sum;
	}

	double fac = S / lambdas[n - 1];
	for (int i = 0; i < n; i++) {
		lambdas[i] *= fac;
	}
}

// Output results to console. Append "best" resuls to specified output file, if given. 
void log_em_results(char* output, double* lambdas, double epsilon, double loglik, double loglik_new, double* lambdas_opt, int iter, double* best_lambdas, double best_eps, double best, const int n) {
	printf("==== PARAMETERS ====\nLambdas: ");
	print_double_array(lambdas, n);
	printf("Epsilon = %g\n", epsilon);
	printf("loglik_new: %g\n", loglik_new);
	printf("loglik: %g\n", loglik);
	printf("norm_diff = %g\n", norm_diff(lambdas, lambdas_opt, n));
	printf("loglik_new - loglik = %g\n", loglik_new - loglik);
	printf("iter = %d\n", iter);
	printf("== best_lambdas: ==\n");
	print_double_array(best_lambdas, n);
	printf("best eps = %g\n", best_eps);
	printf("best = %g\n", best);


	if (strcmp(output, "") != 0) {
		FILE* samplelog = open_file(output, "a");
		fprintf(samplelog, "Best Lambdas:\n");
		fprint_double_array(samplelog, best_lambdas, n);
		fprintf(samplelog, "\nBest epsilon = %g\n", best_eps);
		fclose(samplelog);
	}
}


void em_cleanup(double* prob, double** exp, double** condprob, double** censprob, double* lambdas_opt, double* lambda_exit, double* marginals, double* best_lambdas, const int m, const int n) {
	free_double_array(prob);
	free_double_matrix(exp, n);
	free_double_matrix(condprob, m);
	free_double_matrix(censprob, m);
	free_double_array(lambdas_opt);
	free_double_array(lambda_exit);
	free_double_array(marginals);
	free_double_array(best_lambdas);

	free_int_array(END_STATES);
	free_int_array(LIFT);
	free_int_array(ABS_TO_VALID);
	free_int_array(VALIDS);
}

// Temporary code used for outputting conditional probabilities (Network theta values). Generates a file with all transition probabilities. 
void output_cpt(model* M, int* theta_types, double* lambdas, double* lambda_exit) {
	int m = M->m;
	int n = M->n;
	compute_lambda_exit(M, lambdas, lambda_exit);
	char* file = "cpt.txt";
	FILE* f = open_file(file, "w");

	for (int i = 0;i < m; i++) {
		int abs = M->valids[i];
		int* g = GENOTYPE[abs];
		char parent[12];
		for (int j = 0; j < n; j++) {
			if (g[j]) {
				parent[j] = '1';
			} else {
				parent[j] = '0';
			}
		}
		parent[11] = '\0';
		for (int j = 0; j < M->num_ch[i];j++) {
			char child[12];
			int diff_index = M->ch_diff[i][j];
			g[diff_index] = 1;
			for (int k = 0; k < n; k++) {
				if (g[k]) {
					child[k] = '1';
				} else {
					child[k] = '0';
				}
			}
			g[diff_index] = 0;
			child[11] = '\0';
			double fraction = lambdas[diff_index] / lambda_exit[i];
			fprintf(f, "%s,%s,%g\n", parent, child, fraction);
		}
	}
	fclose(f);
}

// runs lambda EM here, calls epsilon_em to get max likleihood epsilon value for given sets of lambdas. 
double expectation_maximization(model* M, int* theta_types, data* D, int num_unique, double* lambdas, char* output) {
	int m = M->m;
	int n = M->n;
	double S = lambdas[n - 1];

	double* prob = get_double_array(m);
	double** exp = get_double_matrix(n, m);
	double** condprob = get_double_matrix(m, num_unique); // probability of each valid genotype, conditioned on each datapoint. 
	double** censprob = get_double_matrix(m, m);
	double* lambdas_opt = get_double_array(n);
	double* lambda_exit = get_double_array(m);
	double* marginals = get_double_array(num_unique); // one per column in condprob
	double* best_lambdas = get_double_array(n);

	double loglik = 0, loglik_new = 0, best = -100000000.0;

	double epsilon = 0.25, best_eps = 0.25; // epsilon set at 0.25 to give a "warm" start to the epsilon EM
	int iter = 0;
	bool invalid_values = false;
	while((norm_diff(lambdas, lambdas_opt, n) > EM_ACCURACY || (iter < 2)) && (iter < EM_MAX_ITER) && (loglik_new - loglik > 1e-4 || (iter < EM_MAX_ITER / 100))) {
		if (iter > 100 && loglik_new < loglik) { // shouldn't happen, but throw an error and break out of loop if it does.
			fprintf(stderr, "Warning in expectation_maximization: likelihood decreasing.\n");
			break;
		}
		if (!invalid_values) {
			for (int i = 0; i < n; i++) {
				if (isnan(lambdas[i]) || isinf(lambdas[i])) invalid_values = true;
			}
			if (invalid_values) {
				printf("NaN or inf values detected.\n");
				if (strcmp(output, "") != 0) {
					FILE* f = open_file(output, "a");
					fprintf(f, "=== Warning ===\n");
					fprintf(f, "NaN or inf values detected.\n");
					fclose(f);				
				}
			}
		}


		loglik = loglik_new;
		if (iter != 0 && loglik > best) {
			memcpy(best_lambdas, lambdas, sizeof(double) * n);
			best_eps = epsilon;
			best = loglik;
		}
		for (int i = 0; i < n; i++) {
			lambdas_opt[i] = lambdas[i];
		}

		compute_lambda_exit(M, lambdas, lambda_exit);
		compute_all_prob(M, lambdas, lambda_exit, prob);

		// == E-step == 
		// output: new conditional probabilities for each valid genotype given the EM-ed epsilon
		loglik_new = EM_epsilon(M, theta_types, D, num_unique, prob, condprob, &epsilon, marginals);

		// output: censprob, used to compute expected times
		compute_censored_prob(M, theta_types, lambdas, lambda_exit, prob, censprob); 
		compute_all_exp(M, theta_types, lambdas, lambda_exit, prob, exp, censprob);
		compute_condexp(M, theta_types, D, num_unique, prob, condprob, exp);

		// == M-step == update lambda values
		em_mle(M, theta_types, n, D, num_unique, lambdas, S);

		iter++;
		if (iter % 100 == 0) {
			printf("outer loop iteration %d: loglik_new = %g, epsilon = %g\n", iter, loglik_new, epsilon);
		}
	}

	log_em_results(output, lambdas, epsilon, loglik, loglik_new, lambdas_opt, iter, best_lambdas, best_eps, best, n);
	memmove(lambdas, best_lambdas, sizeof(double) * n);

	// output_cpt(M, theta_types, lambdas, lambda_exit);

	em_cleanup(prob, exp, condprob, censprob, lambdas_opt, lambda_exit, marginals, best_lambdas, m, n);
	return best_eps;
}

#endif
