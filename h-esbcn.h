#ifndef BN
#define BN

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdbool.h>
#include <getopt.h>
#include <assert.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_matrix.h>

#include <omp.h>
#include <math.h>

#include "matrix.h"
#include "model.h"
#include "data.h"
#include "vector.h"
#include "a-queue.h"

bool verbose = false;
bool debug = false;

int** GENOTYPE; // GENOTYPE[i] stores the binary expansion of i.
int NUM_SAMPLES_READ; // total number of samples read in from dataset.

// set initial parameter values (seeimingly arbitrary, but makes it easier to tweak them since they're all here)
void initialize_MH_params(unsigned int* seed, int* number_samples, int* record_ith, int* burn_in, char* reg_scheme) {
	*seed = time(NULL); // random seed using time
	*number_samples = 100000;
	*record_ith = 1;
	*burn_in = 0;
	strcpy(reg_scheme, "bic");
}

/*
 * Check file inputs to bn. If no .poset file supplied (f flag), prints error message and exits.
 * Flags: 
 *  
 * 
 * 
 */
void check_file_inputs(int argc, char** argv, char* filestem, char* output, int* number_samples, unsigned int* seed, char* reg_scheme) {
	int c = 0;
	bool input_flag = false;
	bool output_flag = false;
	bool help_flag = false;
	static struct option long_options[] =
        {
          {"dataset",     		required_argument, 0, 'd'},
          {"output",  				required_argument, 0, 'o'},
          {"number_samples",	required_argument, 0, 'n'},
          {"seed", 						required_argument, 0, 's'},
          {"reg", 						required_argument, 0, 'r'},
          {"help", 						no_argument, 			 0, 'h'},
          {NULL, 0, NULL,0,}
        };
	while ((c = getopt_long(argc, argv, "d:o:n:s:r:hvt:e:", long_options, NULL)) != -1 ) {
		switch(c) {
			case 'd': // get input filename, satisfy boolean check for provided input. 
				strcpy(filestem, optarg);
				input_flag = true;
				break;
			case 'o':
				strcpy(output, optarg);
				output_flag = true;
				break;
			case 'n':
				*number_samples = atoi(optarg);
				break;
			case 's':
				*seed = atoi(optarg);
				break;			
			case 'r':
				strcpy(reg_scheme, optarg);
				break;
			case 'h':
				help_flag = true;
				break;
			case 'v':
				verbose = true;
				break;
		}
	}

	if (help_flag) {
		printf("Usage: ./cbn (-(option) (argument) )*\n");
		printf("where dataset input option (-d or --dataset) is required.\n");
		printf("Options include:\n");
		printf("-d | --dataset \tPath to input file with sample data. Please include file extension as well. ex: test1.txt\n");
		printf("-o | --output \tPath to output file to write results to. If not specified, then results not written to file.\n");
		printf("-n | --number_samples \tSet the number of iterations for Metropolis-Hastings to be run. Default value set to 500,000.\n");
		printf("-s | --seed \tSet the seed for random number generators in the program. Time used as default seed.\n");
		printf("-r | --reg \tSet the likelihood regularization scheme for the MCMC algorithm. Valid choices are bic, aic, and loglik.\n");
		printf("-h | --help \tPrint help information for program.\n");
		exit(1);
	}

	if (!input_flag) {
		fprintf(stderr, "Error: No input file specified. Please provide the filename (without filetype) of a file containing formatted sample data using the \"-d\" or \"--dataset\" options.\n");
		exit(1);
	}
	if (!output_flag) {
		strcpy(output, "");
	}
	if (strcmp(reg_scheme, "bic") &&  strcmp(reg_scheme, "aic") && strcmp(reg_scheme, "loglik")) {
		fprintf(stderr, "Error: Invalid regularization scheme. Please choose from bic (default), aic, or loglik.\n");
		exit(1);
	}	
}

/* 
 * Standard open_file helper. Throws error if error while opening file. 
 */
FILE* open_file(char* filename, char* option) {
	FILE* input = fopen(filename, option);
	if (input == NULL) {
		fprintf(stderr, "Error: Could not read %s\n", filename);
		exit(1);
	}
	return input;
}

/* 
 * Base 2 power function. Computes 2^n, for n >= 0
 */
static inline int pow2(int n) {
	assert(n >= 0);
	return 1 << n;
}


// Basic power function implementation. Only accepts integral power values. 
double power(double base, int power) {
	double result = 1.0;
	for (int i = 0; i < power; i++) {
		result *= base;
	}
	return result;
}


/* 
 * Generates a binary compression of pattern data, which is
 * used as a means of quickly determining whether two samples have 
 * the same pattern. Presumably, this is faster (as opposed to
 * comparisons without such a compression scheme), since we compress 
 * all samples once, and then compare the "hash" values. 
 * Parallels index_of in TS code. 
 */
int binary_compression(int* pat_i, int num_events) {
	int index = 0;
	for (int i = 0; i < num_events; i++) {
		if (pat_i[i] ==1)
			index += pow2(num_events - 1 - i);
	}
	return index;
}

// Used to generate binary expansions of genotypes, for easier access later on.
void precompute_binary(const int n) {
	int m = pow2(n);
	GENOTYPE = get_int_matrix(m, n);
	for (int i = 0; i < m; i++) {
		int count = 0;
		for (int j = 0; j < n; j++) {
			if ( ((1 << j) & i) > 0) {
				GENOTYPE[i][n - 1 - j] = 1;
				count++;
			} else {
				GENOTYPE[i][n - 1 - j] = 0;
			}
		}
	}
}

int count_edges(int** poset, int n) {
	int count = 0;
	for (int i=0; i< n; i++) {
		for (int j = 0; j< n; j++) {
			if (poset[i][j])
				count++;
		}
	}
	return count;
}

// TODO: reconsider the use of int* ptr vector here. We've changed the idea, and we're no longer looking
// at the probability of combined events. Instead, using the current poset to determine how children
// are valid or not in find_children_extended. 

/* 
 * Helper function for process_patterns function. Reads pattern data
 * from .txt file. Only 0/1 expected, but places negative values where
 * any deviations occur. Exits if any error is encountered.
 */
int** read_patterns(FILE *input, model* M, char* filename, int* num_samples, int* num_events) {
	size_t temp = 0;
	int n;
	char* line = NULL; 
	int len = getline(&line, &temp, input);
	if (len != -1) {
		n = len / 2;
		free(line);
		line = NULL;
		rewind(input);
	} else {
		fprintf(stderr, "Error reading sample data.\n");
		exit(1);
	}

	M->n = n;
	*num_events = n;

	ptr_vector* v = new_ptr_vector(); // vector takes ints, hopefully won't get mutilated?
	while (getline(&line, &temp, input) != -1) {
		int line_index = 0;
		int* genotype = get_int_array(n);
		for (int i = 0; i< n; i++) { // fill in each event with genotype.
			if (line[line_index] == '0') genotype[i] = 0;
			if (line[line_index] == '1') genotype[i] = 1;
			line_index += 2;
		}
		push_back_ptr(v, genotype);
		free(line);
		line = NULL;
		temp = 0;
	}
	free(line);
	*num_samples = v->size;

	if (debug) {
		printf("Finished reading patterns.\nPatterns matrix: \n");
		print_int_matrix(v->elems, v->size, n);
	}

	return free_ptr_vector(v, false);
}


// counts uniq patterns among samples.
// It seems that pat_idx counts the number of unique samples
// up to that last unique sample.
// ----- variables -----
// index stores the "hashed" patterns
// count stores the number of samples with the same pattern 
// Removed pat_idx necessity, since the same information is stored in index. 

// Note: The binary compression places a limit on the number of loci that can be observed (32 for int, 64 for long int). 
// For a completely scalable implementation, should re-work this scheme, use a hashset for speed. 
void count_uniq_patterns(int* index, int* count, int num_samples, int num_events, int* N_u, int** pat) {

	//binary compression that "hashes" genotypes 
	for (int i = 0; i< num_samples; i++)
		index[i] = binary_compression(pat[i], num_events);

	// count unique patterns. Runs in worst case O(n^2) time. can be optimized to O(nlgn) using sort -> iterate, or linear time using a hashtable/set.
	bool skips; // true if another sample has same pattern ("index" value)
	int num_unique = 0; // counts number of unique patterns
	for (int i = 0; i < num_samples; i++) {
		skips = false;
		for (int j = 0; j < i; j++) {
			if (index[i] == index[j]) { //samples have the same pattern
				count[j]++; //increment the count at first index
				skips = true;
				break;
			}
		}
		if (!skips) { // unique pattern --> NOT previously observed
			count[i]++;
			num_unique++;	
		}
	}
	*N_u = num_unique;

	if (verbose) {
		printf("Number of samples processed = %d\n", num_samples);
		printf("Number of gene loci considered = %d\n", num_events);
		printf("Number of unique genotypes observed = %d\n", num_unique);
	}
	if (debug) {
		printf("num_samples = %d\n", num_samples);
		printf("count array: \n");
		print_int_array(count, num_samples);
		printf("index array: \n");
		print_int_array(index, num_samples);
	}
}

data* construct_data_set(int num_unique, int* count, int* index, int num_samples, int num_events, int** pat) {
	data* D = calloc(num_unique, sizeof(data));
	if (D == NULL) {
		fprintf(stderr, "Error: calloc failure - out of memory.\n");
		exit(1);
	}
	int d_index = 0;
	for (int i = 0; i < num_samples; i++) {
		if (count[i] > 0) {//if was not a repeat
			D[d_index].count = count[i];
			D[d_index].g = get_int_array(num_events);
			for (int j = 0; j < num_events; j++) 
				D[d_index].g[j] = pat[i][j];
			d_index++;
		}
	}

	if (debug) {
		print_data_array(D, num_unique, num_events);
	}
	return D;
}

// N_u = number of unique genotypes in dataset
// n  = number of events
int process_patterns(char* filename, model* M, int* N_u, data** D) {

	FILE* input = open_file(filename, "r");
	int num_samples, num_events;
	int** pat = read_patterns(input, M, filename, &num_samples, &num_events);
	fclose(input);

	//Count unique patterns, construct data set.
	int* index = get_int_array(num_samples);
	int* count = get_int_array(num_samples);
	count_uniq_patterns(index, count, num_samples, num_events, N_u, pat);
	*D = construct_data_set(*N_u, count, index, num_samples, num_events, pat);

	NUM_SAMPLES_READ = num_samples;
	free_int_matrix(pat, num_samples);
	free_int_array(index);
	free_int_array(count);

	return 0;
}

// Suppes condition
// Confirm that the proposed modification (move) is supported by event probabilities in the data set, taking into account
// current theta types. 
bool cumulative_probs_maintained(model* M, data* D, int* theta_types, int num_unique) {
	int n = M->n;
	for (int i = 0; i< n; i++) { // loop through each event
		int num_parents = M->num_pa[i];
		int* parents = M->pa[i];
		int type = theta_types[i];
		int count_parents = 0;
		int count_child = 0;
		int count_parents_child = 0;
		for (int j = 0; j< num_unique; j++) { // loop through all observed genotypes
			int* g = D[j].g;
			int found = 0;
			
			// get child present bool value
			bool child_present = false;
			if (g[i]) child_present = true;

			// get parent set satisfied bool value
			bool parents_present = true;
			switch (type) { // 0 or 1 parents, conjunctive
				case -1 ... 0:
					for (int k = 0; k < num_parents; k++) {
						if (!g[parents[k]]) parents_present = false;
					}
					break;
				case 1: // or
					for (int k = 0; k < num_parents; k++) {
						if (g[parents[k]]) found++;
					}
					if (!found) parents_present = false;
					break;
				case 2: // xor
					for (int k = 0; k < num_parents; k++) {
						if (g[parents[k]]) found++;
					}
					if (found != 1) parents_present = false;
					break;
			}

			if (child_present && parents_present) {
				count_parents_child += D[j].count;
			}

			if (child_present) {
				count_child += D[j].count;
			}

			if (parents_present) {
				count_parents += D[j].count;
			}
		} 
		// check that the checks and additions are correct. 
		if (num_parents > 0) {
			if (count_parents <= count_child || count_parents_child * NUM_SAMPLES_READ <= count_parents * count_child) {
				return false;
			}
		} 
	}
	return true;
}

// Initializes the parent fields in the model so that they can be used to find children more efficiently. 
void set_theta_types(model* M, int* theta_types, int* theta_types_p) {
	int num_events = M->n;
	int** poset = M->P;
	free_parents(M);
	M->num_pa = get_int_array(num_events); // create an array that stores the sizes of parents. 
	M->pa = (int**) malloc(num_events * sizeof(int*)); // create an array to hold arrays of integers

	// for all events, create a vector of all other events that relate to them, according to the current poset. 
	// Then, determine theta type if an event has > 1 parent events. 
	for (int i = 0; i< num_events; i++) {
		vector* v = new_vector();
		for (int j = 0; j< num_events; j++) {
			if (poset[j][i]) {
				push_back(v, j);
			}
		}
		int num_parents = v->size;
		M->num_pa[i] = num_parents;
		M->pa[i] = free_vector(v, false);

		if (num_parents <= 1) { // if 1 parent at most, no need for special care. 
			theta_types_p[i] = -1;
		} else {
			if (theta_types[i] < 0) { // event goes from 0/1 parents to 2+ parents
				theta_types_p[i] = rand() % 3; // randomly pick between 0, 1, and 2. 
			} else { // carry over formula from previous iteration. 
				theta_types_p[i] = theta_types[i];  
			}
		}
	}
}


// imposes a total ordering on the events. Use Fisher-Yates shuffling algorithm to randomize ordering
void new_total_ordering(model* M, int* theta_types, data* D, int num_unique) {
	int num_events = M->n;
	// initialize new array to use.
	int rand_arr[num_events];
	for (int i = 0; i < num_events; i++) {
		rand_arr[i] = i;
	}

	// establish total ordering
	for (int i = num_events - 1; i >= 0; i--) {
		int randy = rand() % (i + 1);
		int temp = rand_arr[i];
		rand_arr[i] = rand_arr[randy];
		rand_arr[randy] = temp;
	}

	int** initial_poset = get_int_matrix(num_events, num_events);
	// mark in poset. 
	for (int i = 0; i < num_events - 1; i++) {
		initial_poset[rand_arr[i]][rand_arr[i+1]] = 1;
	}

	// transitive closure & theta types
	int** T = get_int_matrix(num_events, num_events);
	transitive_closure(initial_poset, T, num_events);

	for (int i = 0; i < num_events; i++) {
		int num_parents = 0;
		for (int j = 0; j < num_events; j++) {
			if (T[j][i]) num_parents++;
		}
		if (num_parents == 0) continue; // don't bother checking if there are no parents, since we're only going to remove excess edges.
		int reduce = num_parents - (rand() % (num_parents + 1));
		for (int j = 0; j < reduce; j++) { // delete a random relation reduce number of times.
			int ith_event = rand() % num_parents;
			int counted = 0;
			for (int k = 0; k < num_events; k++) {
				if (counted == ith_event && T[k][i]) {
					T[k][i] = 0; // delete event
					num_parents--; // decrement number of parents.
					break;
				}
				if (T[k][i]) counted++;
			}
		}
	}

	free_int_matrix(initial_poset, num_events);
	free_int_matrix(M->P, num_events);
	M->P = T;
	
	int* theta_types_p = get_int_array(num_events);
	for (int i = 0; i < num_events; i++) {
		theta_types_p[i] = -1;
	}
	set_theta_types(M, theta_types_p, theta_types);

	model M_p;
	M_p.n = num_events;
	M_p.pa = M->pa;
	M_p.P = T;
	M_p.num_pa = get_int_array(num_events); // initializes each to 0.
	// do suppes checking for each events' parent set. 
	for (int i = 0; i < num_events; i++) {
		// determine suppes prob condition for this event and its parent set. 
		// copy the parent set of this event into a new poset. 
		theta_types_p[i] = theta_types[i];
		M_p.num_pa[i] = M->num_pa[i];
		if (!cumulative_probs_maintained(&M_p, D, theta_types, num_unique)) {// reject entire parent set.
			for (int j = 0; j < num_events; j++) { 
				M->P[j][i] = 0;
			}
			M->num_pa[i] = 0; 
			theta_types[i] = -1;
		}
		theta_types_p[i] = -1; // push to conjunctive so that it doesn't go to 'or' and 'xor' cases.
		M_p.num_pa[i] = 0; // set back to 0.
	}
	free_int_array(M_p.num_pa);
	free_int_array(theta_types_p);
}


// returns the total number of events associated with a given genotype. 
// ex. For n= 4, 0101 returns 2.
int total_events(const int genotype, const int n) {
	int* binary = GENOTYPE[genotype];
	int num_events = 0;
	for (int i=0; i< n; i++) 
		num_events += binary[i];
	return num_events;
}

// Given two integers, returns the number of bitwise differences between the two. 
// If hamming distance > 1, returns index of last (largest) bit. 
int hamming_distance(unsigned int x, unsigned int y, const int num_events, int* diff_index) {
	unsigned int dist = 0;
	unsigned int diff = x ^ y;
	for (int i= 0; i < num_events; i++) {
		if ( (1 << i) & diff ) {
			dist++;
			if (diff_index != NULL)	*diff_index = num_events - 1 - i;
		}
	}
	return dist;
}


// Checks that the parent conditions of the ith event are satisfied in provided genotype g. 
// Returns true if parent conditions (according to poset and theta_types) are fulfilled, false otherwise. 
bool check_parent_set(model* M, int* theta_types, int* g, int i) {
	int num_parents = M->num_pa[i];
	int type = theta_types[i]; // relationship between parents.
	bool compatible = false;
	int found = 0;

	switch (type) { // will only go into one case for each, so this may be faster. 
		case -1 ... 0: // 0/1 parents, or and case if > 1 parents. loop through all parents, confirm that all are present.
			for (int j = 0; j< num_parents; j++) {
				if (!g[M->pa[i][j]]) return false;
			}
			break;
		case 1: // or. 
			for (int j = 0; j< num_parents; j++) {
				if (g[M->pa[i][j]]) {
					compatible = true;
					break;
				}					
			}
			if (!compatible) return false;
			break;
		case 2: // xor. loop through all parents, confirm that only one present at a time.
			for (int j = 0; j< num_parents; j++) {
				if (g[M->pa[i][j]]) { // another event observed that isn't the current
					found++;
				} 
			}
			if (found != 1) return false;
			break;
	}
	return true;
}


// Theta types gives the relationship between "parent events", if an event has more than one direct parent.
// If check_parent_set returns true for all events, then returns true. 
bool compatible_with_poset(int genotype, model* M, int* theta_types) {
	int num_events = M->n;
	int* g = GENOTYPE[genotype];
	for (int i = 0; i< num_events; i++) { // for each event
		if (g[i]) { // if event is observed
			if (!check_parent_set(M, theta_types, g, i)) return false; 
		}
	}
	return true;
}

#endif
