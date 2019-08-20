#ifndef MCMC
#define MCMC

#include "h-esbcn.h"
#include "a-queue.h"
#include "vector.h"
#include "likelihood_computation.c"


#define min(a, b) ( a < b ? a: b) // switch to macro to generalize types


// Percentage weightings for each structural move type. Remaining 0.15 is reserved for relocating theta types. 
double fraction_exchange = 0.07;
double fraction_reincarnation = 0.1;
double fraction_birth = 0.4;
double fraction_death = 0.2;
double fraction_tr_birth = 0.04;
double fraction_tr_death = 0.04;


// upper threshold values for each move type. Remaining 0.15 is reserved for relocating theta types. 
double exchange = 0.07;
double reincarnate = 0.17;
double birth = 0.57;
double death = 0.77;
double tr_birth = 0.81;
double tr_death = 0.85;

// gets two distinct random numbers. 
void get_two_rand(int* x, int* y, int n) {
	*x = rand() % n;
	*y = rand() % n;
	while(x == y) {// ensure that the two numbers are not equal
		*y = rand() % n;
	}
}

// Swaps the relationships of two events. Conceptually, swap positions of 2 nodes in network. 
void propose_event_exchange_move(model *M, double* tp, int* ti, int* tj, int x, int y) {
	int n = M->n;
	int* in = get_int_array(n);
	int* out = get_int_array(n);
	int** P = M->P;

	int xy = P[x][y];
	int yx = P[y][x];
	for (int i = 0; i < n; i++) {
		out[i] = P[x][i]; // out holds everything in row x
		in[i] = P[i][x]; // in holds everything in column x
	}
	for (int i = 0; i< n; i++) { // cannot merge with previous loop since this one changes values
		if (i != x) { // for all except where x= i
			P[x][i] = P[y][i]; // set everything in row x equal to row y
			P[i][x] = P[i][y]; // set everything in column x equal to column y.
		}
	}
	for (int i = 0; i< n; i++) {
		if (y != i) {
			P[y][i] = out[i];
			P[i][y] = in[i];
		}
	}
	P[x][y] = yx;
	P[y][x] = xy;

	*ti = x;
	*tj = y;

	*tp = 2.0 / (double)(n * (n-1));

	free_int_array(in);
	free_int_array(out);
}

// Count number of parents a certain event has. 
int count_event_parents(int** poset, int num_events, int j) {
	int count = 0;
	for (int i = 0; i< num_events; i++) {
		if (poset[i][j]) {
			count++;
		}
	}
	return count;
}


// Propose edge deletion event. Uses random number randy to as the edge to delete. 
bool propose_delete_cover_relation(model* M, double* tp, int randy) {
	int c = 0;
	int num_events = M->n;
	int** P = M->P;
	int valid_c = count_edges(P, num_events);
	if (valid_c <= 0) { // no relations to delete
		return false;
	}

	*tp = 1.0 / (double) valid_c;
	// int new_index = rand() % valid_c; // select a random number in [0, valid_c]
	int new_index = randy % valid_c;
	for (int i = 0; i< num_events; i++) {
		for (int j = 0; j< num_events; j++) {
			if (P[i][j]) {
				if (new_index == c) { //when hit the random number, delete that existing relation. 
					P[i][j] = 0;
					M->P = P;
					return true;
				}
				c++;
			}
		}
	}
	return true;
}

// this is probably correct, but is there a better way to do it? 
// it seems to me that it might be looping through needlessly at the end- is there a DP solution we can employ?
// Return codes: 0 = OK, 1 = redundant edge, 2 = cyclic edge. 
// Go through relations matrix and ensure that only cover relations (direct edges) are present. 
int reduce_to_cover_relations(int** P, int n, int** C) {
	int stat = 0;
	queue* q = new_queue();

	for(int i = 0; i< n; i++) { // for all nodes
		for(int j = 0; j< n; j++) { // queue i's children
			if (P[i][j])
				enqueue(q, j);
		}

		int* visited = get_int_array(n); // keep track of which children have been visited
		while(!is_empty(q)) { 
			int child = dequeue(q);
			for(int k = 0; k< n; k++) { 
				if (P[child][k] && !visited[k]) { // if node is i's grandchild and not visited
					visited[k] = 1;
					enqueue(q, k);

					// remove non-cover relations. 
					// i is not directly related, since we're looking at i's grandchildren. e.g. not a cover relation
					if (P[i][k]) {
						P[i][k] = 0;
						C[i][k] = 1;
						stat = 1; // return code for a redundant edge.
					}

					//if cyclic. Recall that all descendents are queued, with i and j remaining constant.
					if (P[k][i]) { 
						free_int_array(visited);
						free_queue(q);
						return 2; // return code for a cyclic edge.
					}
				}
			}
		}

		free_int_array(visited);
	}
	free_queue(q);
	return stat;
}

// Evaluate the implications of the proposed change (reduce to cover relations), and reverse change if needed. 
// Refactored for readability. 
void check_new_relations(model* M, int** T, int** C, int** V, int i, int j, int n) {
	int stat = reduce_to_cover_relations(T, n, C);
	if (stat == 1 && !C[i][j]) { // if there are redundant edges that are not fresh changes, find and remove.
		for (int k = 0; k< n; k++) {
			for (int l = 0; l < n; l++) {
				if (C[k][l] && M->P[k][l]) 
					T[i][j] = 0;
			}
		}
		V[i][j] = T[i][j];
	} else if (stat== 2) { // made a cycle, delete change. 
		T[i][j] = 0;
	} else if (stat == 0) { // valid change. Mark change in V. 
		V[i][j] = 1;
	}
}

// Fill matrix V with valid new cover relations. 
void valid_new_cover_relation_matrix(model* M, int** V) {
	int n = M->n;
	int** T = get_int_matrix(n, n); // the transitive closure of relations
	int** C = get_int_matrix(n, n); // matrix marking any changes that have been made.

	transitive_closure(M->P, T, n);

	for (int i = 0; i< n; i++) {
		for(int j = 0; j < n; j++) {
			if (i!=j && !T[i][j]) { // if two distinct events are not related
				T[i][j] = 1; // propose new relation

				// check effects of change. Refactored for readability.
				check_new_relations(M, T, C, V, i, j, n);
			}

			// reinitialize T and C for next iteration. 
			for(int l = 0; l<n; l++) {
				for(int k = 0; k< n; k++) {
					T[l][k] = 0;
					C[l][k] = 0;
				}
			}
			transitive_closure(M->P, T, n);
		}
	}
	free_int_matrix(T, n);
	free_int_matrix(C, n);
}

// Propose a new edge to add to poset. 
bool propose_new_cover_relation(model* M, double* tp, data* D, int num_unique, int randy) {
	int n = M->n;
	int** V = get_int_matrix(n, n); // Matrix containing valid relations to try. 
	valid_new_cover_relation_matrix(M, V); // Populate V

	// Count the number of valid relations, pick one at random. 
	int valid_c = count_edges(V, n);

	int c = 0;
	if (valid_c) {
		// int new_index = rand() % valid_c;
		int new_index = randy % valid_c;
		for (int i = 0; i< n; i++) {
			for (int j = 0; j< n; j++) {
				if (V[i][j]) {
					c++;
					if (new_index < c) {
						M->P[i][j] = 1;
						*tp = 1.0 / valid_c;
						free_int_matrix(V, n);
						return true;
					}
				}
			}
		}
	} 
	free_int_matrix(V, n);
	return false;
}

// Propose a death move, then a birth move. Fails if no death move possible. If no birth move possible, death only. 
bool propose_reincarnation_move(model* M, data* D, int num_unique) {
	double tp;
	int x, y;
	get_two_rand(&x, &y, M->n);

	if (!propose_delete_cover_relation(M, &tp, x))  {
		return false;
	}
	propose_new_cover_relation(M, &tp, D, num_unique, y);

	return true;
}

// Get TP value for propose_cover_relation
double get_new_cover_move_tp(model *M, data* D, int num_unique) {
	int n = M->n;
	int** V = get_int_matrix(n, n);
	valid_new_cover_relation_matrix(M, V);
	int valid_c = count_edges(V, n);
	free_int_matrix(V, n);
	return 1.0 / valid_c;
}

// Check that a transitive move is valid- if yes, then mark the index of move with 1 in V. 
void transitive_check_relations(model* M, int** V, int** T2, int** C, int n, int i, int j) {
	int stat = reduce_to_cover_relations(T2, n, C);
	if (stat == 1 && !C[i][j]) {
		if (abs_matrix_diff(T2, M->P, n) > 2) {
			V[i][j] = 1;
		}
	}
}

// Populates matrix V with all valid moves after performing a transitive closure event on current poset. 
void valid_new_transitive_matrix(model* M, int** V) {
	int n = M->n; // num events. e.g. matrix size
	int** T = get_int_matrix(n, n);
	int** T2 = get_int_matrix(n, n);
	int** C = get_int_matrix(n, n);

	transitive_closure(M->P, T, n);
	int orig_n_edges = count_edges(T, n);

	for(int i = 0; i< n; i++) {
		for(int j = 0; j< n; j++) {
			if (!T[i][j] && i != j) { // pick a relation that isn't present and isn't a self-relation
				T[i][j] = 1; // try this relation
				transitive_closure(T, T2, n);
				int new_n_edges = count_edges(T2, n);
				if (new_n_edges == orig_n_edges + 1) {
					transitive_check_relations(M, V, T2, C, n, i, j);
				}
				T[i][j] = 0; // restore original value. 

				// reinitialize matrices. 
				for(int k = 0; k< n; k++) {
					for (int l = 0; l< n; l++) {
						T2[l][k] = 0;
						C[l][k] = 0;
					}
				}
			}
		}
	}
	free_int_matrix(T, n);
	free_int_matrix(T2, n);
	free_int_matrix(C, n);
}

// Propose a new edge after taking a transitive closure. 
bool propose_new_transitive_closure(model* M, double* tp, int randy) {
	int n = M->n;
	int** V = get_int_matrix(n, n);
	int** T = get_int_matrix(n, n);
	int** C = get_int_matrix(n, n);

	valid_new_transitive_matrix(M, V);
	int valid_c = count_edges(V, n);

	if (!valid_c) { // temporary solution for if valid_c comes back with 0. 
		free_int_matrix(V, n);
		free_int_matrix(T, n);
		free_int_matrix(C, n);
		return false;
	}

	*tp = 1.0 / (double) valid_c;
	transitive_closure(M->P, T, n);
	int c = 0;
	// int new_index = rand() % valid_c;
	int new_index = randy % valid_c;
	for(int i = 0; i< n; i++) {
		for(int j = 0; j< n; j++) {
			if(V[i][j]) {
				c++;
				if (new_index < c) {
					T[i][j] = 1;
					reduce_to_cover_relations(T, n, C);
					free_int_matrix(M->P, n);
					M->P = T;
					free_int_matrix(C, n);
					free_int_matrix(V, n);
					return true;
				}
			}
		}
	}
	return false;
}

// Populates matrix V with vali deletion events for the transitive closure of the current poset. 
void valid_delete_transitive_matrix(model* M, int** V) {
	int n = M->n; // num events. e.g. matrix size
	int** T = get_int_matrix(n, n);
	int** T2 = get_int_matrix(n, n);
	int** C = get_int_matrix(n, n);

	transitive_closure(M->P, T, n);

	for(int i = 0; i< n; i++) {
		for(int j = 0; j< n; j++) {
			if (T[i][j] && i != j) { // pick a relation that is present and isn't a self-relation
				T[i][j] = 0; // try this relation
				transitive_closure(T, T2, n);
				if (!T2[i][j]) {
					transitive_check_relations(M, V, T2, C, n, i, j);
				}
				T[i][j] = 1; // restore original value. 

				// reinitialize matrices. 
				for(int k = 0; k< n; k++) {
					for (int l = 0; l< n; l++) {
						T2[l][k] = 0;
						C[l][k] = 0;
					}
				}
			}
		}
	}
	free_int_matrix(T, n);
	free_int_matrix(T2, n);
	free_int_matrix(C, n);
}

// Get TP for transitive_delete event. 
double get_transitive_closure_delete_tp(model* M) {
	int n = M->n; // num events, matrix size.
	int** V = get_int_matrix(n, n);
	int** T = get_int_matrix(n, n);
	int** C = get_int_matrix(n, n);

	valid_delete_transitive_matrix(M, V);
	int valid_c = count_edges(V, n);

	free_int_matrix(V, n);
	free_int_matrix(T, n);
	free_int_matrix(C, n);
	if (valid_c > 0) {
		return 1.0 / (double) valid_c;
	}
	return 0;
}

// Propose an edge deletion event after taking transitive closure of the current poset. 
bool propose_delete_transitive_closure(model* M, double* tp, int randy) {
	int n = M->n;
	int** V = get_int_matrix(n, n);
	int** T = get_int_matrix(n, n);
	int** C = get_int_matrix(n, n);

	valid_delete_transitive_matrix(M, V);
	int valid_c = count_edges(V, n);

	if (valid_c <= 0) {
		free_int_matrix(V, n);
		free_int_matrix(T, n);
		free_int_matrix(C, n);
		return false;
	}

	int c = 0;
	// int new_index = rand() % valid_c;
	int new_index = randy % valid_c;
	for(int i = 0; i< n; i++) {
		for(int j = 0; j< n; j++) {
			if( V[i][j] ) {
				c++;
				if (new_index < c) {
					T[i][j] = 1;
					reduce_to_cover_relations(T, n, C);
					*tp = 1.0 / (double) valid_c;
					free_int_matrix(M->P, n);
					M->P = T;
					free_int_matrix(C, n);
					free_int_matrix(V, n);
					return true;
				}
			}
		}
	}
	return true;
}

// Get TP for transitive_birth event.
double get_transitive_closure_new_tp(model* M) {
	int n = M->n;
	int** V = get_int_matrix(n, n);

	valid_new_transitive_matrix(M, V);
	int valid_c = count_edges(V, n);
	free_int_matrix(V, n);
	if (valid_c > 0) {
		return 1.0 / (double) valid_c;
	}
	return 0;
}

// Updates current model to reflect the proposed move and prints logging information to console. 
void make_move(long double* posterior_o, long double posterior_p, model* M, model* M_p, int* theta_types, int* theta_types_p, double ms, int* num_edges) {
	*posterior_o = posterior_p;
	int n = M->n;
	memcpy(theta_types, theta_types_p, sizeof(int) * n);
	free_model(*M, false);
	model_copy(M_p, M);
	if ( ms < exchange) {
		// printf("Success: EXCHANGE move. ");
	}
	else if ( ms < reincarnate) {
		// printf("Success: REINCARNATION move. ");
	}
	else if ( ms < birth) {
		*num_edges += 1;
		// printf("Success: BIRTH move. ");
	}
	else if ( ms < death ) {
		*num_edges -= 1;
		// printf("Success: DEATH move. ");
	} else {
		*num_edges = count_edges(M->P, n);
		if ( ms < tr_birth ) {
			// printf("Success: transitive BIRTH move. ");
		} else if (ms < tr_death) {
			// printf("Success: transitive DEATH move. ");
		} else if (ms > 1) {
			// printf("Aggressive move taken: New poset. ");
		}
	} 
	// printf("posterior %.10Lg. number edges %d.\n", posterior_p, *num_edges);
}

// Relocates theta type for a random event to one of the other two possible event types. 
static inline void relocate_theta_type(int* theta_types_p, int randy) {
	int curr = theta_types_p[randy];
	if (curr >= 0) {
		curr += rand() % 2;
		theta_types_p[randy] = curr % 3;
	}
}

// Log poset from best (most likely) model to output file. Also prints out some aditional information. 
void log_results(int** P, int* best_theta_types, int accepted, int actual_iterations, FILE* samplelog, int suboptimal, int largest, int n) {
	
	if (samplelog != NULL) {
		fprintf(samplelog, "best poset: \n");
		fprint_int_matrix(samplelog, P, n, n);
		fprintf(samplelog, "\n accepted: %d, ratio %g\n", accepted, (float) accepted / actual_iterations);
		fprintf(samplelog, "best theta types = \n");
		fprint_int_array(samplelog, best_theta_types, n);		
	}

	printf("best poset:\n");
	print_int_matrix(P, n, n);
	printf("\n accepted: %d, ratio %g\n", accepted, (float)accepted/actual_iterations);
	printf("best theta types = \n");
	print_int_array(best_theta_types, n);
	printf("largest number of consecutive failures = %d\n", largest);
	printf("suboptimal: %d, percentage : %g\n", suboptimal, (double) suboptimal / accepted);

}


	/**
	 * suffixes of variables in this method
	 * *_p ... proposal
	 * *_o ... old value
	 * *_r ... reverse
	 */
void metropolis_hastings(model* M, data* D, int num_unique, int burn_in, int number_samples, int record_ith, unsigned int rand_seed, char* output, int* theta_types, char* scheme) {

	int n = M->n; // number of events
	int num_edges = count_edges(M->P, n); // count existing edges
	int accepted = 0;

	// initialize RNG to use in iterations
	srand(rand_seed);
	gsl_rng *RNG = gsl_rng_alloc(gsl_rng_taus); 
	gsl_rng_set(RNG, rand_seed);  // seed rng

	long double posterior_o = likelihood_computation(scheme, M, D, num_unique);
	// printf("posterior_o = %.10Lg\n", posterior_o);
	
	//make a copy of model, thetas, and error for use in iterations
	model M_p;
	model_copy(M, &M_p);
	int* theta_types_p = get_int_array(n);
	memcpy(theta_types_p, theta_types, sizeof(int) * n);

	// open a log file to write to
	FILE *samplelog = NULL;
	if (strcmp(output, "") != 0) samplelog = open_file(output, "w");

	double msp, msp_r, tp, tp_r;
	msp = msp_r = tp = tp_r = 1;
	int total_iterations = burn_in + number_samples;
	int count_failures = 0;
	int largest = 0;
	int failure_threshold = total_iterations / 1000; // consider lowering this for future use with 100k iterations
	int actual_iterations = 0;
	int actual_iterations_threshold = total_iterations * 1000;
	int suboptimal = 0;
	long double best = -FLT_MAX;
	int** best_poset = get_int_matrix(n, n);
	int* best_theta_types = get_int_array(n);

	for (int i = 0; i< total_iterations; i++) {
		actual_iterations++;
		if (actual_iterations >= actual_iterations_threshold) {
			printf("Actual number of MCMC iterations reached threshold of %d iterations. Stopping.\n", actual_iterations_threshold);
			break;
		}
		double ms;
		if (count_failures >= failure_threshold) { // do an aggressive move
			ms = 2;
			new_total_ordering(&M_p, theta_types_p, D, num_unique);
			printf("Number of consecutive failures = %d: exceeds threshold. On iteration i = %d.\n", count_failures, i);
		} else {
			ms = gsl_ran_flat(RNG, 0, 1); // generate random number in [0, 1] from uniform distribution.

			if (ms < exchange) { // exchange
				if (count_edges(M->P, n) == 0) {
					i--;
					count_failures++;
					continue;
				}
				int ti, tj, x, y;
				get_two_rand(&x, &y, n);
				propose_event_exchange_move(&M_p, &tp, &ti, &tj, x, y);
				msp = msp_r = fraction_exchange;
				tp = tp_r = 1;
			} else if (ms < reincarnate) {
				if (!propose_reincarnation_move(&M_p, D, num_unique)){
					i--;
					continue;
				}
			
				msp = msp_r = fraction_reincarnation;
				tp = tp_r = 1;
			} else if (ms < birth) {
				int x = rand();
				if (!propose_new_cover_relation(&M_p, &tp, D, num_unique, x)) {
					i--; 
					count_failures++;
					continue;
				}

				tp_r = 1.0 / (double) count_edges(M_p.P, M_p.n); // tp_r for delete cover move
				msp = fraction_birth;
				msp_r = fraction_death;
			} else if (ms < death) {
				int x = rand();
				if (!propose_delete_cover_relation(&M_p, &tp, x)) {
					i--;
					count_failures++;
					continue;
				}
				tp_r = get_new_cover_move_tp(&M_p, D, num_unique);
				msp = fraction_death;
				msp_r = fraction_birth;
			} else if (ms < tr_birth) {
				int x = rand();
				if (!propose_new_transitive_closure(&M_p, &tp, x)) {
					i--;
					count_failures++;
					continue;
				}
				tp_r = get_transitive_closure_delete_tp(&M_p);
				msp = fraction_tr_birth;
				msp_r = fraction_tr_death;
			} else if (ms < tr_death) {
				int x = rand();
				if (!propose_delete_transitive_closure(&M_p, &tp, x)) {
					i--;
					count_failures++;
					continue;
				}
				tp_r = get_transitive_closure_new_tp(&M_p);
				msp = fraction_tr_death;
				msp_r = fraction_tr_birth;
			} else { // relocate theta types
				bool valid = false;
				for (int j = 0; j < n; j++) {
					if (M_p.num_pa[j] >= 2) {
						valid = true;
					}
				}
				if (!valid) {
					i--;
					count_failures++;
					continue;
				}
			 	int x = rand() % n;
			 	relocate_theta_type(theta_types_p, x);
				tp = tp_r = 1;
				msp = msp_r = 1;
			} 
			set_theta_types(&M_p, theta_types, theta_types_p);	

			// check that new model supports a cumulative model of cancer. If not, retry. 
			if (!cumulative_probs_maintained(&M_p, D, theta_types_p, num_unique)) {
				free_model(M_p, false);
				model_copy(M, &M_p);
				memcpy(theta_types_p, theta_types, sizeof(int)*n);
				i--;
				count_failures++;
				continue;
			}
		}

		long double posterior_p = likelihood_computation(scheme, &M_p, D, num_unique);
		long double MH_ratio = posterior_o * msp_r * tp_r / (posterior_p * msp * tp);
		long double alpha = min(1, MH_ratio);
		double randy = gsl_ran_flat (RNG, 0, 1);

		// decide if take the move if not. It would be significantly more efficient if we could pre-compute alpha, since there is a chance that 
		// we do a bunch of computations and then don't make the move. 
		if (alpha > randy) {
			if (posterior_o > posterior_p) suboptimal++;
			if (posterior_p > best) {
				best = posterior_p;
				for(int j = 0; j < n; j++) {
					for(int k = 0; k < n; k++) {
						best_poset[j][k] = M_p.P[j][k];
					}
				}
				memcpy(best_theta_types, theta_types_p, n * sizeof(int));
			}
			make_move(&posterior_o, posterior_p, M, &M_p, theta_types, theta_types_p, ms, &num_edges);

			if (count_failures > largest) largest = count_failures;
			if (i > burn_in) {
				accepted++;
				// if (!(accepted % record_ith)) {
				// 	printf("sample number %d: posterior %.10Lf number edges %d\n", i, posterior_o, num_edges);
				// }
			}
		} else { //restore original model, theta_types
			free_model(M_p, false);
			model_copy(M, &M_p);
			memcpy(theta_types_p, theta_types, sizeof(int)*n);
		}

		count_failures = 0;
		if(i > burn_in && !(i % record_ith) ) { // log samples
			if (!(i % 1000)) {
				printf("sample number %d. posterior %.10Lg. number edges %d.\n", i, posterior_o, num_edges);
			}
		}

	}

	log_results(best_poset, best_theta_types, accepted, actual_iterations, samplelog, suboptimal, largest, n);
		
	if (samplelog != NULL) fclose(samplelog);

	free_int_matrix(M->P, n);
	M->P = best_poset;
	memmove(theta_types, best_theta_types, sizeof(int) * n);

	
  gsl_rng_free (RNG);
	free_model(M_p, false);
	free_int_array(best_theta_types);
	free_int_array(theta_types_p);
}


#endif
