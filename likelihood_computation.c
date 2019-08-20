#ifndef LIKE_COMP
#define LIKE_COMP

#include "stdio.h"
#include "math.h"
#include "model.h"
#include "h-esbcn.h"


// num_params contains the number of events that occurred in the genotype, which allows for penalties to be applied given the number of events observed. 
long double dlik(int i, data* D, int num_unique) {
	long double res = 0;

	int count = 0;
	for (int j = 0; j < num_unique; j++) {
		if (D[j].g[i]) count += D[j].count;
	}
	if (count) res = (long double) count * log( (long double) count / (long double) NUM_SAMPLES_READ);
	return res;
}

bool parent_combo_present(int* parents, int num_parents, int* genotype, int* combo) {
	for (int i = 0; i < num_parents; i++) {
		if (genotype[parents[i]] != combo[i]) return false;
	}
	return true;
}

int* generate_combination(int num_parents, int* parents, int curr) {
	int* ret = get_int_array(num_parents);
	for (int i = 0; i < num_parents; i++) {
		if ((curr & (1 << (num_parents - 1 - i)))) ret[i] = 1;
	}
	return ret;
}

// Compute conditional likelihood of a node, using supplied data. 
long double cdlik(int num_parents, int* parents, data* D, int num_unique, int index, int* nparams) {
	long double res = 0;

	// represent all combinations of parents, as ordered in the parent set.
	int num_parent_combos = pow2(num_parents);
	int** parent_combos = malloc(num_parent_combos * sizeof(int*));
	int* parent_combo_counts = get_int_array(num_parent_combos); // marginals
	int** joint_freqs = get_int_matrix(2, num_parent_combos);

	// Generate combinations of parent set
	for (int i = 0; i < num_parent_combos; i++) {
		parent_combos[i] = generate_combination(num_parents, parents, i);
	}

	// Count marginal and joint frequencies
	for (int i = 0; i < num_unique; i++) {
		int* g = D[i].g;
		for (int j = 0; j < num_parent_combos; j++) {
			if (parent_combo_present(parents, num_parents, g, parent_combos[j])) {
				parent_combo_counts[j] += D[i].count;
				joint_freqs[g[index]][j] += D[i].count;
				break;
			}
		}
	}

	// calculate log-likelihood
	for (int i = 0; i < num_parent_combos; i++) {
		if (parent_combo_counts[i] != joint_freqs[0][i] + joint_freqs[1][i]) {
			printf("something's a bit off with the counting.\n");
			exit(1);
		}
		if (joint_freqs[1][i] != 0) {
			res += (long double) joint_freqs[1][i] * log( (long double) joint_freqs[1][i] / (long double) parent_combo_counts[i]);
		}
	}
	if (nparams != NULL)	{
		*nparams = num_parent_combos;
	}
	free_int_matrix(parent_combos, num_parent_combos);
	free_int_matrix(joint_freqs, 2);
	free_int_array(parent_combo_counts);
	return res;
}

// Use unconditioned likelihood if no parents, conditioned likelihood otherwise. 
long double loglik_dnode(model* M, int i, int* nparams, data* D, int num_unique, int num_parents) {
	long double loglik = 0;

	if (num_parents == 0) {
		loglik = dlik(i, D, num_unique);
	} else {
		loglik = cdlik(num_parents, M->pa[i], D, num_unique, i, nparams);
	}

	return loglik;
}

// Evaluates the likelihood of the supplied model, given the provided dataset and specified regularization scheme. 
long double likelihood_computation(char* scheme, model* M, data* D, int num_unique) {
	int num_events = M->n;
	long double prob = 0;

	// factored this way to prevent multiple checks on the scheme.
	if (scheme[0] == 'b') {
		double penalty_factor = log(NUM_SAMPLES_READ);
		int num_parameters = 0;
		for (int i = 0; i< num_events; i++) {
			int nparams = 1;
			long double this =loglik_dnode(M, i, &nparams, D, num_unique, M->num_pa[i]);
			// printf("%Lf ", this); 
			prob += 2 * this;
			num_parameters += nparams;
		}
		// prob *= 2;
		// printf("total prob = %Lf. total penalty = %f\n\n", prob, num_parameters * penalty_factor);
		prob -= num_parameters * penalty_factor; 
	} else if (scheme[0] == 'a') {
		double penalty_factor = 2;
		int num_parameters = 0;
		for (int i = 0; i< num_events; i++) {
			int nparams = 1;
			prob += loglik_dnode(M, i, &nparams, D, num_unique, M->num_pa[i]);
			num_parameters += nparams;
		}
		prob *= 2;
		prob -= num_parameters * penalty_factor; 
	}	else {
		int nparams;
		for (int i = 0; i < num_events; i++) {
			prob += loglik_dnode(M, i, &nparams, D, num_unique, M->num_pa[i]);
		}
	}

	return prob;
}

#endif
