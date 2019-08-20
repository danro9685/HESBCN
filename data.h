#ifndef DATA_H_
#define DATA_H_


typedef struct {
	int *g;  // int array storing genotype. Consider using the binary compression? however, this limits to the number of bits avaialble. 
	bool is_compatible;  // true if compatible with model
	double* t; // time of each event.
	int** P;  // poset used for induced refinement of the poset. 
	int* subset;  // Binary array indicating compatibility of valid genotypes with this data point. 
	int count;  // number of observations of this type
} data;

void free_data_array(data* D, int N_u, const int n) {
	for (int i=0; i< N_u; i++) {
		free(D[i].g);
		free_int_matrix(D[i].P, n);
		free_double_array(D[i].t);
		free_int_array(D[i].subset);
	}
	free(D);
}

void print_data(data D, const int n) {
	printf("Data struct:\n");
	printf("Genotype: ");
	for (int i=0; i < n; i++) {
		printf("%d ", D.g[i]);
	}
	printf("\n");
	printf("Count = %d\n", D.count);
	printf("is compatible? ... ");
	if (D.is_compatible) {
		printf("true.\n");
	} else {
		printf("false.\n");
	}
}

void print_data_array(data* D, int N_u, const int n) {
	for (int i=0; i< N_u; i++) {
		print_data(D[i], n);
	}
}

#endif
