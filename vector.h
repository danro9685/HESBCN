/**
 * Basic implementation for dynamically resizable integer array class (vector/arraylist).
 * Kevin Chen, July 7th, 2016.
 */

#ifndef VECTOR_H
#define VECTOR_H

#include <stdbool.h>
#include "matrix.h"

#define DEFAULT_START_SIZE 128 // Arbitrarily chosen power of 2

typedef struct {
	int size; // number of elements in array
	int size_cap; // current cap for number of elements
	int* elems; //pointer to elements
	int pad; //padding for memory friendly size.
} vector;

typedef struct {
	int size;
	int size_cap;
	int** elems;
	int pad;
} ptr_vector;


int get_size(vector* v) {
	return v->size;
}

vector* new_vector() {
	vector* v = malloc(sizeof(vector));
	if (v == NULL) {
		fprintf(stderr, "Error: Out of memory.\n");
		exit(1);
	}
	v->size = 0;
	v->size_cap = DEFAULT_START_SIZE;
	v->elems = calloc(DEFAULT_START_SIZE, sizeof(int));
	if (v->elems == NULL) {
		fprintf(stderr, "Error: Out of memory.\n");
		exit(1);
	}
	return v;
}

ptr_vector* new_ptr_vector() {
	ptr_vector* v = malloc(sizeof(ptr_vector));
	if (v == NULL) {
		fprintf(stderr, "Error: Out of memory.\n");
		exit(1);
	}
	v->size = 0;
	v->size_cap = DEFAULT_START_SIZE;
	v->elems = calloc(DEFAULT_START_SIZE, sizeof(int*));
	if (v->elems == NULL) {
		fprintf(stderr, "Error: Out of memory.\n");
		exit(1);
	}
	return v;
}


int* free_vector(vector* v, bool free_elems) {
	if (!free_elems) {
		int* ret = v->elems;
		free(v);
		return ret;
	}
	free(v->elems);
	free(v);
	return NULL;
}

int** free_ptr_vector(ptr_vector* v, bool free_elems) {
	if (!free_elems) {
		int** ret = v->elems;
		free(v);
		return ret;
	}
	for (int i = 0; i< v->size; i++) {
		free_int_array(v->elems[i]);
	}
	free(v->elems);
	free(v);
	return NULL;
}

void increase_size(vector* v) {
	v->size_cap *= 2;
	v->elems = realloc(v->elems, v->size_cap * sizeof(int));
	if (v->elems == NULL) {
		fprintf(stderr, "Error: out of memory.\n");
		exit(1);
	}
}

void increase_size_ptr(ptr_vector* v) {
	v->size_cap *=2;
	v->elems = realloc(v->elems, v->size_cap * sizeof(int*));
	if (v->elems == NULL) {
		fprintf(stderr, "Error: out of memory.\n");
		exit(1);
	}
}

void push_back(vector* v, int payload) {
	if (v->size >= v->size_cap) {
		increase_size(v);
	}
	v->elems[v->size] = payload;
	v->size++;
}

void push_back_ptr(ptr_vector* v, int* payload) {
	if(v->size >= v->size_cap) {
		increase_size_ptr(v);
	}
	v->elems[v->size] = payload;
	v->size++;
}


//TODO: insert, remove functions


#endif
