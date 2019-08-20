/**
 * Array implementation of FIFO queue.
 */


#ifndef QUEUE_H
#define QUEUE_H

#include <stdbool.h>
#define MAX_QUEUESIZE 32768

typedef struct {
	int head;
	int tail;
	int size;
	int q[MAX_QUEUESIZE];
} queue;

queue* new_queue() {
	queue* q = malloc(sizeof(queue));
	if (q == NULL) {
		fprintf(stderr, "Error: Out of memory.\n");
		exit(1);
	}
	q->head = 0;
	q->tail = MAX_QUEUESIZE - 1;
	q->size = 0;
	return q;
}

bool is_empty(queue* q) {
	return q->size <= 0;
}

void free_queue(queue *q) {
	free(q);
}

void enqueue(queue *q, int payload) {
	if (q->size >= MAX_QUEUESIZE) {
		fprintf(stderr, "Error: Queue overflow.\n");
		exit(1);
	}
	q->tail = (q->tail + 1) % MAX_QUEUESIZE; // mod to handle loop-around edge case
	q->q[q->tail] = payload;
	q->size++;
}

int dequeue(queue* q) {
	int ret = 0; //arbitrary initailization value
	if (is_empty(q)) {
		fprintf(stderr, "Attempt to dequeue empty queue");
	} else {
		ret = q->q[q->head];
		q->head = (q->head + 1) % MAX_QUEUESIZE; // mod to handle loop-around edge case
		q->size--;
	}
	return ret;
}

#endif
