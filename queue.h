/**
 * Linked-list implementation of FIFO queue. 
 * List nodes have integer payloads.
 * Kevin Chen. July 6th, 2016
 */

#ifndef QUEUE_H
#define QUEUE_H

#include <stdbool.h>


typedef struct {
	void* next;
	int payload;
} node;

typedef struct {
    node* head; // front of queue
    node* tail; // back of queue
} queue;

queue* new_queue() {
	queue* q = malloc(sizeof(queue));
	if (q == NULL) {
		fprintf(stderr, "Error: out of memory.\n");
		exit(1);
	}
	q->head = NULL;
	q->tail = NULL;
	return q;
}

bool is_empty(queue* q) {
	return q->head == NULL;
}

void enqueue(queue *q, int payload) {
	node* curr = malloc(sizeof(node));
	if (curr == NULL) {
		fprintf(stderr, "Error: out of memory.\n");
		exit(1);
	}
	curr->next = NULL;
	curr->payload = payload;
	if (is_empty(q)) {
		q->head = curr;
	} else {
		q->tail->next = curr;
	}
	q->tail = curr;
}

int dequeue(queue* q) {
	node* first = q->head;
	q->head = first->next;
	return first->payload;
}

#endif
