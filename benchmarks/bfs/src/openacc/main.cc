/***************************************************************************
 *cr
 *cr						(C) Copyright 2007 The Board of Trustees of the
 *cr												University of Illinois
 *cr												 All Rights Reserved
 *cr
 ***************************************************************************/
/*
	Implementing Breadth first search on CUDA using algorithm given in DAC'10
	paper "An Effective GPU Implementation of Breadth-First Search"

	Copyright (c) 2010 University of Illinois at Urbana-Champaign.
	All rights reserved.

	Permission to use, copy, modify and distribute this software and its documentation for
	educational purpose is hereby granted without fee, provided that the above copyright
	notice and this permission notice appear in all copies of this software and that you do
	not sell the software.

	THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,EXPRESS, IMPLIED OR
	OTHERWISE.

	Author: Lijiuan Luo (lluo3@uiuc.edu)
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <parboil.h>

#define INFINITY 2147483647 // 2^31-1
#define UP_LIMIT 16677216   // 2^24
#define WHITE 16677217
#define GRAY 16677218
#define BLACK 16677221

typedef struct node {
	int edges_start;
	int edges_count;
} Node;

typedef struct edge {
	int node;
	int cost;
} Edge;

#pragma acc routine seq
extern int atomicExch ( int *, int );

#pragma acc routine seq
extern int atomicMin ( int *, int );

int breadthFirstSearch ( Node *nodes, Edge *edges, int *colors, int *costs, int source, int number_of_nodes, int number_of_edges ) {
	int	counter = 0;

	int *queue0	 = ( int * ) calloc ( number_of_nodes, sizeof ( int ) );
	int *queue1	 = ( int * ) calloc ( number_of_nodes, sizeof ( int ) );

	if ( queue0 == NULL || queue1 == NULL ) {
		return 1;
	}

	// switch to know which queue to use
	int k = 0;

	// enqueue source node
    queue0[counter] = source;

    counter = (counter + 1) % number_of_nodes;

    colors[source] = GRAY;

	#pragma acc data copyin ( queue0 [0:number_of_nodes], nodes [0:number_of_nodes], edges [0:number_of_edges], colors [0:number_of_nodes] ) copy ( costs [0:number_of_nodes] ) create ( queue1 [0:number_of_nodes] )
	{
		while ( counter > 0 ) {
			int new_counter = 0;

			#pragma acc parallel loop copy ( new_counter )
			for (  int i = 0; i < counter; ++i ) {
				int *input_queue  = ( k & 1 ) == 0 ? queue0 : queue1;
				int *output_queue = ( k & 1 ) == 0 ? queue1 : queue0;

				int n = input_queue[i];
				
				colors[n] = BLACK;
				for ( int e = nodes[n].edges_start; e < ( nodes[n].edges_start + nodes[n].edges_count ); ++e ) {
					int id = edges[e].node;

					int current_cost = edges[e].cost + costs[n];

					int original_cost = atomicMin ( &costs[id], current_cost );

					if ( original_cost > current_cost ) {
						int current_color = atomicExch ( &colors[id], GRAY );

						if ( current_color == WHITE ) {
							costs[id] = current_cost;

							int next = 0;

							#pragma acc atomic capture
							{
								next = new_counter;
								new_counter += 1;
							}

							output_queue[next] = id;
						}
					}
				}
			}

			counter = new_counter;
			k++;
		}
	}

	free ( queue0 );
	free ( queue1 );

	return 0;
}

int main ( int argc, char **argv ) {
    struct pb_TimerSet timers;
    pb_InitializeTimerSet ( &timers );

	struct pb_Parameters *parameters = pb_ReadParameters ( &argc, argv );
	if ( ( parameters->inpFiles[0] == NULL ) || ( parameters->inpFiles[1] != NULL ) ) {
		fprintf ( stderr, "Expecting one input filename\n" );
		return 1;
	}

	pb_SwitchToTimer ( &timers, pb_TimerID_IO );

	FILE *input_file = fopen ( parameters->inpFiles[0], "r" );
	if ( input_file == NULL ) {
		fprintf ( stderr, "Error Reading graph file\n" );
		return 2;
	}

	// read number of nodes
    int number_of_nodes = 0;
    fscanf ( input_file, "%d", &number_of_nodes );

	// allocate node array and color (visisted) array
    Node *nodes = ( Node * ) malloc ( sizeof ( Node ) * number_of_nodes );
    int  *color = ( int  *) malloc ( sizeof ( int ) * number_of_nodes );


	for ( unsigned int i = 0; i < number_of_nodes; i++ ) {
    	int edges_start = 0, edges_count = 0;
		fscanf ( input_file, "%d %d", &edges_start, &edges_count );

	    nodes[i].edges_start = edges_start;
	    nodes[i].edges_count = edges_count;

	    color[i] = WHITE;
    }

	// read the source node from the file
	int source;
	fscanf ( input_file, "%d", &source );

	// read the size of the edge array
    int number_of_edges = 0;
    fscanf ( input_file, "%d", &number_of_edges );

	// allocate array for edges
	Edge *edges = ( Edge * ) malloc ( sizeof ( Edge ) * number_of_edges );

	for ( int i = 0; i < number_of_edges; ++i ) {
		int id	 = 0;
		int cost = 0;

		fscanf ( input_file, "%d", &id );
		fscanf ( input_file, "%d", &cost );

	    edges[i].node = id;
	    edges[i].cost = cost;
    }

	fclose ( input_file );

	// allocate memory for the costs of paths (result)
	int *costs = ( int * ) malloc ( sizeof ( int ) * number_of_nodes );

	for ( int i = 0; i < number_of_nodes; ++i ) {
		costs[i] = INFINITY;
	}

	costs[source] = 0;

	// printf("start cpu version\n");
	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );
	breadthFirstSearch ( nodes, edges, color, costs, source, number_of_nodes, number_of_edges );
	pb_SwitchToTimer ( &timers, pb_TimerID_IO );

	if ( parameters->outFile != NULL ) {
		FILE *output_file = fopen ( parameters->outFile, "w" );

		fprintf ( output_file, "%d\n", number_of_nodes );

		for ( unsigned int i = 0; i < number_of_nodes; ++i ) {
			fprintf ( output_file, "%d %d\n", i, costs[i] );
		}

		fclose ( output_file );
	}

	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

	// cleanup memory
	free ( nodes );
	free ( edges );
	free ( color );
	free ( costs );
	pb_SwitchToTimer ( &timers, pb_TimerID_NONE );
	pb_PrintTimerSet ( &timers );
	pb_FreeParameters ( parameters );

	return 0;
}