/***********************************************************************************
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
	Revised for Parboil 2 Benchmark Suite by: Geng Daniel Liu (gengliu2@illinois.edu)
 ************************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <parboil.h>
#include <deque>
#include <iostream>

#include "config.h"

typedef int2 Node;
typedef int2 Edge;

#include "kernel.cu"

const int h_top = 1;
const int zero	= 0;

int main ( int argc, char** argv ) {
	struct pb_TimerSet timers;
	pb_InitializeTimerSet ( &timers );

	struct pb_Parameters* parameters = pb_ReadParameters ( &argc, argv );
	if ( ( parameters->inpFiles[0] == NULL ) || ( parameters->inpFiles[1] != NULL ) ) {
	fprintf ( stderr, "Expecting one input filename\n" );
	exit ( -1 );
	}

	pb_SwitchToTimer ( &timers, pb_TimerID_IO );

	// Read graph from a file
	FILE *file = fopen ( parameters->inpFiles[0], "r" );
	if ( !file ) {
		printf ( "Error Reading graph file\n" );
		return 0;
	}

	// Read the number of nodes in the graph
	int num_of_nodes = 0;
	fscanf ( file, "%d", &num_of_nodes );

	// Allocate host memory
    Node* h_graph_nodes = ( Node* ) malloc ( sizeof ( Node ) * num_of_nodes );
    int*  color         = ( int* ) malloc ( sizeof ( int ) * num_of_nodes );
	// Initalize the memory
	for ( int i = 0; i < num_of_nodes; ++i ) {
    	int	 start, edgeno;
		fscanf ( file, "%d %d", &start, &edgeno );

		h_graph_nodes[i].x = start;
		h_graph_nodes[i].y = edgeno;

	    color[i] = WHITE;
    }
	// Read the source node and the number of edges in graph from the file
    int source       = 0;
    int num_of_edges = 0;
    fscanf ( file, "%d", &source );
	fscanf ( file, "%d", &num_of_edges );

	Edge* h_graph_edges = ( Edge* ) malloc ( sizeof ( Edge ) * num_of_edges );
		for ( int i = 0; i < num_of_edges; ++i ) {
		int	 id, cost;
		fscanf ( file, "%d", &id );
		fscanf ( file, "%d", &cost );

		h_graph_edges[i].x = id;
		h_graph_edges[i].y = cost;
	}

	if ( file ) {
		fclose ( file );
	}

	// Allocate memory for the result on host side
	int* h_cost = ( int* ) malloc ( sizeof ( int ) * num_of_nodes );
	for ( int i = 0; i < num_of_nodes; ++i ) {
		h_cost[i] = INFINITY;
	}
	h_cost[source] = 0;

	pb_SwitchToTimer ( &timers, pb_TimerID_COPY );

	// Copy the Node list to device memory
	Node* d_graph_nodes;
	cudaMalloc ( ( void** ) &d_graph_nodes, sizeof ( Node ) * num_of_nodes );
	cudaMemcpy ( d_graph_nodes, h_graph_nodes, sizeof ( Node ) * num_of_nodes, cudaMemcpyHostToDevice );

	// Copy the Edge List to device Memory
	Edge* d_graph_edges;
	cudaMalloc ( ( void** ) &d_graph_edges, sizeof ( Edge ) * num_of_edges );
	cudaMemcpy ( d_graph_edges, h_graph_edges, sizeof ( Edge ) * num_of_edges, cudaMemcpyHostToDevice );

	int* d_color;
	cudaMalloc ( ( void** ) &d_color, sizeof ( int ) * num_of_nodes );
	cudaMemcpy ( d_color, color, sizeof ( int ) * num_of_nodes, cudaMemcpyHostToDevice );

	int* d_cost;
	cudaMalloc ( ( void** ) &d_cost, sizeof ( int ) * num_of_nodes );
	cudaMemcpy ( d_cost, h_cost, sizeof ( int ) * num_of_nodes, cudaMemcpyHostToDevice );

	int* d_q1;
	int* d_q2;
	cudaMalloc ( ( void** ) &d_q1, sizeof ( int ) * num_of_nodes );
	cudaMalloc ( ( void** ) &d_q2, sizeof ( int ) * num_of_nodes );

	int* tail;
	cudaMalloc ( ( void** ) &tail, sizeof ( int ) );

	int* front_cost_d;
	cudaMalloc ( ( void** ) &front_cost_d, sizeof ( int ) );

	// Bind the texture memory with global memory
	cudaBindTexture ( 0, g_graph_node_ref, d_graph_nodes, sizeof ( Node ) * num_of_nodes );
	cudaBindTexture ( 0, g_graph_edge_ref, d_graph_edges, sizeof ( Edge ) * num_of_edges );

	printf ( "Starting GPU kernel\n" );
	cudaThreadSynchronize ( );
	pb_SwitchToTimer ( &timers, pb_TimerID_KERNEL );

	cudaMemcpy ( tail, &h_top, sizeof ( int ), cudaMemcpyHostToDevice );
	cudaMemcpy ( &d_cost[source], &zero, sizeof ( int ), cudaMemcpyHostToDevice );

	cudaMemcpy ( &d_q1[0], &source, sizeof ( int ), cudaMemcpyHostToDevice );
	int num_t; // number of threads
	int k = 0; // BFS level index

	do {
		cudaMemcpy ( &num_t, tail, sizeof ( int ), cudaMemcpyDeviceToHost );
		cudaMemcpy ( tail, &zero, sizeof ( int ), cudaMemcpyHostToDevice );

		// frontier is empty
		if ( num_t == 0 ) { 			
			break;
		}

	    int num_of_blocks            = 1;
	    int num_of_threads_per_block = num_t;

	    if ( num_of_threads_per_block < NUM_BIN ) {
			num_of_threads_per_block = NUM_BIN;
		}

		if ( num_t > MAX_THREADS_PER_BLOCK ) {
	        num_of_blocks            = ( int ) ceil ( num_t / ( double ) MAX_THREADS_PER_BLOCK );
	        num_of_threads_per_block = MAX_THREADS_PER_BLOCK;
	    }

		// will call "BFS_in_GPU_kernel"
		if ( num_of_blocks == 1 ) {
			num_of_threads_per_block = MAX_THREADS_PER_BLOCK;
		}

		// will call "BFS_kernel_multi_blk_inGPU"
		if ( num_of_blocks > 1 && num_of_blocks <= NUM_SM ) {
			num_of_blocks = NUM_SM;
		}

		dim3 grid ( num_of_blocks, 1, 1 );
		dim3 threads ( num_of_threads_per_block, 1, 1 );

		if ( k % 2 == 0 ) {
			BFS_kernel<<<grid, threads>>> ( 
				d_q1,
				d_q2,
				d_graph_nodes,
				d_graph_edges,
				d_color,
				d_cost,
				num_t,
				tail,
				GRAY0,
				k 
			);
		} else {
			BFS_kernel<<<grid, threads>>> ( 
				d_q2,
				d_q1,
				d_graph_nodes,
				d_graph_edges,
				d_color,
				d_cost,
				num_t,
				tail,
				GRAY1,
				k 
			);
		}

		k++;
	} while ( 1 );

	cudaThreadSynchronize ( );
	pb_SwitchToTimer ( &timers, pb_TimerID_COPY );
	printf ( "GPU kernel done\n" );

	// Copy result from device to host
	cudaMemcpy ( h_cost, d_cost, sizeof ( int ) * num_of_nodes, cudaMemcpyDeviceToHost );
	cudaMemcpy ( color, d_color, sizeof ( int ) * num_of_nodes, cudaMemcpyDeviceToHost );
	cudaUnbindTexture ( g_graph_node_ref );
	cudaUnbindTexture ( g_graph_edge_ref );

	cudaFree ( d_graph_nodes );
	cudaFree ( d_graph_edges );
	cudaFree ( d_color );
	cudaFree ( d_cost );
	cudaFree ( tail );
	cudaFree ( front_cost_d );

	// Store the result into a file
	pb_SwitchToTimer ( &timers, pb_TimerID_IO );

	file = fopen ( parameters->outFile, "w" );

	fprintf ( file, "%d\n", num_of_nodes );
	for ( int i = 0; i < num_of_nodes; ++i ) {
		fprintf ( file, "%d %d\n", i, h_cost[i] );
	}
	fclose ( file );

	// cleanup memory
	free ( h_graph_nodes );
	free ( h_graph_edges );
	free ( color );
	free ( h_cost );
	pb_SwitchToTimer ( &timers, pb_TimerID_NONE );
	pb_PrintTimerSet ( &timers );
	pb_FreeParameters ( parameters );
	return 0;
}
