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

	Permission to use, copy, modify and distribute this software and its
	documentation for educational purpose is hereby granted without fee, provided
	that the above copyright notice and this permission notice appear in all
	copies of this software and that you do not sell the software.

	THE SOFTWARE IS PROVIDED "AS IS" AND WITHOUT WARRANTY OF ANY KIND,EXPRESS,
	IMPLIED OR OTHERWISE.

	Author: Lijiuan Luo (lluo3@uiuc.edu)
*/
#include <parboil.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <deque>
#include <iostream>

#define INFINITY 2147483647  // 2^31-1
#define WHITE 16677217
#define GRAY 16677218
#define BLACK 16677221

struct Node {
	int egdes_start;
	int edges_count;
};

struct Edge {
	int node;
	int cost;
};

void bfs ( Node* nodes, Edge* edges, int* colors, int* costs, int source ) {
	std::deque<int> wavefront;

	wavefront.push_back ( source );

	colors[source] = GRAY;

	while ( !wavefront.empty ( ) ) {
		int index = wavefront.front ( );
		wavefront.pop_front ( );

		#pragma omp parallel for
		for ( int i = nodes[index].egdes_start; i < ( nodes[index].edges_count + nodes[index].egdes_start ); ++i ) {
			int id = edges[i].node;

			if ( colors[id] == WHITE ) {
				costs[id] = costs[index] + 1;

				#pragma omp critical
				wavefront.push_back ( id );

				colors[id] = GRAY;
			}
		}

		colors[index] = BLACK;
	}
}

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
    FILE* file = fopen ( parameters->inpFiles[0], "r" );
    if ( !file ) {
	    printf ( "Error Reading graph file\n" );
	    return;
    }

    // Read number of nodes from file
    int no_of_nodes = 0;
    fscanf ( file, "%d", &no_of_nodes );

    // allocate host memory
    Node* nodes  = ( Node* ) malloc ( sizeof ( Node ) * no_of_nodes );
    int*  colors = ( int* ) malloc ( sizeof ( int ) * no_of_nodes );

    // initalize the memory
    for ( int i = 0; i < no_of_nodes; ++i ) {
	    int start, edgeno;
	    fscanf ( file, "%d %d", &start, &edgeno );

	    nodes[i].egdes_start = start;
	    nodes[i].edges_count = edgeno;

        colors[i] = WHITE;
    }

    // read the source node and edge list size from the file
    int source         = 0;
    int edge_list_size = 0;
    fscanf ( file, "%d", &source );
    fscanf ( file, "%d", &edge_list_size );

    Edge* edges = ( Edge* ) malloc ( sizeof ( Edge ) * edge_list_size );
    for ( int i = 0; i < edge_list_size; ++i ) {
	    int id, cost;
	    fscanf ( file, "%d", &id );
	    fscanf ( file, "%d", &cost );

	    edges[i].node = id;
	    edges[i].cost = cost;
    }

    if ( file ) {
		fclose ( file );
	}

    // allocate memory for the result
    int* costs = ( int* ) malloc ( sizeof ( int ) * no_of_nodes );
    for ( int i = 0; i < no_of_nodes; ++i ) {
	    costs[i] = INFINITY;
    }
    costs[source] = 0;

    pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );
    bfs ( nodes, edges, colors, costs, source );
    pb_SwitchToTimer ( &timers, pb_TimerID_IO );

    if ( parameters->outFile != NULL ) {
	    FILE* file = fopen ( parameters->outFile, "w" );
	    fprintf ( file, "%d\n", no_of_nodes );

	    for ( int i = 0; i < no_of_nodes; i++ ) {
			fprintf ( file, "%d %d\n", i, costs[i] );
		}

	    fclose ( file );
    }

    pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

    // cleanup memory
    free ( nodes );
    free ( edges );
    free ( colors );
    free ( costs );
    pb_SwitchToTimer ( &timers, pb_TimerID_NONE );
    pb_PrintTimerSet ( &timers );
    pb_FreeParameters ( parameters );

    return 0;
}
