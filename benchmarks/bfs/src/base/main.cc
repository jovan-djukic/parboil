#include <deque>
#include <iostream>
#include <parboil.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INFINITY 2147483647 // 2^31-1
#define UP_LIMIT 16677216   // 2^24
#define WHITE 16677217
#define GRAY 16677218
#define BLACK 16677221

typedef struct {
	int edges_start;
	int edges_count;
} Node;

typedef struct {
	int node;
	int cost;
} Edge;

void bfs ( Node *nodes, Edge *edges, int *colors, int *costs, int source ) {
	std::deque<int> wavefront;

	wavefront.push_back ( source );

	colors[source] = GRAY;

	while ( !wavefront.empty ( ) ) {
		int index = wavefront.front ( );
		wavefront.pop_front ( );

		for ( int i = nodes[index].edges_start; i < ( nodes[index].edges_count + nodes[index].edges_start ); ++i ) {
			int id = edges[i].node;
			if ( colors[id] == WHITE ) {
				costs[id] = costs[index] + 1;
				wavefront.push_back ( id );
				colors[id] = GRAY;
			}
		}

		colors[index] = BLACK;
	}
}

int main ( int argc, char **argv ) {
    struct pb_TimerSet timers;
    pb_InitializeTimerSet ( &timers );

    struct pb_Parameters *parameters = pb_ReadParameters ( &argc, argv );
	if ( ( parameters->inpFiles[0] == NULL ) || ( parameters->inpFiles[1] != NULL ) ) {
		fprintf ( stderr, "Expecting one input filename\n" );
		exit ( -1 );
	}

	pb_SwitchToTimer ( &timers, pb_TimerID_IO );
	FILE *file = fopen ( parameters->inpFiles[0], "r" );
	if ( !file ) {
		printf ( "Error Reading graph file\n" );
		return 1;
	}

    int number_of_nodes = 0;
	fscanf ( file, "%d", &number_of_nodes );

	// allocate memory
    Node *nodes  = ( Node* ) malloc ( sizeof ( Node ) * number_of_nodes );
    int  *colors = ( int* ) malloc ( sizeof ( int ) * number_of_nodes );

	// initalize the memory
	for ( int i = 0; i < number_of_nodes; ++i ) {
		int	 start, edgeno;
		fscanf ( file, "%d %d", &start, &edgeno );

		nodes[i].edges_start = start;
		nodes[i].edges_count = edgeno;

		colors[i] = WHITE;
	}

	// read the source node from the file
	int source = 0;
	fscanf ( file, "%d", &source );

	// read edges from file
    int edge_list_size = 0;
	fscanf ( file, "%d", &edge_list_size );

	Edge *edges = ( Edge* ) malloc ( sizeof ( Edge ) * edge_list_size );

	for ( int i = 0; i < edge_list_size; ++i ) {
		int	 id, cost;
		fscanf ( file, "%d", &id );
		fscanf ( file, "%d", &cost );

		edges[i].node = id;
		edges[i].cost = cost;
	}

	if ( file ) {
		fclose ( file );
	}

	// allocate memory for the result	
	int *costs = ( int* ) malloc ( sizeof ( int ) * number_of_nodes );

	for ( int i = 0; i < number_of_nodes; ++i ) {
		costs[i] = INFINITY;
	}
	costs[source] = 0;

	// run bfs
	pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );
	bfs ( nodes, edges, colors, costs, source );
	pb_SwitchToTimer ( &timers, pb_TimerID_IO );

	if ( parameters->outFile != NULL ) {
		FILE *file = fopen ( parameters->outFile, "w" );
		fprintf ( file, "%d\n", number_of_nodes );

		for ( int i = 0; i < number_of_nodes; ++i ) {
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

