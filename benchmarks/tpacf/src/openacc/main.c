/***************************************************************************
 *cr
 *cr            (C) Copyright 2007 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#include "args.h"
#include "model.h"

int main ( int argc, char **argv ) {
    struct pb_TimerSet timers;
    pb_InitializeTimerSet ( &timers );

    struct pb_Parameters *params = pb_ReadParameters ( &argc, argv );

    options args;
    parse_args ( argc, argv, &args );

    pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );
    int nbins = ( int ) floor ( bins_per_dec *  ( log10 ( max_arcmin ) - log10 ( min_arcmin ) ) );
    size_t memsize = ( nbins + 2 ) * sizeof ( long long );

    // memory for bin boundaries
    float *binb =  ( float* ) malloc ( ( nbins + 1 ) * sizeof ( float ) );
    if ( binb == NULL ) {
        fprintf ( stderr, "Unable to allocate memory\n" );
        exit ( -1 );
    }

    for ( int k = 0; k < nbins + 1; ++k ) {
        binb[k] = cos ( pow ( 10, log10 ( min_arcmin ) + k * 1.0 / bins_per_dec ) / 60.0 * D2R );
    }

    // memory for DD
    long long *DD = ( long long* ) malloc ( memsize );
    if ( DD == NULL ) {
        fprintf ( stderr, "Unable to allocate memory\n" );
        exit ( -1 );
    }
    bzero ( DD, memsize );

    // memory for RR
    long long *RRS = ( long long* ) malloc ( memsize * args.random_count );
    if ( RRS == NULL ) {
        fprintf ( stderr, "Unable to allocate memory\n" );
        exit ( -1 );
    }
    bzero ( RRS, memsize * args.random_count );

    // memory for DR
    long long *DRS = ( long long* ) malloc ( memsize * args.random_count );
    if ( DRS == NULL ) {
        fprintf ( stderr, "Unable to allocate memory\n" );
        exit ( -1 );
    }
    bzero ( DRS, memsize * args.random_count );

    // memory for input data
    struct cartesian *data = ( struct cartesian* ) malloc ( args.npoints * sizeof ( struct cartesian ) );
    if ( data == NULL ) {
        fprintf ( stderr, "Unable to allocate memory for % data points  ( #1 )\n", args.npoints );
        return ( 0 );
    }

    // allocate memory for all random data
    struct cartesian *random = ( struct cartesian* ) malloc ( args.random_count * args.npoints * sizeof ( struct cartesian ) );
    if ( random == NULL ) {
        fprintf ( stderr, "Unable to allocate memory for % data points  ( #2 )\n", args.npoints );
        return ( 0 );
    }

    printf ( "Min distance: %f arcmin\n", min_arcmin );
    printf ( "Max distance: %f arcmin\n", max_arcmin );
    printf ( "Bins per dec: %i\n", bins_per_dec );
    printf ( "Total bins  : %i\n", nbins );

    // read data file
    pb_SwitchToTimer ( &timers, pb_TimerID_IO );
    int npd = readdatafile ( params->inpFiles[0], data, args.npoints );

    // read all random files
    for ( int rf = 0; rf < args.random_count; ++rf ) {
        struct cartesian *currentRandom = random + rf * args.npoints;
        int npr = readdatafile ( params->inpFiles[rf + 1], currentRandom, args.npoints );
        if ( npr != args.npoints ) {
            fprintf ( stderr, "Error: read %i random points out of %i in file %s\n", npr, args.npoints, params->inpFiles[rf + 1] );
            return  ( 0 );
        }
    }

    pb_SwitchToTimer ( &timers, pb_TimerID_COMPUTE );

    // #pragma acc data copyin(data[0:npd], random[0:(args.npoints * args.random_count)], binb[0:(nbins + 1)]) copy(DD[0:(nbins + 2)], RRS[0:(nbins + 2)], DRS[0:(nbins + 2)])
    // {
    //     // compute DD
    //     doCompute ( data, npd, NULL, 0, 1, DD, nbins, binb );

    //     // loop through random data files
    //     for ( int rf = 0; rf < args.random_count; ++rf ) {
    //         struct cartesian *currentRandom = random + rf * args.npoints;

    //         // compute RR
    //         doCompute ( currentRandom, args.npoints, NULL, 0, 1, RRS, nbins, binb );

    //         // compute DR
    //         doCompute ( data, npd, currentRandom, args.npoints, 0, DRS, nbins, binb );
    //     }
    // }
    doComputeWrapper ( args, data, random, DD, RRS, DRS, nbins, binb );


    // compute and output results
    FILE *outfile = fopen ( params->outFile, "w" );
    if ( outfile == NULL ) {
        fprintf ( stderr, "Unable to open output file %s for writing, assuming stdout\n", params->outFile );
        outfile = stdout;
    }

    pb_SwitchToTimer ( &timers, pb_TimerID_IO );
    for ( int k = 1; k < nbins + 1; ++k ) {
        fprintf ( outfile, "%d\n%d\n%d\n", DD[k], DRS[k], RRS[k] );
    }

    if ( outfile != stdout ) {
        fclose ( outfile );
    }

    // free memory
    free ( data );
    free ( random );
    free ( binb );
    free ( DD );
    free ( RRS );
    free ( DRS );

    pb_SwitchToTimer ( &timers, pb_TimerID_NONE );
    pb_PrintTimerSet ( &timers );
    pb_FreeParameters ( params );

    return ( 0 );
}
