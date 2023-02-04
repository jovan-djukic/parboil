/***************************************************************************
 *cr
 *cr            (C) Copyright 2008-2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "atom.h"
#include "cutoff.h"

#undef DEBUG_PASS_RATE
#define CHECK_CYLINDER_CPU

#define CELLEN 4.f
#define INV_CELLEN (1.f / CELLEN)

extern int cpu_compute_cutoff_potential_lattice (
    Lattice *lattice, /* the lattice */
    float cutoff,     /* cutoff distance */
    Atoms *atoms      /* array of atoms */
)
{
    int nx = lattice->dim.nx;
    int ny = lattice->dim.ny;
    int nz = lattice->dim.nz;
    float xlo = lattice->dim.lo.x;
    float ylo = lattice->dim.lo.y;
    float zlo = lattice->dim.lo.z;
    float gridspacing = lattice->dim.h;
    int natoms = atoms->size;
    Atom *atom = atoms->atoms;

    const float a2 = cutoff * cutoff;
    const float inv_a2 = 1.f / a2;
    const float inv_gridspacing = 1.f / gridspacing;
    const int radius = ( int ) ceilf ( cutoff * inv_gridspacing ) - 1;
    /* lattice point radius about each atom */
    float inv_cellen = INV_CELLEN;
    Vec3 minext, maxext; /* Extent of atom bounding box */

    /* find min and max extent */
    get_atom_extent ( &minext, &maxext, atoms );

    /* number of cells in each dimension */
    int nxcell = ( int ) floorf ( ( maxext.x - minext.x ) * inv_cellen ) + 1;
    int nycell = ( int ) floorf ( ( maxext.y - minext.y ) * inv_cellen ) + 1;
    int nzcell = ( int ) floorf ( ( maxext.z - minext.z ) * inv_cellen ) + 1;
    int ncell = nxcell * nycell * nzcell;

    /* allocate for cursor link list implementation */
    int *first = ( int* ) malloc ( ncell * sizeof ( int ) );
    for ( int gindex = 0; gindex < ncell; ++gindex ) {
        first[gindex] = -1;
    }

    int *next = ( int* ) malloc ( natoms * sizeof ( int ) );
    for ( int n = 0; n < natoms; ++n ) {
        next[n] = -1;
    }

    /* geometric hashing */
    // this determines the index of a cell within the grid, two attoms may be in the same cell
    // this links the attoms in  the same cell, this must be attomic
    for ( int n = 0; n < natoms; ++n ) {
        if (0 == atom[n].q) {
            continue; /* skip any non-contributing atoms */
        }

        int i = ( int ) floorf ( ( atom[n].x - minext.x ) * inv_cellen );
        int j = ( int ) floorf ( ( atom[n].y - minext.y ) * inv_cellen );
        int k = ( int ) floorf ( ( atom[n].z - minext.z ) * inv_cellen );

        int gindex = (k * nycell + j) * nxcell + i;
        next[n] = first[gindex];
        first[gindex] = n;
    }

    float *lattices = lattice->lattice;
    int size = ( ( lattice->dim.nx * lattice->dim.ny * lattice->dim.nz ) + 7 ) & ~7;

    /* traverse the grid cells */
    #pragma acc parallel loop copyin(atom[0:natoms], first[0:ncell], next[0:natoms]) copy(lattices[0:size])
    for ( int gindex = 0; gindex < ncell; ++gindex ) {
        for ( int n = first[gindex]; n != -1; n = next[n] ) {
            float x = atom[n].x - xlo;
            float y = atom[n].y - ylo;
            float z = atom[n].z - zlo;
            float q = atom[n].q;

            /* find closest grid point with position less than or equal to atom */
            int ic = ( int ) ( x * inv_gridspacing );
            int jc = ( int ) ( y * inv_gridspacing );
            int kc = ( int ) ( z * inv_gridspacing );

            /* find extent of surrounding box of grid points */
            int ia = ic - radius;
            int ib = ic + radius + 1;
            int ja = jc - radius;
            int jb = jc + radius + 1;
            int ka = kc - radius;
            int kb = kc + radius + 1;

            /* trim box edges so that they are within grid point lattice */
            if ( ia < 0   ) { ia = 0;      }
            if ( ib >= nx ) { ib = nx - 1; }
            if ( ja < 0   ) { ja = 0;      }
            if ( jb >= ny ) { jb = ny - 1; }
            if ( ka < 0   ) { ka = 0;      }
            if ( kb >= nz ) { kb = nz - 1; }

            /* loop over surrounding grid points */
            float dxStart = ia * gridspacing - x;
            float dyStart = ja * gridspacing - y;
            float dzStart = ka * gridspacing - z;

            for ( int k = ka; k <= kb; ++k ) {
                for ( int j = ja; j <= jb; ++j ) {
                    for ( int i = ia; i <= ib; ++i ) {
                        float dx = dxStart + ( i - ia ) * gridspacing;
                        float dy = dyStart + ( j - ja ) * gridspacing;
                        float dz = dzStart + ( k - ka ) * gridspacing;
                        float r2 = dx * dx + dy * dy + dz * dz; 

                        if ( r2 >= a2 ) {
                            continue;
                        }

                        int index = ( k * ny + j ) * nx + i;

                        float s = ( 1.f - r2 * inv_a2 );
                        float e = q * ( 1 / sqrtf ( r2 ) ) * s * s;

                        #pragma acc atomic update
                        lattices[index] += e;
                    }
                }
            } /* end loop over surrounding grid points */

        } /* end loop over atoms in a gridcell */
    }     /* end loop over gridcells */

    /* free memory */
    free ( next );
    free ( first );

    return 0;
}
