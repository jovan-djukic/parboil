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

#define CELL_LENGTH 4.f
#define INVERSE_CELL_LENGTH ( 1.f / CELL_LENGTH )

extern int cpu_compute_cutoff_potential_lattice(
	Lattice *lattice,	  /* the lattice */
	float cutoff,		  /* cutoff distance */
	Atoms *atomsStructure /* array of atoms */
) {
	int numberOfGridCellsX = lattice->dim.nx;
	int numberOfGridCellsY = lattice->dim.ny;
	int numberOfGridCellsZ = lattice->dim.nz;

	float xLow = lattice->dim.lo.x;
	float yLow = lattice->dim.lo.y;
	float zLow = lattice->dim.lo.z;

	float gridspacing = lattice->dim.h;

	int numberOfAtoms = atomsStructure->size;
	Atom *atoms = atomsStructure->atoms;

	const float cutoff2 = cutoff * cutoff;
	const float inverseCuttof2 = 1.f / cutoff2;
	const float inverseGridspacing = 1.f / gridspacing;

	const int radius = (int)ceilf(cutoff * inverseGridspacing) - 1;

	Vec3 minimumExtent, maximumExtent; /* Extent of atom bounding box */
	/* find min and max extent */
	get_atom_extent ( &minimumExtent, &maximumExtent, atomsStructure );

	/* number of cells in each dimension */
	int numberOfAtomCellsX = ( int ) floorf ( ( maximumExtent.x - minimumExtent.x ) * INVERSE_CELL_LENGTH ) + 1;
	int numberOfAtomCellsY = ( int ) floorf ( ( maximumExtent.y - minimumExtent.y ) * INVERSE_CELL_LENGTH ) + 1;
	int numberOfAtomCellsZ = ( int ) floorf ( ( maximumExtent.z - minimumExtent.z ) * INVERSE_CELL_LENGTH ) + 1;
	int numberOfAtomCells  = numberOfAtomCellsX * numberOfAtomCellsY * numberOfAtomCellsZ;

	/* allocate for cursor link list implementation */
	int *first = ( int* ) malloc ( numberOfAtomCells * sizeof ( int ) );

	for ( int i = 0; i < numberOfAtomCells; ++i ) {
		first[i] = -1;
	}

	int *next = ( int* ) malloc ( numberOfAtoms * sizeof ( int ) );
	for ( int i = 0; i < numberOfAtoms; ++i ) {
		next[i] = -1;
	}

	/* geometric hashing */
	for ( int i = 0; i < numberOfAtoms; ++i ) {
		if ( 0 == atoms[i].q ) {
			continue; /* skip any non-contributing atoms */
		}

		int atomCellX = ( int ) floorf ( ( atoms[i].x - minimumExtent.x ) * INVERSE_CELL_LENGTH );
		int atomCellY = ( int ) floorf ( ( atoms[i].y - minimumExtent.y ) * INVERSE_CELL_LENGTH );
		int atomCellZ = ( int ) floorf ( ( atoms[i].z - minimumExtent.z ) * INVERSE_CELL_LENGTH );

		int gridIndex = ( atomCellZ * numberOfAtomCellsY + atomCellY ) * numberOfAtomCellsX + atomCellX;

		next[i] = first[gridIndex];
		first[gridIndex] = i;
	}

	/* traverse the grid cells */
	for ( int gridIndex = 0; gridIndex < numberOfAtomCells; ++gridIndex ) {
		for ( int atomIndex = first[gridIndex]; atomIndex != -1; atomIndex = next[atomIndex] ) {
			float atomX = atoms[atomIndex].x - xLow;
			float atomY = atoms[atomIndex].y - yLow;
			float atomZ = atoms[atomIndex].z - zLow;
			float atomQ = atoms[atomIndex].q;

			/* find closest grid point with position less than or equal to atom */
			int gridCellX = ( int ) ( atomX * inverseGridspacing );
			int gridCellY = ( int ) ( atomY * inverseGridspacing );
			int gridCellZ = ( int ) ( atomZ * inverseGridspacing );

			/* find extent of surrounding box of grid points */
			int gridCellXStart = gridCellX - radius;
			int gridCellXEnd   = gridCellX + radius + 1;
			int gridCellYStart = gridCellY - radius;
			int gridCellYEnd   = gridCellY + radius + 1;
			int gridCellZStart = gridCellZ - radius;
			int gridCellZEnd   = gridCellZ + radius + 1;

			/* trim box edges so that they are within grid point lattice */
			if ( gridCellXStart < 0 ) { gridCellXStart = 0; }
			if ( gridCellYStart < 0 ) { gridCellYStart = 0; }
			if ( gridCellZStart < 0 ) { gridCellZStart = 0; }

			if ( gridCellXEnd >= numberOfGridCellsX ) { gridCellXEnd = numberOfGridCellsX - 1; }
			if ( gridCellYEnd >= numberOfGridCellsY ) { gridCellYEnd = numberOfGridCellsY - 1; }
			if ( gridCellZEnd >= numberOfGridCellsZ ) { gridCellZEnd = numberOfGridCellsZ - 1; }

			/* loop over surrounding grid points */
			for ( int k = gridCellZStart; k <= gridCellZEnd; ++k ) {
				for ( int j = gridCellYStart; j <= gridCellYEnd; ++j ) {
					for ( int i = gridCellXStart; i <= gridCellXEnd; ++i ) {
						float dz = k * gridspacing - atomZ;
						float dy = j * gridspacing - atomY;
						float dx = i * gridspacing - atomX;
						float r2 = dx * dx + dy * dy + dz * dz;

						if ( r2 >= cutoff2 ) {
							continue;
						}

						int index = ( k * numberOfGridCellsY + j ) * numberOfGridCellsX + i;

						float s = (1.f - r2 * inverseCuttof2);
						float e = atomQ * (1 / sqrtf(r2)) * s * s;

						lattice->lattice[index] += e;
					}
				}
			}
		}	
	}

	free ( next );
	free ( first );

	return 0;
}

// #define BIN_LENGTH 8.f
// #define ATOMS_PER_BIN 128
#define BIN_LENGTH 4.f
#define ATOMS_PER_BIN 16
#define INVERSE_BIN_LENGTH ( 1.f / BIN_LENGTH )
#define FLOATS_PER_ATOM 4

extern int openacc_compute_cutoff_potential_lattice(
	Lattice *lattice,	  /* the lattice */
	float cutoff,		  /* cutoff distance */
	Atoms *atomsStructure /* array of atoms */
) {
	int numberOfGridCellsX = lattice->dim.nx;
	int numberOfGridCellsY = lattice->dim.ny;
	int numberOfGridCellsZ = lattice->dim.nz;

	float xLow = lattice->dim.lo.x;
	float yLow = lattice->dim.lo.y;
	float zLow = lattice->dim.lo.z;

	float gridspacing = lattice->dim.h;

	int numberOfAtoms = atomsStructure->size;
	Atom *atoms = atomsStructure->atoms;

	const float cutoff2 = cutoff * cutoff;
	const float inverseCuttof2 = 1.f / cutoff2;
	const float inverseGridspacing = 1.f / gridspacing;

	const int radius = ( int ) ceilf ( cutoff * inverseGridspacing ) - 1;

	Vec3 minimumExtent, maximumExtent; /* Extent of atom bounding box */
	/* find min and max extent */
	get_atom_extent ( &minimumExtent, &maximumExtent, atomsStructure );

	/* number of cells in each dimension */
	int numberOfAtomBinsX = ( int ) floorf ( ( maximumExtent.x - minimumExtent.x ) * INVERSE_BIN_LENGTH ) + 1;
	int numberOfAtomBinsY = ( int ) floorf ( ( maximumExtent.y - minimumExtent.y ) * INVERSE_BIN_LENGTH ) + 1;
	int numberOfAtomBinsZ = ( int ) floorf ( ( maximumExtent.z - minimumExtent.z ) * INVERSE_BIN_LENGTH ) + 1;
	int numberOfAtomBins  = numberOfAtomBinsX * numberOfAtomBinsY * numberOfAtomBinsZ;

	/* allocate for cursor link list implementation */
	Atoms *extraAtoms = ( Atoms* ) calloc ( 1, sizeof ( Atoms ) );
	extraAtoms->size = 0;
	extraAtoms->atoms = ( Atom* ) calloc ( numberOfAtoms, sizeof ( Atom ) );

	float *bins = ( float* ) calloc ( numberOfAtomBins * ATOMS_PER_BIN * FLOATS_PER_ATOM, sizeof ( float ) );
	int *binCounters = ( int* ) calloc ( numberOfAtomBins, sizeof ( int ) );

	/* bin hashing */
	for ( int i = 0; i < numberOfAtoms; ++i ) {
		if ( 0 == atoms[i].q ) {
			continue; /* skip any non-contributing atoms */
		}

		int atomBinX = ( int ) floorf ( ( atoms[i].x - minimumExtent.x ) * INVERSE_BIN_LENGTH );
		int atomBinY = ( int ) floorf ( ( atoms[i].y - minimumExtent.y ) * INVERSE_BIN_LENGTH );
		int atomBinZ = ( int ) floorf ( ( atoms[i].z - minimumExtent.z ) * INVERSE_BIN_LENGTH );

		int binIndex = ( atomBinZ * numberOfAtomBinsY + atomBinY ) * numberOfAtomBinsX + atomBinX;

		int binSize = binCounters[binIndex];
		if ( binSize < ATOMS_PER_BIN ) {
			int index = ( binIndex * ATOMS_PER_BIN + binSize ) * FLOATS_PER_ATOM;

			bins[index + 0] = atoms[i].x;
			bins[index + 1] = atoms[i].y;
			bins[index + 2] = atoms[i].z;
			bins[index + 3] = atoms[i].q;

			binCounters[binIndex]++;
		} else {
			// TODO: this might be a problem
			extraAtoms->atoms[extraAtoms->size] = atoms[i];
			extraAtoms->size++;
		}
	}

    int binLatticesSize = ( ( numberOfGridCellsX * numberOfGridCellsY * numberOfGridCellsZ ) + 7 ) & ~7;
	float *binLattices = ( float* ) calloc ( binLatticesSize, sizeof ( float ) );

	/* traverse bins */
	#pragma acc parallel loop collapse(2) copy(binLattices[0:binLatticesSize]) copyin(xLow, yLow, zLow, radius, numberOfGridCellsX, numberOfGridCellsY, numberOfGridCellsZ, gridspacing, bins[0:numberOfAtomBins*ATOMS_PER_BIN*FLOATS_PER_ATOM], cutoff2, inverseCuttof2) async
	for ( int binIndex = 0; binIndex < numberOfAtomBins; ++binIndex ) {
		for ( int atomIndex = 0; atomIndex < ATOMS_PER_BIN; ++atomIndex ) {
			int index = ( binIndex * ATOMS_PER_BIN + atomIndex ) * FLOATS_PER_ATOM;

			float atomX = bins[index + 0] - xLow;
			float atomY = bins[index + 1] - yLow;
			float atomZ = bins[index + 2] - zLow;
			float atomQ = bins[index + 3];

			/* find closest grid point with position less than or equal to atom */
			int gridCellX = ( int ) ( atomX * inverseGridspacing );
			int gridCellY = ( int ) ( atomY * inverseGridspacing );
			int gridCellZ = ( int ) ( atomZ * inverseGridspacing );

			/* find extent of surrounding box of grid points */
			int gridCellXStart = gridCellX - radius;
			int gridCellXEnd   = gridCellX + radius + 1;
			int gridCellYStart = gridCellY - radius;
			int gridCellYEnd   = gridCellY + radius + 1;
			int gridCellZStart = gridCellZ - radius;
			int gridCellZEnd   = gridCellZ + radius + 1;

			/* trim box edges so that they are within grid point lattice */
			if ( gridCellXStart < 0 ) { gridCellXStart = 0; }
			if ( gridCellYStart < 0 ) { gridCellYStart = 0; }
			if ( gridCellZStart < 0 ) { gridCellZStart = 0; }

			if ( gridCellXEnd >= numberOfGridCellsX ) { gridCellXEnd = numberOfGridCellsX - 1; }
			if ( gridCellYEnd >= numberOfGridCellsY ) { gridCellYEnd = numberOfGridCellsY - 1; }
			if ( gridCellZEnd >= numberOfGridCellsZ ) { gridCellZEnd = numberOfGridCellsZ - 1; }

			/* loop over surrounding grid points */
			for ( int gridZ = gridCellZStart; gridZ <= gridCellZEnd; ++gridZ ) {
				for ( int gridY = gridCellYStart; gridY <= gridCellYEnd; ++gridY ) {
					for ( int gridX = gridCellXStart; gridX <= gridCellXEnd; ++gridX ) {
						float dz = gridZ * gridspacing - atomZ;
						float dy = gridY * gridspacing - atomY;
						float dx = gridX * gridspacing - atomX;
						float r2 = dx * dx + dy * dy + dz * dz;

						if ( r2 >= cutoff2 ) {
							continue;
						}

						int gridIndex = ( gridZ * numberOfGridCellsY + gridY ) * numberOfGridCellsX + gridX;

						float s = ( 1.f - r2 * inverseCuttof2 );
						float e = atomQ * ( 1 / sqrtf ( r2 ) ) * s * s;

						#pragma acc atomic update
						binLattices[gridIndex] += e;
					}
				}
			}
		}	
	}

	cpu_compute_cutoff_potential_lattice ( lattice, cutoff, extraAtoms );

	#pragma acc wait

	for ( int i = 0; i < binLatticesSize; ++i ) {
		lattice->lattice[i] += binLattices[i];
	}

	free ( binLattices );
	free ( extraAtoms->atoms );
	free ( extraAtoms );
	free ( bins );
	free ( binCounters );

	return 0;
}
