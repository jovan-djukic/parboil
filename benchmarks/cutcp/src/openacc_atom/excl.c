/***************************************************************************
 *cr
 *cr            (C) Copyright 2008-2010 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <math.h>
// #include "atom.h"
// #include "cutoff.h"
// 
// #define CELL_LENGTH 4.f
// #define INVERSE_CELL_LENGTH ( 1.f / CELL_LENGTH )
// 
// extern int cpu_remove_exclusions(
// 	Lattice *lattice,	  /* the lattice */
// 	float cutoff,		  /* cutoff distance */
// 	Atoms *atomsStructure /* array of atoms */
// ) {
// 	int numberOfGridCellsX = lattice->dim.nx;
// 	int numberOfGridCellsY = lattice->dim.ny;
// 	int numberOfGridCellsZ = lattice->dim.nz;
// 
// 	float xLow = lattice->dim.lo.x;
// 	float yLow = lattice->dim.lo.y;
// 	float zLow = lattice->dim.lo.z;
// 
// 	float gridspacing = lattice->dim.h;
// 
// 	int numberOfAtoms = atomsStructure->size;
// 	Atom *atoms = atomsStructure->atoms;
// 
// 	const float cutoff2 = cutoff * cutoff;
// 	const float inverseCuttof2 = 1.f / cutoff2;
// 	const float inverseGridspacing = 1.f / gridspacing;
// 
// 	const int radius = (int)ceilf(cutoff * inverseGridspacing) - 1;
// 
// 	Vec3 minimumExtent, maximumExtent; /* Extent of atom bounding box */
// 	/* find min and max extent */
// 	get_atom_extent ( &minimumExtent, &maximumExtent, atomsStructure );
// 
// 	/* number of cells in each dimension */
// 	int numberOfAtomCellsX = ( int ) floorf ( ( maximumExtent.x - minimumExtent.x ) * INVERSE_CELL_LENGTH ) + 1;
// 	int numberOfAtomCellsY = ( int ) floorf ( ( maximumExtent.y - minimumExtent.y ) * INVERSE_CELL_LENGTH ) + 1;
// 	int numberOfAtomCellsZ = ( int ) floorf ( ( maximumExtent.z - minimumExtent.z ) * INVERSE_CELL_LENGTH ) + 1;
// 	int numberOfAtomCells  = numberOfAtomCellsX * numberOfAtomCellsY * numberOfAtomCellsZ;
// 
// 	/* allocate for cursor link list implementation */
// 	int *first = ( int* ) malloc ( numberOfAtomCells * sizeof ( int ) );
// 
// 	for ( int i = 0; i < numberOfAtomCells; ++i ) {
// 		first[i] = -1;
// 	}
// 
// 	int *next = ( int* ) malloc ( numberOfAtoms * sizeof ( int ) );
// 	for ( int i = 0; i < numberOfAtoms; ++i ) {
// 		next[i] = -1;
// 	}
// 
// 	/* geometric hashing */
// 	for ( int i = 0; i < numberOfAtoms; ++i ) {
// 		if ( 0 == atoms[i].q ) {
// 			continue; /* skip any non-contributing atoms */
// 		}
// 
// 		int atomCellX = ( int ) floorf ( ( atoms[i].x - minimumExtent.x ) * INVERSE_CELL_LENGTH );
// 		int atomCellY = ( int ) floorf ( ( atoms[i].y - minimumExtent.y ) * INVERSE_CELL_LENGTH );
// 		int atomCellZ = ( int ) floorf ( ( atoms[i].z - minimumExtent.z ) * INVERSE_CELL_LENGTH );
// 
// 		int gridIndex = ( atomCellZ * numberOfAtomCellsY + atomCellY ) * numberOfAtomCellsX + atomCellX;
// 
// 		next[i] = first[gridIndex];
// 		first[gridIndex] = i;
// 	}
// 
// 	/* traverse the grid cells */
// 	for ( int gridIndex = 0; gridIndex < numberOfAtomCells; ++gridIndex ) {
// 		for ( int atomIndex = first[gridIndex]; atomIndex != -1; atomIndex = next[atomIndex] ) {
// 			float atomX = atoms[atomIndex].x - xLow;
// 			float atomY = atoms[atomIndex].y - yLow;
// 			float atomZ = atoms[atomIndex].z - zLow;
// 			float atomQ = atoms[atomIndex].q;
// 
// 			/* find closest grid point with position less than or equal to atom */
// 			int gridCellX = ( int ) ( atomX * inverseGridspacing );
// 			int gridCellY = ( int ) ( atomY * inverseGridspacing );
// 			int gridCellZ = ( int ) ( atomZ * inverseGridspacing );
// 
// 			/* find extent of surrounding box of grid points */
// 			int gridCellXStart = gridCellX - radius;
// 			int gridCellXEnd   = gridCellX + radius + 1;
// 			int gridCellYStart = gridCellY - radius;
// 			int gridCellYEnd   = gridCellY + radius + 1;
// 			int gridCellZStart = gridCellZ - radius;
// 			int gridCellZEnd   = gridCellZ + radius + 1;
// 
// 			/* trim box edges so that they are within grid point lattice */
// 			if ( gridCellXStart < 0 ) { gridCellXStart = 0; }
// 			if ( gridCellYStart < 0 ) { gridCellYStart = 0; }
// 			if ( gridCellZStart < 0 ) { gridCellZStart = 0; }
// 
// 			if ( gridCellXEnd >= numberOfGridCellsX ) { gridCellXEnd = numberOfGridCellsX - 1; }
// 			if ( gridCellYEnd >= numberOfGridCellsY ) { gridCellYEnd = numberOfGridCellsY - 1; }
// 			if ( gridCellZEnd >= numberOfGridCellsZ ) { gridCellZEnd = numberOfGridCellsZ - 1; }
// 
// 			/* loop over surrounding grid points */
// 			for ( int k = gridCellZStart; k <= gridCellZEnd; ++k ) {
// 				for ( int j = gridCellYStart; j <= gridCellYEnd; ++j ) {
// 					for ( int i = gridCellXStart; i <= gridCellXEnd; ++i ) {
// 						float dz = k * gridspacing - atomZ;
// 						float dy = j * gridspacing - atomY;
// 						float dx = i * gridspacing - atomX;
// 						float r2 = dx * dx + dy * dy + dz * dz;
// 
// 						if ( r2 >= cutoff2 ) {
// 							continue;
// 						}
// 
// 						int index = ( k * numberOfGridCellsY + j ) * numberOfGridCellsX + i;
// 
// 						lattice->lattice[index] += 0;
// 					}
// 				}
// 			}
// 		}	
// 	}
// 
// 	free ( next );
// 	free ( first );
// 
// 	return 0;
// }
// 
// // #define BIN_LENGTH 8.f
// // #define ATOMS_PER_BIN 128
// #define BIN_LENGTH 4.f
// #define ATOMS_PER_BIN 16
// #define INVERSE_BIN_LENGTH ( 1.f / BIN_LENGTH )
// #define FLOATS_PER_ATOM 4
// 
// extern int remove_exclusions(
// 	Lattice *lattice,	  /* the lattice */
// 	float cutoff,		  /* cutoff distance */
// 	Atoms *atomsStructure /* array of atoms */
// ) {
// 	int numberOfGridCellsX = lattice->dim.nx;
// 	int numberOfGridCellsY = lattice->dim.ny;
// 	int numberOfGridCellsZ = lattice->dim.nz;
// 
// 	float xLow = lattice->dim.lo.x;
// 	float yLow = lattice->dim.lo.y;
// 	float zLow = lattice->dim.lo.z;
// 
// 	float gridspacing = lattice->dim.h;
// 
// 	int numberOfAtoms = atomsStructure->size;
// 	Atom *atoms = atomsStructure->atoms;
// 
// 	const float cutoff2 = cutoff * cutoff;
// 	const float inverseCuttof2 = 1.f / cutoff2;
// 	const float inverseGridspacing = 1.f / gridspacing;
// 
// 	const int radius = ( int ) ceilf ( cutoff * inverseGridspacing ) - 1;
// 
// 	Vec3 minimumExtent, maximumExtent; /* Extent of atom bounding box */
// 	/* find min and max extent */
// 	get_atom_extent ( &minimumExtent, &maximumExtent, atomsStructure );
// 
// 	/* number of cells in each dimension */
// 	int numberOfAtomBinsX = ( int ) floorf ( ( maximumExtent.x - minimumExtent.x ) * INVERSE_BIN_LENGTH ) + 1;
// 	int numberOfAtomBinsY = ( int ) floorf ( ( maximumExtent.y - minimumExtent.y ) * INVERSE_BIN_LENGTH ) + 1;
// 	int numberOfAtomBinsZ = ( int ) floorf ( ( maximumExtent.z - minimumExtent.z ) * INVERSE_BIN_LENGTH ) + 1;
// 	int numberOfAtomBins  = numberOfAtomBinsX * numberOfAtomBinsY * numberOfAtomBinsZ;
// 
// 	/* allocate for cursor link list implementation */
// 	Atoms *extraAtoms = ( Atoms* ) calloc ( 1, sizeof ( Atoms ) );
// 	extraAtoms->size = 0;
// 	extraAtoms->atoms = ( Atom* ) calloc ( numberOfAtoms, sizeof ( Atom ) );
// 
// 	float *bins = ( float* ) calloc ( numberOfAtomBins * ATOMS_PER_BIN * FLOATS_PER_ATOM, sizeof ( float ) );
// 	int *binCounters = ( int* ) calloc ( numberOfAtomBins, sizeof ( int ) );
// 
// 	/* bin hashing */
// 	for ( int i = 0; i < numberOfAtoms; ++i ) {
// 		if ( 0 == atoms[i].q ) {
// 			continue; /* skip any non-contributing atoms */
// 		}
// 
// 		int atomBinX = ( int ) floorf ( ( atoms[i].x - minimumExtent.x ) * INVERSE_BIN_LENGTH );
// 		int atomBinY = ( int ) floorf ( ( atoms[i].y - minimumExtent.y ) * INVERSE_BIN_LENGTH );
// 		int atomBinZ = ( int ) floorf ( ( atoms[i].z - minimumExtent.z ) * INVERSE_BIN_LENGTH );
// 
// 		int binIndex = ( atomBinZ * numberOfAtomBinsY + atomBinY ) * numberOfAtomBinsX + atomBinX;
// 
// 		int binSize = binCounters[binIndex];
// 		if ( binSize < ATOMS_PER_BIN ) {
// 			int index = ( binIndex * ATOMS_PER_BIN + binSize ) * FLOATS_PER_ATOM;
// 
// 			bins[index + 0] = atoms[i].x;
// 			bins[index + 1] = atoms[i].y;
// 			bins[index + 2] = atoms[i].z;
// 			bins[index + 3] = atoms[i].q;
// 
// 			binCounters[binIndex]++;
// 		} else {
// 			// TODO: this might be a problem
// 			extraAtoms->atoms[extraAtoms->size] = atoms[i];
// 			extraAtoms->size++;
// 		}
// 	}
// 
//     int binLatticesSize = ( ( numberOfGridCellsX * numberOfGridCellsY * numberOfGridCellsZ ) + 7 ) & ~7;
// 	float *binLattices = ( float* ) calloc ( binLatticesSize, sizeof ( float ) );
// 	for ( int i = 0; i < binLatticesSize; ++i ) {
// 		binLattices[i] = -1;
// 	}
// 
// 	/* traverse bins */
// 	#pragma acc parallel loop collapse(2) copy(binLattices[0:binLatticesSize]) copyin(xLow, yLow, zLow, radius, numberOfGridCellsX, numberOfGridCellsY, numberOfGridCellsZ, gridspacing, bins[0:numberOfAtomBins*ATOMS_PER_BIN*FLOATS_PER_ATOM], cutoff2, inverseCuttof2) async
// 	for ( int binIndex = 0; binIndex < numberOfAtomBins; ++binIndex ) {
// 		for ( int atomIndex = 0; atomIndex < ATOMS_PER_BIN; ++atomIndex ) {
// 			int index = ( binIndex * ATOMS_PER_BIN + atomIndex ) * FLOATS_PER_ATOM;
// 
// 			float atomX = bins[index + 0] - xLow;
// 			float atomY = bins[index + 1] - yLow;
// 			float atomZ = bins[index + 2] - zLow;
// 			float atomQ = bins[index + 3];
// 
// 			/* find closest grid point with position less than or equal to atom */
// 			int gridCellX = ( int ) ( atomX * inverseGridspacing );
// 			int gridCellY = ( int ) ( atomY * inverseGridspacing );
// 			int gridCellZ = ( int ) ( atomZ * inverseGridspacing );
// 
// 			/* find extent of surrounding box of grid points */
// 			int gridCellXStart = gridCellX - radius;
// 			int gridCellXEnd   = gridCellX + radius + 1;
// 			int gridCellYStart = gridCellY - radius;
// 			int gridCellYEnd   = gridCellY + radius + 1;
// 			int gridCellZStart = gridCellZ - radius;
// 			int gridCellZEnd   = gridCellZ + radius + 1;
// 
// 			/* trim box edges so that they are within grid point lattice */
// 			if ( gridCellXStart < 0 ) { gridCellXStart = 0; }
// 			if ( gridCellYStart < 0 ) { gridCellYStart = 0; }
// 			if ( gridCellZStart < 0 ) { gridCellZStart = 0; }
// 
// 			if ( gridCellXEnd >= numberOfGridCellsX ) { gridCellXEnd = numberOfGridCellsX - 1; }
// 			if ( gridCellYEnd >= numberOfGridCellsY ) { gridCellYEnd = numberOfGridCellsY - 1; }
// 			if ( gridCellZEnd >= numberOfGridCellsZ ) { gridCellZEnd = numberOfGridCellsZ - 1; }
// 
// 			/* loop over surrounding grid points */
// 			for ( int gridZ = gridCellZStart; gridZ <= gridCellZEnd; ++gridZ ) {
// 				for ( int gridY = gridCellYStart; gridY <= gridCellYEnd; ++gridY ) {
// 					for ( int gridX = gridCellXStart; gridX <= gridCellXEnd; ++gridX ) {
// 						float dz = gridZ * gridspacing - atomZ;
// 						float dy = gridY * gridspacing - atomY;
// 						float dx = gridX * gridspacing - atomX;
// 						float r2 = dx * dx + dy * dy + dz * dz;
// 
// 						if ( r2 >= cutoff2 ) {
// 							continue;
// 						}
// 
// 						int gridIndex = ( gridZ * numberOfGridCellsY + gridY ) * numberOfGridCellsX + gridX;
// 
// 						#pragma acc atomic write
// 						binLattices[gridIndex] = 0;
// 					}
// 				}
// 			}
// 		}	
// 	}
// 
// 	cpu_remove_exclusions ( lattice, cutoff, extraAtoms );
// 
// 	#pragma acc wait
// 
// 	for ( int i = 0; i < binLatticesSize; ++i ) {
// 		if ( binLattices[i] == 0 ) {
// 			lattice->lattice[i] = 0;
// 		}
// 	}
// 
// 	free ( binLattices );
// 	free ( extraAtoms->atoms );
// 	free ( extraAtoms );
// 	free ( bins );
// 	free ( binCounters );
// 
// 	return 0;
// }
// 

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

#define CELLEN      4.f
#define INV_CELLEN  (1.f/CELLEN)

extern int remove_exclusions(
    Lattice *lattice,                  /* the lattice */
    float cutoff,                      /* exclusion cutoff distance */
    Atoms *atoms                       /* array of atoms */
    )
{
  int nx = lattice->dim.nx;
  int ny = lattice->dim.ny;
  int nz = lattice->dim.nz;
  float xlo = lattice->dim.lo.x;
  float ylo = lattice->dim.lo.y;
  float zlo = lattice->dim.lo.z;
  float gridspacing = lattice->dim.h;
  Atom *atom = atoms->atoms;

  const float a2 = cutoff * cutoff;
  const float inv_gridspacing = 1.f / gridspacing;
  const int radius = (int) ceilf(cutoff * inv_gridspacing) - 1;
    /* lattice point radius about each atom */

  int n;
  int i, j, k;
  int ia, ib, ic;
  int ja, jb, jc;
  int ka, kb, kc;
  int index;
  int koff, jkoff;

  float x, y, z, q;
  float dx, dy, dz;
  float dz2, dydz2, r2;
  float e;
  float xstart, ystart;

  float *pg;

  int gindex;
  int ncell, nxcell, nycell, nzcell;
  int *first, *next;
  float inv_cellen = INV_CELLEN;
  Vec3 minext, maxext;

  /* find min and max extent */
  get_atom_extent(&minext, &maxext, atoms);

  /* number of cells in each dimension */
  nxcell = (int) floorf((maxext.x-minext.x) * inv_cellen) + 1;
  nycell = (int) floorf((maxext.y-minext.y) * inv_cellen) + 1;
  nzcell = (int) floorf((maxext.z-minext.z) * inv_cellen) + 1;
  ncell = nxcell * nycell * nzcell;

  /* allocate for cursor link list implementation */
  first = (int *) malloc(ncell * sizeof(int));
  for (gindex = 0;  gindex < ncell;  gindex++) {
    first[gindex] = -1;
  }
  next = (int *) malloc(atoms->size * sizeof(int));
  for (n = 0;  n < atoms->size;  n++) {
    next[n] = -1;
  }

  /* geometric hashing */
  for (n = 0;  n < atoms->size;  n++) {
    if (0==atom[n].q) continue;  /* skip any non-contributing atoms */
    i = (int) floorf((atom[n].x - minext.x) * inv_cellen);
    j = (int) floorf((atom[n].y - minext.y) * inv_cellen);
    k = (int) floorf((atom[n].z - minext.z) * inv_cellen);
    gindex = (k*nycell + j)*nxcell + i;
    next[n] = first[gindex];
    first[gindex] = n;
  }

  /* traverse the grid cells */
  for (gindex = 0;  gindex < ncell;  gindex++) {
    for (n = first[gindex];  n != -1;  n = next[n]) {
      x = atom[n].x - xlo;
      y = atom[n].y - ylo;
      z = atom[n].z - zlo;
      q = atom[n].q;

      /* find closest grid point with position less than or equal to atom */
      ic = (int) (x * inv_gridspacing);
      jc = (int) (y * inv_gridspacing);
      kc = (int) (z * inv_gridspacing);

      /* find extent of surrounding box of grid points */
      ia = ic - radius;
      ib = ic + radius + 1;
      ja = jc - radius;
      jb = jc + radius + 1;
      ka = kc - radius;
      kb = kc + radius + 1;

      /* trim box edges so that they are within grid point lattice */
      if (ia < 0)   ia = 0;
      if (ib >= nx) ib = nx-1;
      if (ja < 0)   ja = 0;
      if (jb >= ny) jb = ny-1;
      if (ka < 0)   ka = 0;
      if (kb >= nz) kb = nz-1;

      /* loop over surrounding grid points */
      xstart = ia*gridspacing - x;
      ystart = ja*gridspacing - y;
      dz = ka*gridspacing - z;
      for (k = ka;  k <= kb;  k++, dz += gridspacing) {
        koff = k*ny;
        dz2 = dz*dz;

        dy = ystart;
        for (j = ja;  j <= jb;  j++, dy += gridspacing) {
          jkoff = (koff + j)*nx;
          dydz2 = dy*dy + dz2;

          dx = xstart;
          index = jkoff + ia;
          pg = lattice->lattice + index;

          for (i = ia;  i <= ib;  i++, pg++, dx += gridspacing) {
            r2 = dx*dx + dydz2;

	    /* If atom and lattice point are too close, set the lattice value
	     * to zero */
            if (r2 < a2) *pg = 0;
          }
        }
      } /* end loop over surrounding grid points */

    } /* end loop over atoms in a gridcell */
  } /* end loop over gridcells */

  /* free memory */
  free(next);
  free(first);

  return 0;
}
