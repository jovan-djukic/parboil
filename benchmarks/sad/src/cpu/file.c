/***************************************************************************
 *cr
 *cr            (C) Copyright 2007 The Board of Trustees of the
 *cr                        University of Illinois
 *cr                         All Rights Reserved
 *cr
 ***************************************************************************/

#include <stdio.h>
#include "file.h"

unsigned short read16u ( FILE *file ) {
	int value = fgetc ( file );
	value += fgetc ( file ) << 8;

	return value;
}

short read16i ( FILE *file ) {
	int value = fgetc ( file );
	value += fgetc ( file ) << 8;

	return value;
}

void write32u ( FILE *file, unsigned int value ) {
	putc ( value, file );
	putc ( value >> 8, file );
	putc ( value >> 16, file );
	putc ( value >> 24, file );
}

void write16u ( FILE *file, unsigned short value ) {
	putc ( value, file );
	putc ( value >> 8, file );
}

void write16i ( FILE *file, short value ) {
	putc ( value, file );
	putc ( value >> 8, file );
}
