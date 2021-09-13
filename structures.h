// file: structures.h
#ifndef _STRUCTURE
#define _STRUCTURE
#include <stdlib.h>
#include "def.h"

typedef struct
{
   	//double x, y;			// minimum coordinates
	int ctag;				// id of occupying cell, 0 if no cell
} VOX;

typedef struct
{
	//double x, y;			// coordinates
	double fx, fy;			// predefined forces
	double ux, uy;			// calculated displacements
	BOOL restrictx, restricty; // nodal restrictions
} NOD;


#endif
