// def.h
#ifndef _DEF
#define _DEF

#include "myparameters.h"

#define NULL 0
#define FALSE 0
#define TRUE 1
typedef int BOOL;

#define NV  ((par.NVX)*par.NVY)
#define NNX ((par.NVX)+1)
#define NNY (par.NVY+1)
#define NN  (NNX*NNY)
#define NDOF (2*NN)

// loading
#define FORCE ((par.LOAD)*(par.VOXSIZE))

#define JCM ((par.NOSTICKJCM)*(par.VOXSIZE))  // cell-medium
#define JCC ((par.NOSTICKJCC)*(par.VOXSIZE)) // cell-cell

#define SQ05 sqrt(0.5)
#endif
