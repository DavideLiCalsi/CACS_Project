#ifndef DEBUG_HEADER
#define DEBUG_HEADER

#include "../Matrix/BinaryMatrix.h"
#include "../Utilities/dataStructures.h"
#include "../Tree/VectorList.h"

// Set it to 1 to print debugging info
#define DEBUGGING 1

#if DEBUGGING==1
#define PRINTF(...) printf(__VA_ARGS__)
#define PRINTMATRIX(x) printMatrix(x)
#define PRINTMATRIX_PTR(x) printMatrix(*x)
#define PRINTSET(s) printSet(s)
#define PRINTVLIST(l) VectorList_print(v)
#else
#define PRINTF(...)
#define PRINTMATRIX(x)
#define PRINTMATRIX_PTR(x)
#define PRINTSET(s)
#define PRINTVLIST(l)
#endif

#endif