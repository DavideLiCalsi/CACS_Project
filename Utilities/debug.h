#ifndef DEBUG_HEADER
#define DEBUG_HEADER

#include "../Matrix/BinaryMatrix.h"
#include "../Utilities/dataStructures.h"
#include "../Tree/VectorList.h"

#define DEBUG_PRINT 1UL
#define DEBUG_MATRICES 1UL<<1
#define DEBUG_SETS 1UL<<2
#define DEBUG_VLISTS 1UL<<3

// Or the desired flags to configure the debug level
#define DEBUGGING 0

#if (DEBUGGING & DEBUG_PRINT) == DEBUG_PRINT
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

#if (DEBUGGING & DEBUG_MATRICES) == DEBUG_MATRICES
#define PRINTMATRIX(x) printMatrix(x)
#define PRINTMATRIX_PTR(x) printMatrix(*x)
#else
#define PRINTMATRIX(x)
#define PRINTMATRIX_PTR(x)
#endif

#if (DEBUGGING & DEBUG_SETS) == DEBUG_SETS
#define PRINTSET(s) printSet(s)
#else
#define PRINTSET(s)
#endif

#if (DEBUGGING & DEBUG_VLISTS) == DEBUG_VLISTS
#define PRINTVLIST(l) VectorList_print(v)
#else
#define PRINTVLIST(l)
#endif

#endif