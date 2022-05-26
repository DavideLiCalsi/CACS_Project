//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SplitSyndrome.h"

#define N 15

int main(){

    int d=7;
    
    int A_array[N*N] = {

        1,1,0,1,1,1,1,1,1,0,0,1,0,0,1,
0,1,0,0,1,1,0,1,1,1,0,1,1,1,0,
0,0,1,0,1,1,0,1,0,0,0,0,0,1,0,
0,1,1,0,1,1,0,1,0,1,1,0,1,1,0,
1,0,0,0,0,1,1,0,0,0,0,0,0,1,1,
0,0,1,1,1,1,1,0,1,0,1,1,0,0,0,
1,0,0,1,0,1,1,0,0,0,0,0,0,0,0,
0,0,1,0,1,0,0,0,0,0,0,1,0,0,0,
1,0,1,0,0,0,1,1,1,0,0,1,0,0,1,
0,1,1,1,0,0,0,0,0,1,1,0,1,0,0,
1,1,1,1,1,1,0,0,1,0,1,0,1,0,0,
1,1,1,1,0,1,1,1,0,0,0,1,0,0,0,
1,1,1,0,1,1,0,0,1,0,1,0,1,1,0,
0,1,1,1,1,0,1,0,0,0,0,1,1,1,0,
0,0,0,0,0,0,0,1,1,0,0,1,1,0,0
    };

    int s_array[N]={
        1,1,0,1,0,1,0,0,1,0,0,0,0,1,1
    };

    BinMatrix* A = buildMatrix(A_array,N,N);
    BinMatrix* A_t = transpose(*A);
    BinMatrix* I = identityMatrix(N);
    
    BinMatrix* H = concat(*I,*A_t,0); // H=[I|A]
    destroyMatrix(A);
    destroyMatrix(A_t);
    destroyMatrix(I);

    BinMatrix* s = buildMatrix(s_array,1,N);

    BinMatrix* e=NULL;

    printMatrix(*H);
    SplitSyndrome(*H,*s,d,&e);

    if (e == NULL)
        printf("FAILURE! ERROR NOT FOUND\n");
    else
        printMatrix(*e);

    int vector[10];

    //iterateOverM_Vectors(10,3);

    
    return 0;
}