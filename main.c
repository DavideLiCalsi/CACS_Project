//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SplitSyndrome.h"

int main(){

    int d=4;
    
    int A_array[25] = {

        0,0,1,0,0,
        0,0,0,1,0,
        1,0,0,0,0,
        0,1,0,1,1,
        0,0,1,1,1
    };

    int s_array[5]={
        0,1,1,0,0
    };

    BinMatrix* A = buildMatrix(A_array,5,5);
    BinMatrix* I = identityMatrix(5);
    BinMatrix* H = concat(*A,*I,0);
    destroyMatrix(A);
    destroyMatrix(I);

    BinMatrix* s = buildMatrix(s_array,1,5);

    BinMatrix* e=NULL;

    SplitSyndrome(*H,*s,d,&e);

    if (e == NULL)
        printf("FAILURE\n");
    else
        printMatrix(*e);
    
    return 0;
}