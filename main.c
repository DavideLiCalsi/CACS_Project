//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SplitSyndrome.h"


int main(){
    
    BST Xl=NULL;
    BST Xr=NULL;
    int H_arr[25]={
        1,1,0,1,1,
        1,1,1,1,0,
        0,1,0,0,1,
        0,1,0,0,1,
        1,0,1,1,1
    };

    BinMatrix* I5 = identityMatrix(5);
    BinMatrix* A = buildMatrix(H_arr,5,5);
    BinMatrix* H = concat(*transpose(*A),*I5,0);
    destroyMatrix(I5);
    destroyMatrix(A);

    buildLeftTable(3,2,*H,&Xl);
    buildRightTable(3,1,5,*H,&Xr);

    printTree(Xl,BST_COMPARISON_BINMATRIX);
    printTree(Xr,BST_COMPARISON_BINMATRIX);

    findEqualSize_u_m(8,2,1);
    
    return 0;
}