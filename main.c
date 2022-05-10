//#include "Matrix.h"
#include "Matrix/BinaryMatrix.h"

int main(){
    
    BinMatrix* H;
    BinMatrix* a;
    int H_vect[105]={1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,
                    1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,
                    1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,
                    1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1,
                    1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1};
    int a_vect[7]={0,1,0,1,1,1,0};

    H = buildMatrix(H_vect,15,7);
    a = buildMatrix(a_vect,7,1);

    BinMatrix* prod = product(*H,*a);
    printMatrix(*prod);
    printMatrix(*H);

    swapRows(H,1,14);

    printMatrix(*H);

    return 0;
}