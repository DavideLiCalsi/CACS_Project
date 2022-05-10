//#include "Matrix.h"
#include "Matrix/BinaryMatrix.h"

int main(){
    
    BinMatrix* H;
    int H_vect[16]=
    {
        1,0,1,1,
        0,1,0,1,
        0,0,1,0,
        0,0,0,0
    };

    H = buildMatrix(H_vect,4,4);
    putElement(H,3,3,1);

    printf("Det(H): %d\n", determinant(*H));
    BinMatrix* inv = inverse(*H);

    printMatrix(*product(*inv,*H));

    return 0;
}