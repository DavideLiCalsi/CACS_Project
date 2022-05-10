//#include "Matrix.h"
#include "Matrix/BinaryMatrix.h"

int main(){
    
    BinMatrix* H;
    BinMatrix* a;
    int H_vect[21]={1,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0,0,1,1,1,1};
    int a_vect[7]={0,1,0,1,1,1,0};

    H = buildMatrix(H_vect,3,7);
    a = buildMatrix(a_vect,7,1);

    BinMatrix* prod = product(*H,*a);
    printMatrix(*prod);
    printMatrix(*transpose(*H));

    return 0;
}