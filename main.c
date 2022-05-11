//#include "Matrix.h"
#include "Matrix/BinaryMatrix.h"

int main(){
    
    BinMatrix* H;
    int H_vect[128]=
    {
        0,0,1,1,1,0,1,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,
        1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1
    };

    H = buildMatrix(H_vect,2,64);
    printMatrix(*H);

    addRows(H,0,1);

    printMatrix(*H);
    int rows[3]={0,1,3};
    BinMatrix* code = sampleFromMatrix(rows,1,*H, MATRIX_SAMPLE_ROWS);
    printf("%d\n",codeWeight(*code));

    return 0;
}