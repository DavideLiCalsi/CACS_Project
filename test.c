#include <stdio.h>
#include "Matrix/BinaryMatrix.h"
#include "Utilities/utilities.h"


int main(){
    
    BinMatrix* G;
    
    int G_vect[28]=
    {
        0, 1, 1, 0, 0, 0, 1,
        1, 1, 1, 1, 1, 1, 1,
        1, 0, 0, 0, 1, 0, 1,
        1, 1, 0, 0, 0, 1, 0
    };
    G = buildMatrix(G_vect,4,7);

    //printMatrix(*G);
    //printMatrix(*standardizeGeneratorMatrix(G));


    BinMatrix* H;
    
    int H_vect[28]=
    {
        1, 1, 1, 0, 0, 1, 0,
        0, 1, 0, 0, 1, 1, 1,
        1, 0, 0, 1, 0, 1, 1
    };
    H = buildMatrix(H_vect,3,7);

    printMatrix(*H);
    printMatrix(*standardizeParityMatrix(H));

    return 0;
}