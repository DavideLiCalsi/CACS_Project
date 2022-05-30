//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
#include "Utilities/dataReader.h"
#include "SplitSyndrome/SplitSyndrome.h"

#define N 15

int main(){

    char path[100] = "./Utilities/info.txt";
    Info *info = readData(path);

    BinMatrix* A = transpose(*info->H_t);
    BinMatrix* I = identityMatrix((info->n)/2);
    
    BinMatrix* H = concat(*I,*A,0); // H=[I|A]
    destroyMatrix(A);
    destroyMatrix(I);

    BinMatrix* e=NULL;

    printMatrix(*H);
    SplitSyndrome(*H,*info->s,info->w,&e);

    if (e == NULL)
        printf("FAILURE! ERROR NOT FOUND\n");
    else
        printMatrix(*e);

    int vector[10];

    //iterateOverM_Vectors(10,3);

    
    return 0;
}