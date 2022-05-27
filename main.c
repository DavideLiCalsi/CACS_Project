//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SplitSyndrome.h"
#include "Utilities/dataReader.h"

int main(){

    char path[100]="./Utilities/info.txt";
    Info* info=readData(path);

    BinMatrix* A = transpose(*info->H_t);
    BinMatrix* I = identityMatrix((info->n)/2);

    BinMatrix* H = concat(*I,*A,0);
    destroyMatrix(A);
    destroyMatrix(I);

    BinMatrix* e=NULL;
    SplitSyndrome(*H,*info->s,8,&e);

    
    return 0;
}