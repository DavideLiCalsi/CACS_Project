//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SplitSyndrome.h"
#include "Utilities/dataReader.h"
#include "Supercode/Supercode.h"

int main(){

    char path[100]="./Utilities/info.txt";
    Info* info=readData(path);

    BinMatrix* A = transpose(*info->H_t);
    BinMatrix* I = identityMatrix((info->n)/2);

    BinMatrix* H = concat(*I,*A,0);
    destroyMatrix(A);
    destroyMatrix(I);

    VectorList l;
    VectorList r;

    SplitSyndrome(*H,*info->s,info->w,&l,&r);
    VectorList_print(l);

    BinMatrix b=*transpose(*zeroVector(info->n));
    int e=4,y=4;
    SupercodeDecoding(*H,b,info->n, 10,e,y,info->w);
    
    return 0;
}