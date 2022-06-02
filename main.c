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

    int b_vect[20]={ 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1};
    BinMatrix* b= buildMatrix(b_vect,20,1);
    /*BinMatrix err=*transpose(*zeroVector(info->n));
    putElement(&err,0,0,1);
    putElement(&err,0,3,1);*/

    int e=4,y=10;
    SupercodeDecoding(*H,*b,info->n,(info->n)/2 ,e,y,info->w);
    
    return 0;
}