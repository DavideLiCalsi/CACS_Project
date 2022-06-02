//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
#include "Utilities/dataReader.h"
#include "Utilities/utilities.h"
#include "SplitSyndrome/SplitSyndrome.h"


int main(){

    int seed = 6;   // fixed for repeatability
    char path[100] = "./Utilities/info.txt";
    Info *info = readData(path);

    BinMatrix* A = transpose(*info->H_t);
    BinMatrix* I = identityMatrix((info->n)/2);

    // H=[I|A]
    BinMatrix* H = concat(*I,*A,0);
    BinMatrix* H_t = transpose(*H);

    //G*(H^t)=0 --> G=[A^t|I]
    BinMatrix* G = concat(*info->H_t,*I,0);

    destroyMatrix(A);
    destroyMatrix(I);

    printMatrix(*G);
    printMatrix(*H);

    /*
    // check that G*(H^t)=0
    BinMatrix* res = product(*G,*transpose(*H));
    printMatrix(*res);
    */

    BinMatrix *codeword = generateCodeword(G, seed);
    printMatrix(*codeword);

    /*
    // check that codeword*(H^t)=0
    BinMatrix* res = product(*codeword,*transpose(*H));
    printMatrix(*res);
    */

    BinMatrix *error = generateError(info->n, info->w, seed);
    printMatrix(*error);

    BinMatrix *receivedCodeword = vectorSum(*codeword, *error);
    printMatrix(*receivedCodeword);

    BinMatrix *syndrome = product(*receivedCodeword, *H_t);

    generateAllCodeword(G);


    /*BinMatrix* e=NULL;

    SplitSyndrome(*H,*syndrome,info->w,&e);

    if (e == NULL)
    printf("FAILURE! ERROR NOT FOUND\n");
    else
    printMatrix(*e);
    printMatrix(*error);*/

    /*
    printMatrix(*product(*vectorSum(*codeword,*error),*H_t));
    printMatrix(*product(*vectorSum(*codeword,*e),*H_t));

    printMatrix(*product(*vectorSum(*receivedCodeword,*error), *H_t));
    printMatrix(*product(*vectorSum(*receivedCodeword,*e), *H_t));*/

    int vector[10];

    //iterateOverM_Vectors(10,3);
    destroyMatrix(G);
    destroyMatrix(H);
    destroyMatrix(H_t);
    destroyMatrix(codeword);
    destroyMatrix(error);
    destroyMatrix(receivedCodeword);
    destroyMatrix(syndrome);
    destroyInfo(info);


    return 0;
}