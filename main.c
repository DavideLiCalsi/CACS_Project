//#include "Matrix.h"
#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SplitSyndrome.h"
#include "Utilities/dataReader.h"
#include "Utilities/utilities.h"
#include "Supercode/Supercode.h"

int main(){

    srand(time(NULL));
    int seed = 24;   // fixed for repeatability
    char path[100] = "./Utilities/info.txt";
    Info *info = readData(path);

    BinMatrix* A = transpose(*info->H_t);
    BinMatrix* I = identityMatrix((info->n)/2);

    BinMatrix* H = concat(*I,*A,0);
    BinMatrix* H_t = transpose(*H);
    BinMatrix* G = concat(*info->H_t,*I,0);

    //printMatrix(*H);
    printMatrix(*G);
    destroyMatrix(A);
    destroyMatrix(I);

    BinMatrix *codeword = generateCodeword(G, seed);
    //printMatrix(*codeword);
    BinMatrix *error = generateError(info->n, info->w, seed);
    //printMatrix(*error);
    BinMatrix *receivedCodeword = vectorSum(*codeword, *error);
    //printMatrix(*receivedCodeword);

    BinMatrix *syndrome = product(*receivedCodeword, *H_t);
    //printMatrix(*syndrome);

    VectorList l = NULL;
    VectorList r = NULL;

    SplitSyndrome(*H,*syndrome,info->w,&l,&r);

    /*puts("INIT---l");
    VectorList_print(l);
    puts("INIT----r");
    VectorList_print(r);
    BinMatrix *con = concat(*(vectorList_pop(&l)->v),*(vectorList_pop(&r)->v),0);
    puts("INIT----CONCAT");
    printMatrix(*con);
    puts("SYNDROME");
    printMatrix(*product(*con, *H_t));
    scanf("%d", &seed);*/

    //int b_vect[20]={ 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1};
    //BinMatrix* b= buildMatrix(b_vect,20,1);
    /*BinMatrix err=*transpose(*zeroVector(info->n));
    putElement(&err,0,0,1);
    putElement(&err,0,3,1);*/

    precomputeBinCoefficients(info->n,(info->n)/2);

    int e=4,y=4;
    BinMatrix *decoded = SupercodeDecoding(*G,*H,*receivedCodeword,info->n,(info->n)/2 ,e,y,info->w);

    printf("\nDist(received codeword,guess): %d\n",HammingDistance(*decoded,*receivedCodeword) );
    printf("Dist(original codeword,received codeword): %d\n\n",HammingDistance(*receivedCodeword,*codeword) );

    printf("Supercode Decoding guess: ");
    printMatrix(*decoded);
    printf("Original codeword: ");
    printMatrix(*codeword);
    printf("Error: ");
    printMatrix(*error);
    printf("Received codeword: ");
    printMatrix(*receivedCodeword);

    printf("Syndrome of original codeword ");
    printMatrix(*product(*codeword, *H_t) );

    printf("Syndrome of decoded codeword ");
    printMatrix(*product(*decoded,*H_t) );

    printf("Syndrome of received codeword ");
    printMatrix(*product(*receivedCodeword, *H_t) );

    //printf("%d--%d--%d", HammingWeight(*codeword), HammingWeight(*error), HammingWeight(*receivedCodeword));

    destroyInfo(info);
    destroyMatrix(H);
    destroyMatrix(H_t);
    destroyMatrix(G);
    destroyMatrix(codeword);
    destroyMatrix(error);
    destroyMatrix(receivedCodeword);
    destroyMatrix(syndrome);
    destroyMatrix(decoded);
    VectorList_destroy(&l);
    VectorList_destroy(&r);
    
    return 0;
}