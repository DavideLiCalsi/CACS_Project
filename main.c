#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SupercodeSplitSyndrome.h"
#include "Utilities/dataReader.h"
#include "Utilities/utilities.h"
#include "Supercode/Supercode.h"

double run_test(Info* info,BinMatrix* G, BinMatrix* H, BinMatrix* H_t, int e, int y, int b, int iterations){

    int failed=0, success=0, acceptable=0;
    BinMatrix* zero=zeroVector(info->n/2);

    for (int i=0;i<iterations;++i){

        int seed=69;
        BinMatrix *codeword = generateCodeword(G, seed);
        BinMatrix *error = generateError(info->n, info->w, seed);
        BinMatrix *receivedCodeword = vectorSum(*codeword, *error);

        printf("Error weight: %d\n",HammingWeight(*error));
        puts("Original code");
        printMatrix(*codeword);

        BinMatrix *decoded = SupercodeDecoding(*G,*H,*receivedCodeword,info->n,(info->n)/2 ,e,y,b,info->w);

        BinMatrix* syn=product(*decoded,*H_t);

        if ( compareVectors(*syn,*zero) !=0){
            failed++;
            continue;
        }

        if ( compareVectors(*decoded,*codeword) == 0)
            success++;
        else{
            if (HammingDistance(*decoded,*receivedCodeword)<=HammingDistance(*receivedCodeword,*codeword)){
                acceptable++;
                printf("Found dist %d, true dist %d\n",HammingDistance(*decoded,*receivedCodeword),HammingDistance(*receivedCodeword,*codeword));
            }
            else{
                printf("Found dist %d, true dist %d\n",HammingDistance(*decoded,*receivedCodeword),HammingDistance(*receivedCodeword,*codeword));
                failed++;
            }
        }
        // free memory
        destroyMatrix(codeword);
        destroyMatrix(error);
        destroyMatrix(receivedCodeword);
        destroyMatrix(decoded);
        destroyMatrix(syn);                
    }

    destroyMatrix(zero);

    printf("Successful decodings: %d/%d\n",success,iterations);
    printf("Acceptable decodings: %d/%d\n",acceptable,iterations);
    printf("Failed decodings: %d/%d\n",failed,iterations);

    return success*1.0/iterations;
}

int main(){

    srand(time(NULL));
    int seed = rand();   
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
    BinMatrix *error = generateError(info->n, info->w, seed);
    //printMatrix(*error);
    BinMatrix *receivedCodeword = vectorSum(*codeword, *error);
    //printMatrix(*receivedCodeword);

    BinMatrix *syndrome = product(*receivedCodeword, *H_t);
    //printMatrix(*syndrome);

    VectorList l = NULL;
    VectorList r = NULL;

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

    int e=1,y=1,b=21;
    run_test(info,G,H,H_t,e,y,b,1);

    //printf("%d--%d--%d", HammingWeight(*codeword), HammingWeight(*error), HammingWeight(*receivedCodeword));

    destroyInfo(info);
    destroyMatrix(H);
    destroyMatrix(H_t);
    destroyMatrix(G);
    destroyMatrix(syndrome);
    destroyMatrix(codeword);
    destroyMatrix(receivedCodeword);
    destroyMatrix(error);
    VectorList_destroy(&l);
    VectorList_destroy(&r);
    
    return 0;
}