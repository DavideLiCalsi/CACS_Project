#include "Matrix/BinaryMatrix.h"
#include "SplitSyndrome/SupercodeSplitSyndrome.h"
#include "Utilities/dataReader.h"
#include "Utilities/utilities.h"
#include "Supercode/Supercode.h"

double run_test(Info* info,BinMatrix* G, BinMatrix* H, BinMatrix* H_t, int e, int y, int b, int iterations){

    int failed=0, success=0, acceptable=0;
    BinMatrix* zero=zeroVector(info->n/2);

    for (int i=0; i<iterations; ++i){

        int seed=rand();

        // generate original codeword, error and received codeword
        BinMatrix *codeword = generateCodeword(G, seed);
        BinMatrix *error = generateError(info->n, info->w-2, seed);
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

    printf("Successful decodings: %d/%d\n", success, iterations);
    printf("Acceptable decodings: %d/%d\n", acceptable, iterations);
    printf("Failed decodings: %d/%d\n", failed, iterations);

    return success*1.0/iterations;
}

int main(){

    // set random seed
    srand(time(NULL));
    int seed = rand();  

    // get general information: n, seed, w (minimum distance), H_t (without identity matrix)
    char path[100] = "./Utilities/info.txt";
    Info *info = readData(path);

    // compute actual parity matrix H and related generator matrix G
    BinMatrix* A = transpose(*info->H_t);
    BinMatrix* I = identityMatrix((info->n)/2);

    BinMatrix* H = concat(*I,*A,0);
    BinMatrix* H_t = transpose(*H);
    BinMatrix* G = concat(*info->H_t, *I,0);

    // compute needed binomial coefficients only once
    precomputeBinCoefficients(info->n,(info->n)/2);

    // set parameters found through "find_param.py"
    int e=1,y=1,b=2;

    // run test "iterations" times
    int iterations = 1;
    run_test(info, G, H, H_t, e, y, b, iterations);


    // free memory
    destroyMatrix(A);
    destroyMatrix(I);
    destroyInfo(info);
    destroyMatrix(H);
    destroyMatrix(H_t);
    destroyMatrix(G);
    
    return 0;
}