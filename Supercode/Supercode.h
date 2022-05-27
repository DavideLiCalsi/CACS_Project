#include "../SplitSyndrome/SplitSyndrome.h"
#include "../Utilities/randomSelector.h"
#include "../Utilities/utilities.h"

void SupercodeDecoding(BinMatrix H, BinMatrix b, int n,int k, int e,int y){

    // Compute the syndrome
    BinMatrix* synd = product(H,b);

    // Build the set {1,...,n}
    int i;
    int* set_n_array = (int*) malloc(sizeof(int)*n);

    for (i=0;i<n;++i)
        set_n_array[i]=i;

    Set* n_set = buildSet(set_n_array,n);
    Set* gamma=NULL;

    // Iterate Ln(k,e) times
    for (i=0;i<computeLn(n,k,e);++i){

        // generate the information set
        gamma=getRandomElements(n_set,k);

        BinMatrix* I_n_k=identityMatrix(n-k);
        BinMatrix* A = sampleFromMatrix(n_set->data,n_set->length,H,MATRIX_SAMPLE_COLUMNS);
        BinMatrix* H_prime = concat(*A,*I_n_k,0);
        destroyMatrix(I_n_k);
        destroyMatrix(A);
    }

    int s = (n-k)/y;
    int j;

    // Iterate s times
    for (j=0;j<s;++j){

        BinMatrix* Iy=identityMatrix(y);
        
    }

}