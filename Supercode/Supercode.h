#include "../SplitSyndrome/SplitSyndrome.h"
#include "../Utilities/randomSelector.h"
#include "../Utilities/utilities.h"

VectorList examined=NULL;

/**
 * @brief Searches a vector in some lists and counts
 * the occurances of the vector in the lsts.
 * 
 * @param lists 
 * @param vector
 * @param len 
 * @param counter The counter that needs to be updated
 */
void searchInLists(VectorList* lists, BinMatrix vector, int len,int* counter){

    int i;

    for(i=0;i<len;++i){

        VectorList list=lists[i];

        VectorList temp=list;

        while (temp != NULL)
        {   
            // If you found the target, update counter
            if (compareMatrices(vector,*(temp->v))==0){

                *counter = *counter +1;
            }

            temp=temp->next;
        }
        

    }
}

/**
 * @brief Given an array of vectorlists Ki, returns
 * the list of the vectors that appear in at least l
 * lists
 * 
 * @param lists 
 * @param len 
 * @return VectorList 
 */
void buildKGamma(VectorList* lists, int len, int l, VectorList* KGamma){

    int i,counter;
    VectorList list;

    for (i=0;i<len;++i){

        list=lists[i];
        VectorList vector=list;
        counter=0;

        while (vector != NULL)
        {
            vector=VectorList_pop(&list);

            if (!VectorList_search(*vector->v,examined))
                searchInLists(lists+i,*vector->v,len-i,&counter);
            
            VectorList_addHead(vector->v,&examined);

            if (counter>=l)
                VectorList_addHead(vector->v,KGamma);
        }
        



    }

}

void SupercodeDecoding(BinMatrix H, BinMatrix b, int n,int k, int e,int y, int d){

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
        gamma=getRandomElements(n_set,k,9);

        // Find H'= [A|I_n-k]
        BinMatrix* I_n_k=identityMatrix(n-k);
        BinMatrix* A = sampleFromMatrix(n_set->data,n_set->length,H,MATRIX_SAMPLE_COLUMNS);
        BinMatrix* H_prime = concat(*A,*I_n_k,0);
        destroyMatrix(I_n_k);
        //destroyMatrix(A);

        int s = (n-k)/y;
        int j;
        int* indexes = (int*) malloc(sizeof(int)*y);

        // Iterate s times
        for (j=0;j<s;++j){

            BinMatrix* Iy=identityMatrix(y);
            
            for(int k=0;k<y;++k){
                indexes[k]=k+y*j;
            }

            // Build Hi = [Ai|Iy]
            BinMatrix* A_i = sampleFromMatrix(indexes,y,*A,MATRIX_SAMPLE_COLUMNS);
            BinMatrix* H_i=concat(*A_i,*Iy,0);
            destroyMatrix(A_i);
            destroyMatrix(Iy);

            // Build si
            BinMatrix* si = sampleFromMatrix(indexes,y,*synd,MATRIX_SAMPLE_COLUMNS);

            BinMatrix* e=NULL;
            SplitSyndrome(*H_i,*si,d, &e);
        }

    }

}