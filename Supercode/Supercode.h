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
 * @param lists List of the lists K1,K2,...,Kn
 * @param len The number of lists
 * @return VectorList 
 */
void buildKGamma(VectorList* lists, int len, int l, VectorList* KGamma){

    int i,counter;
    VectorList* list;

    for (i=0;i<len;++i){

        
        // Retrieve Ki
        printf("Working on K(%d)\n",i);
        lists[i];
        VectorList vector=lists[i];
        counter=0;
        
        // iterate over all the vectors in Ki
        while (1)
        {
            vector=VectorList_pop(lists+i);

            if (vector==NULL)
                break;

            // Search the vector in all the remaining lists
            if (!VectorList_search(*vector->v,examined))
                searchInLists(lists+i,*vector->v,len-i,&counter);
            
            // Store the vector in the list of the vectors that you already searched
            VectorList_addHead(vector->v,&examined);

            // If it appeared in more than l lists, add it to KGamma
            if (counter>=l)
                VectorList_addHead(vector->v,KGamma);
        }

    }

}

/**
 * @brief Interprets the row vector v as a binary number and increments it
 * 
 * @param v Ptr to the vector to increment
 */
void incrementVector(BinMatrix* v){

    char k=1;
    int len=v->cols;
    int i;

    for (i=len-1;i>=0;--i){
        char curr=getElement(*v,0,i);
        char new_digit=curr+k;
        k=new_digit/2;
        new_digit=new_digit%2;
        putElement(v,0,i,new_digit);
    }

}

/**
 * @brief Loops over all the vectors c' whose projection
 * on gamma equals m
 * 
 * @param gamma The information set
 * @param m The ref vector
 * @param b The vector to decode
 * @param guess The current best guess
 * @param guess_dist The best distance so far
 */
void loopOverProjectionOnGamma(Set* gamma,BinMatrix m,BinMatrix H,BinMatrix b,BinMatrix* guess, int* guess_dist){

    int* indexes=gamma->data;
    int gamma_len=gamma->length;
    int n=b.rows;

    quicksort(indexes,0,gamma_len);

    BinMatrix* full_vector;

    //We will try all possible values for this vector,starting from [0,...,0]
    //Its entries represent the non-fixed entries of the new vector
    BinMatrix* curr_trial=zeroVector(n-gamma_len);

    // We will stop when we reach final, i.e. a vector full of ones
    BinMatrix* final =oneVector(n-gamma_len);

    int i;
    int gamma_index, curr_trial_index;

    // TODO: add check that the generated vector is a valid codeword
    while (compareVectors(*curr_trial,*final)!=0)
    {
        // Initialize the full vector to all zeros
        full_vector=zeroVector(n);

        /*Code that updates curr_trial*/
        incrementVector(curr_trial);

        // Merge m and curr_trial into full_vector
        gamma_index=0;
        curr_trial_index=0;
        for (i=0;i<n;++i){
            
            // This position is in gamma
            if (i==indexes[gamma_index]){
                //puts("In gamma");
                char entry=getElement(m,0,gamma_index);

                if (entry!=0)
                    putElement(full_vector,0,i,entry);
                
                gamma_index++;
            }
            else{ //This position is not in gamma
                //puts("NOT gamma");
                char entry=getElement(*curr_trial,0,curr_trial_index);

                if (entry!=0)
                    putElement(full_vector,0,i,entry);
                
                curr_trial_index++;

            }
        }

        BinMatrix* syndrome = product(H,*full_vector);
        BinMatrix* zero = zeroVector(gamma_len);

        if (compareVectors(*syndrome,*zero)==0)
            continue;

        // Now check if you improved the distance
        int new_dist = HammingDistance(*transpose(b),*full_vector);

        //printf("NEW DIST: %d\n",new_dist);

        // If you improved the distance, update your guess
        if (new_dist<*guess_dist){
            printf("Improved: new dist %d\n",new_dist);
            printMatrix(*full_vector);
            printMatrix(b);
            BinMatrix* old=guess;
            *guess=*full_vector;
            *guess_dist=new_dist;
            destroyMatrix(old);
        }
    }
    


}

/**
 * @brief Get an information set for matrix H
 * 
 * @param H Parity check matrix
 * @param n_set The set {1,...,n}
 * @param k Lenght of the information set
 * @return Set* The information set
 */
Set* getInformationSet(BinMatrix H,Set* n_set,int k){

    Set* gamma;
    BinMatrix* H_sub;

    int i=0;

    puts("Searching for information set...");
    do{
        gamma=getRandomElements(n_set,k,clock());
        H_sub=sampleFromMatrix(gamma->data,k,H,MATRIX_SAMPLE_COLUMNS);

        //printMatrix(*H_sub);
        int det = determinant(*H_sub);

        //printf("Determinant %d\n",det);

        if (det==1 ){
            printf("Stopped at %d\n",i);
            quicksort(gamma->data,0,k-1);
            destroyMatrix(H_sub);
            return gamma;
        }

        destroyMatrix(H_sub);
        destroySet(gamma);
    }while (true);
    
}

void SupercodeDecoding(BinMatrix H, BinMatrix b, int n,int k, int e,int y, int d){

    // The decoded vector, to return at the end
    BinMatrix* decoded=zeroVector(n);

    // The best distance so far
    int best_dist=HammingDistance(*transpose(b),*decoded);

    printf("Best distance starts at %d\n",best_dist);

    // Compute the syndrome
    BinMatrix* synd = product(H,b);

    printf("Found syndrome\n");
    //printMatrix(*synd);

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
        gamma=getInformationSet(H,n_set,k);
        printf("Generated subset gamma: ");
        printSet(gamma);

        // Find H'= [A|I_n-k]
        int* gamma_indexes=gamma->data;//gamma->data;
        BinMatrix* I_n_k=identityMatrix(n-k);
        BinMatrix* A = sampleFromMatrix(gamma_indexes,n-k,H,MATRIX_SAMPLE_COLUMNS);
        BinMatrix* H_prime = concat(*A,*I_n_k,0);
        destroyMatrix(I_n_k);

        printf("Found H'\n");
        printMatrix(*H_prime);
        /*BinMatrix* H_stand=standardizeParityMatrix(&H);
        printMatrix(*H_stand);*/

        int s = ceil( (n-k)*1.0/y );
        int* indexes = (int*) malloc(sizeof(int)*y);

        // The list of vectors KGamma
        VectorList KGamma=NULL;

        // An array of VectorLists that will contain the results of SplitSyndrome
        VectorList* K= (VectorList*) malloc(sizeof(VectorList)*s);

        // Iterate s times
        for (int j=0;j<s;++j){

            // We are splitting for the last time
            if (j==s-1 && (n-k)%y != 0){
                y=(n-k)%y;
            }

            BinMatrix* Iy=identityMatrix(y);
            
            for(int k=0;k<y;++k){
                indexes[k]=k+y*j;
            }

            // Build Hi = [Ai|Iy]
            BinMatrix* A_i = sampleFromMatrix(indexes,y,*A,MATRIX_SAMPLE_ROWS);
            BinMatrix* H_i=concat(*A_i,*Iy,0);
            destroyMatrix(A_i);
            destroyMatrix(Iy);

            // Build s_i
            BinMatrix* s_i_transp = sampleFromMatrix(indexes,y,*synd,MATRIX_SAMPLE_ROWS);
            BinMatrix* s_i = transpose(*s_i_transp);
            destroyMatrix(s_i_transp);

            VectorList u_1=NULL;
            VectorList u_2=NULL;
            SplitSyndrome(*H_i,*s_i,2, &u_1,&u_2);

            if (u_1==NULL)
                puts("U1 null");
            puts("U_1");
            VectorList_print(u_1);
            K[j]=u_1;
        }

        printf("Building K(gamma), containing %d lists\n",s);
        buildKGamma(K,s,2,&KGamma);
        printf("Gamma successfully built\n");

        // Inspect the list KGamma
        VectorList temp=KGamma;

        while (temp != NULL)
        {
            BinMatrix* ul = temp->v;
            BinMatrix* ur=zeroVector(gamma->length-ul->cols);
            BinMatrix* u_full=concat(*ul,*ur,0);
            BinMatrix* b_Gamma=sampleFromMatrix(gamma->data,gamma->length,b,MATRIX_SAMPLE_ROWS);
            BinMatrix* b_Gamma_t=transpose(*b_Gamma);
            printMatrix(*b_Gamma);
            BinMatrix* m = vectorSum(*b_Gamma_t,*u_full);

            printf("Found m:\n");
            printMatrix(*m);

            /*
            Iterate over the vectors c' whose projection on Gamma 
            equals m
            */
           loopOverProjectionOnGamma(gamma,*m,*H_prime,b,decoded,&best_dist);

           // Cleaning
           destroyMatrix(ul);
           destroyMatrix(b_Gamma);
           destroyMatrix(m);

           temp=temp->next;
        }
        

    }

    printf("FINAL GUESS: ");
    printMatrix(*decoded);
    printf("%d\n",HammingDistance(*decoded,*transpose(b)) );

}