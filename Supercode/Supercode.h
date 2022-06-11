#include "../SplitSyndrome/SupercodeSplitSyndrome.h"
#include "../Utilities/randomSelector.h"
#include "../Utilities/utilities.h"
#include "../Utilities/debug.h"
#include <stdbool.h>

#define DELTA_0 0.11002786443835955

clock_t begin,end;

/**
 * @brief Computes the e2 parameters
 * 
 * @param n 
 * @param k 
 * @param y 
 * @param e1 
 * @param b 
 * @param delta_0 
 * @return int 
 */
int compute_e2(int n, int k, int y, int e1, int b, double delta_0){

    double alpha = e1*1.0/n;
    double v=y*1.0/n;
    double R=k*1.0/n;
    return ceil( (delta_0-alpha)*1.0/(1-R-(b-1)*v) );
}


/**
 * @brief Searches a vector in some lists and counts
 * the occurances of the vector in the lists.
 * This function is called to build the list K(gamma)
 * of vectors appearing in at least b lists K_i
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
            /*puts("GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
            printMatrix(vector);
            printMatrix(*(temp->v));
            printf("%d\n", compareVectors(vector,*(temp->v)));
            puts("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN");*/

            // If you found the target, update counter
            if (compareVectors(vector,*(temp->v)) == 0)
                *counter = *counter + 1;

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
    VectorList examined=NULL;


    for (i=0;i<len;++i){

        
        // Retrieve K_i
        //printf("Working on K(%d)\n",i);
        VectorList vector;
        int index = 0;
        
        // iterate over all the vectors in K_i
        while (1)
        {   
            counter=0;
            vector=vectorList_get(&lists[i], index);
            index++;

            if (vector==NULL)
                break;

            // Search the vector in all the remaining lists
            if (!VectorList_search(*vector->v,examined))
                searchInLists(lists,*vector->v,len-i,&counter);
            
            // Store the vector in the list of the vectors that you already searched
            VectorList_addHead(vector->v,&examined);

            // If it appeared in more than l lists, add it to KGamma
            if (counter>=l)
                VectorList_addHead(vector->v,KGamma);

            /* 
                Free only the node of the vector list, but do not destroy
                the BinMatrix in it because you might need it later
            */
            free(vector);
        }
        /*puts("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        VectorList_print(*KGamma);
        puts("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ");*/
    }

    //VectorList_destroy(&examined);
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
 * on gamma equals m. If you improved the best Hamming distance
 * from the codeword to decode, update the distance and the guess
 * 
 * @param gamma The information set
 * @param m The ref vector
 * @param b The vector to decode
 * @param guess The pointer to the pointer to the current best guess
 * @param guess_dist The best distance so far
 */
void loopOverProjectionOnGamma(Set* gamma,BinMatrix m,BinMatrix H,BinMatrix b,BinMatrix** guess, int* guess_dist){

    int* indexes=gamma->data;
    int gamma_len=gamma->length;
    int n=b.cols;
    BinMatrix *H_t = transpose(H);
    BinMatrix* H_t_gamma= sampleFromMatrix(gamma->data,n/2,*H_t,MATRIX_SAMPLE_ROWS);

    int* complement_indexes=malloc(sizeof(int)*n/2);

    for (int i=0, gamma_index=0, compl_index=0; i<n;++i ){
        
        if (i < gamma->data[gamma_index] || gamma_index==n/2 ){
            complement_indexes[compl_index]=i;
            compl_index++;
        }
        else
            gamma_index++;
    }

    int effective_dimension=0;

    for(int i=0;i<n/2;++i){
        if (complement_indexes[i]<n/2)
            effective_dimension++;
    }
    Set* gamma_compl=buildSet(complement_indexes,n/2);

    BinMatrix* H_t_complement =sampleFromMatrix(complement_indexes,n/2,*H_t,MATRIX_SAMPLE_ROWS);

    effective_dimension=n/2-effective_dimension;

    int* effective_indexes=(int*) malloc(sizeof(int)*effective_dimension);

    for(int i=0;i<effective_dimension;++i)
        effective_indexes[i]=n/2-(effective_dimension-i);

    BinMatrix* H_t_complement_eff = sampleFromMatrix(effective_indexes,effective_dimension,*H_t_complement,MATRIX_SAMPLE_ROWS);
    
    /*Compute the syndrome of vector m w.r.t. gamma*/
    BinMatrix* m_syndrome_gamma=product(m,*H_t_gamma);

    BinMatrix* full_vector;

    //We will try all possible values for this vector,starting from [0,...,0]
    //Its entries represent the non-fixed entries of the new vector
    BinMatrix* curr_trial=zeroVector(effective_dimension);

    // We will stop when we reach final, i.e. a vector full of ones
    BinMatrix* final = oneVector(effective_dimension);

    int i;
    int gamma_index, curr_trial_index;
    // TODO: add check that the generated vector is a valid codeword
    
    int curr_as_int=0;
    int final_int= 1 << effective_dimension;
    while (compareVectors(*curr_trial,*final)!=0)
    {   

        // Initialize the full vector to all zeros
        full_vector=zeroVector(n);

        for(int i=0, gamma_index=0, compl_index=0,effective_index=0;i<n;++i){
            
            // The current position is in gamma, copy from m
            if (gamma->data[gamma_index]==i){
                putElement(full_vector,0,i,getElement(m,0,gamma_index));
                gamma_index++;
                continue;
            }

            // Curr position not in gamma
            if (gamma_compl->data[compl_index]==i ){
                
                // This position is not multiplied by the identity, copy from the current trial
                if (i>=n/2 ){
                    //printf("EffIndex %d, I index %d, elem %d\n",effective_index,i,getElement(*curr_trial,0,effective_index));
                    putElement(full_vector,0,i,getElement(*curr_trial,0,effective_index));
                    effective_index++;
                }
                else{
                    putElement(full_vector,0,i,0);
                }

                compl_index++;
            }
        }

        BinMatrix* target_syndrome=product(*full_vector,*H_t);

        bool failed=false;

        for(int i=0, compl_index=0;i<target_syndrome->cols;++i){
            
            if ( getElement(*target_syndrome,0,i)==1 ){

                if (  i != complement_indexes[compl_index ] ) {
                    failed=true;
                    break;
                }
                else{
                    //printf("Correcting %d\n",i);
                    char curr_elem = getElement(*full_vector,0,i);
                    putElement(full_vector,0,i,curr_elem==0?1:0);
                    ++compl_index;
                }
                    
            }
            else{
                if (  i == complement_indexes[compl_index ] ) 
                    compl_index++; 
            }
                
        }
        /*printMatrix(*curr_trial);*/
        /*BinMatrix* syndrome = product(*full_vector, *H_t);
        printMatrix(*full_vector);
        printMatrix(*target_syndrome);
        printMatrix(*syndrome);
        printSet(gamma_compl);
        puts("end");*/

        if (failed){
            
            destroyMatrix(target_syndrome);
            destroyMatrix(full_vector);
            incrementVector(curr_trial);
            continue;
        }

        /*BinMatrix* zero = zeroVector(gamma_len);
        BinMatrix* syndrome = product(*full_vector, *H_t);
        //printMatrix(*syndrome);

        if (compareVectors(*syndrome,*zero)!=0){
            puts("Error: codeword does not belong to code");
            /*printSet(gamma);
            printSet(gamma_compl);
            printMatrix(*full_vector);
            printMatrix(*target_syndrome);
            printMatrix(*syndrome);
            exit(1);
            incrementVector(curr_trial);
            continue;
        }*/

        // Now check if you improved the distance
        int new_dist = HammingDistance(b,*full_vector);

        PRINTMATRIX(*full_vector);

        // If you improved the distance, update your guess
        if (new_dist<*guess_dist){
            printf("\nImproved: new distance %d\nnew_vector: ",new_dist);
            printMatrix(*full_vector);
            printf("\n");
            BinMatrix* old=*guess;
            *guess=full_vector;
            *guess_dist=new_dist;
            destroyMatrix(old);
        }
        else
            destroyMatrix(full_vector);

        destroyMatrix(target_syndrome);
        /*Code that updates curr_trial*/
        incrementVector(curr_trial);       
    }
    
    destroyMatrix(H_t);
    destroyMatrix(final);
    destroyMatrix(curr_trial);
    free(complement_indexes);
    free(effective_indexes);

}

/**
 * @brief Get an information set for matrix H
 * 
 * @param H Parity check matrix
 * @param n_set The set {1,...,n}
 * @param k Lenght of the information set
 * @param iteration iteration of the algorithm (used to add randomicity to the subset choice)
 * @return Set* The information set
 */
Set* getInformationSet(BinMatrix H,Set* n_set,int k, int iteration){

    Set* gamma;
    BinMatrix* H_sub;

    int i=0;

    //puts("Searching for information set...");
    do{
        gamma=getRandomElements(n_set,k,clock()*(i+iteration+7) + time(NULL)*(iteration+i+3));
        H_sub=sampleFromMatrix(gamma->data,k,H,MATRIX_SAMPLE_COLUMNS);

        //printMatrix(*H_sub);
        //int det = determinant(*H_sub);

        //printf("Determinant %d\n",det);

        
            //printf("Stopped at %d\n",i);
        quicksort(gamma->data,0,k-1);
        destroyMatrix(H_sub);
        return gamma;
        
    }while (true);
    
}

BinMatrix *SupercodeDecoding(BinMatrix G, BinMatrix H, BinMatrix b, int n, int k, int e, int y, int b_param, int d){

    BinMatrix *H_t = transpose(H);

    // The decoded vector, to return at the end
    BinMatrix* decoded=zeroVector(n);

    // The best distance so far
    int best_dist=HammingDistance(b,*decoded);

    int e2=1;compute_e2(n,k,y,e,b_param,DELTA_0);
    int tot_iterations=computeLn(n,k,e);

    printf("--------------------------\n"
           "Best distance starts at %d\n",best_dist);
    printf("Parameters:\n");
    printf("\t-y:%d\n",y);
    printf("\t-e1:%d\n",e);
    printf("\t-e2:%d\n",e2);
    printf("\t-b:%d\n",b_param);
    printf("Total iterations:%d\n",tot_iterations);

    // Compute the syndrome
    BinMatrix* synd = product(b,*H_t);

    printf("Syndrome is: \n");
    printMatrix(*synd);

    // Build the set {1,...,n}
    int i;
    int* set_n_array = (int*) malloc(sizeof(int)*n);

    for (i=0;i<n;++i)
        set_n_array[i]=i;

    Set* n_set = buildSet(set_n_array,n);
    Set* gamma=NULL;

    // Iterate Ln(k,e) times
    for (i=0;i<tot_iterations;++i){
        
        printf("\r--------------------------ITERATION: %d", i);
        fflush(stdout);
        // generate the information set

        begin=clock();
        gamma=getInformationSet(H,n_set,k, i);
        end=clock();

        //printf("Generate subset: %ld\n",end-begin);
        PRINTF("Generated subset gamma: ");
        PRINTSET(gamma);

        // Find H'= [A|I_n-k]
        int* gamma_indexes=gamma->data;

        // Find the set complementary to gamma, i.e., {"whole set"}\{gamma}
        int *compl_gamma_indexes = (int *)malloc(sizeof(int)*(n-k));
        int curr_gamma_index = 0;
        for (int i=0; i<n-k; i++){
            for (int j=0; j<n;){
                if (curr_gamma_index == k){
                    compl_gamma_indexes[i]=j;
                    i++; j++;
                }
                else if ((curr_gamma_index < k) && (j < gamma->data[curr_gamma_index])){
                    compl_gamma_indexes[i]=j;
                    i++; j++;
                }else{
                    curr_gamma_index++;
                    j++;
                }
            }
        }


        BinMatrix* I_n_k=identityMatrix(n-k);
        BinMatrix* A = sampleFromMatrix(gamma_indexes,n-k,H,MATRIX_SAMPLE_COLUMNS);
        BinMatrix* H_prime = concat(*A,*I_n_k,0);

        int s = ceil( (n-k)*1.0/y );
        int* indexes = (int*) malloc(sizeof(int)*y);

        // An array of VectorLists that will contain the results of SplitSyndrome
        VectorList* K = (VectorList*)malloc(sizeof(VectorList)*s);

        // Iterate s times
        for (int j=0;j<s;++j){

            PRINTF("\nMICRO-ITERATION %d/%d\n",j+1,s);

            // Corrections needed to check if we are splitting for the last time
            bool cond = j==s-1 && (n-k)%y != 0;
            int samplingLen = cond ? (n-k)%y : y;

            BinMatrix* Iy=identityMatrix(samplingLen);
            
            for(int z=0; z<samplingLen; ++z){
                indexes[z]=z+samplingLen*j;
            }

            // Build Hi = [Ai|Iy]
            BinMatrix* A_i = sampleFromMatrix(indexes, samplingLen,*A,MATRIX_SAMPLE_ROWS);
            BinMatrix* H_i=concat(*A_i,*Iy,0);
            PRINTMATRIX_PTR(H_i);
            destroyMatrix(A_i);
            destroyMatrix(Iy);

            // Build s_i
            BinMatrix* s_i = sampleFromMatrix(indexes,samplingLen,*synd,MATRIX_SAMPLE_COLUMNS);

            VectorList u_1=NULL;
            VectorList u_2=NULL;

            PRINTF("length of H_i: %d\n",H_i->cols);
            SplitSyndrome(*H_i,*s_i,e+e2, &u_1,&u_2,e,e2,k,samplingLen);
            
            K[j]=u_1;

        }

        // The list of vectors KGamma
        VectorList KGamma=NULL;
        PRINTF("Building K(gamma), processind s=%d lists\n",s);

        begin=clock();
        buildKGamma(K,s,b_param,&KGamma);
        end=clock();

        PRINTF("K(gamma) built in %ld ticks\n",end-begin);

        //VectorList_print(KGamma);
        if (KGamma == NULL)
            PRINTF("WARNING: null K(gamma) iteration %d\n", i);


        // Inspect the list KGamma
        VectorList temp=KGamma;
        BinMatrix* b_Gamma=sampleFromMatrix(gamma->data,gamma->length,b,MATRIX_SAMPLE_COLUMNS);

        while (temp != NULL)
        {
            BinMatrix* ul = temp->v;
        
            BinMatrix* m = vectorSum(*b_Gamma,*ul);

            // Iterate over the vectors c' whose projection on Gamma equals m

            begin=clock();
            loopOverProjectionOnGamma(gamma,*m,H,b,&decoded,&best_dist);
            end=clock();
            PRINTF("Looping over Projection on KGamma took %ld ticks\n",end-begin);

            // Cleaning
            destroyMatrix(ul);
            destroyMatrix(m);

            temp=temp->next;
        }
        
        destroyMatrix(b_Gamma);
        free(K);

    }

    return decoded;

}