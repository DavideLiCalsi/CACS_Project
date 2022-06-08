#include "../SplitSyndrome/SupercodeSplitSyndrome.h"
#include "../Utilities/randomSelector.h"
#include "../Utilities/utilities.h"
#include <stdbool.h>

int compute_e2(int n, int k, int y, int e1, int b, double delta_0){

    double alpha = e1*1.0/n;
    double v=y*1.0/n;
    double R=k*1.0/n;
    return ceil( (delta_0-alpha)/(1-R-(b-1)*v) );
}


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
        printf("Working on K(%d)\n",i);
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
 * on gamma equals m
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

    for (int i=0; i<n/2;++i ){
        int j=gamma->data[i];
        complement_indexes[i]=j;
    }
    BinMatrix* H_t_complement =sampleFromMatrix(complement_indexes,n/2,*H_t,MATRIX_SAMPLE_ROWS);
    free(complement_indexes);
    
    /*Compute the syndrome of vector m w.r.t. gamma*/
    BinMatrix* m_syndrome_gamma=product(m,*H_t_gamma);

    // quicksort(indexes,0,gamma_len-1); likely to be useless (indeces already sorted)

    BinMatrix* full_vector;

    //We will try all possible values for this vector,starting from [0,...,0]
    //Its entries represent the non-fixed entries of the new vector
    BinMatrix* curr_trial=zeroVector(n-gamma_len);

    // We will stop when we reach final, i.e. a vector full of ones
    BinMatrix* final = oneVector(n-gamma_len);

    int i;
    int gamma_index, curr_trial_index;
    // TODO: add check that the generated vector is a valid codeword
    

    while (compareVectors(*curr_trial,*final)!=0)
    {   

        
        // Initialize the full vector to all zeros
        full_vector=zeroVector(n);

        /*Code that updates curr_trial*/
        incrementVector(curr_trial);

        BinMatrix* compl_syndrome=product(*curr_trial,*H_t_complement);

        if (compareVectors(*compl_syndrome,*m_syndrome_gamma)!=0)
            continue;

        // Merge m and curr_trial into full_vector
        gamma_index=0;
        curr_trial_index=0;
        for (i=0;i<n;++i){

            // This position is in gamma
            if (gamma_index < gamma_len && i==indexes[gamma_index]){

                char entry=getElement(m,0,gamma_index);

                if (entry!=0)
                    putElement(full_vector,0,i,entry);
                
                gamma_index++;
            }
            else{ //This position is not in gamma

                char entry=getElement(*curr_trial,0,curr_trial_index);

                if (entry!=0)
                    putElement(full_vector,0,i,entry);
                
                curr_trial_index++;
            }
        }

        /*BinMatrix* syndrome = product(*full_vector, *H_t);
        BinMatrix* zero = zeroVector(gamma_len);
        //printMatrix(*syndrome);

        if (compareVectors(*syndrome,*zero)!=0){
            destroyMatrix(syndrome);
            destroyMatrix(zero);
            destroyMatrix(full_vector);
            continue;
        }*/

        // Now check if you improved the distance
        int new_dist = HammingDistance(b,*full_vector);

        printf("NEW DIST: %d\n",new_dist);
        printMatrix(*full_vector);

        // If you improved the distance, update your guess
        if (new_dist<*guess_dist){
            printf("Improved: new dist %d\nnew_vector: ",new_dist);
            printMatrix(*full_vector);
            BinMatrix* old=*guess;
            *guess=full_vector;
            *guess_dist=new_dist;
            destroyMatrix(old);
        }
        else
            destroyMatrix(full_vector);
    }
    
    destroyMatrix(H_t);
    destroyMatrix(final);
    destroyMatrix(curr_trial);

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

    puts("Searching for information set...");
    do{
        gamma=getRandomElements(n_set,k,clock()*(i+iteration+7) + time(NULL)*(iteration+i+3));
        H_sub=sampleFromMatrix(gamma->data,k,H,MATRIX_SAMPLE_COLUMNS);

        //printMatrix(*H_sub);
        int det = determinant(*H_sub);

        //printf("Determinant %d\n",det);

        if (1 || det==1 ){
            //printf("Stopped at %d\n",i);
            quicksort(gamma->data,0,k-1);
            destroyMatrix(H_sub);
            return gamma;
        }

        destroyMatrix(H_sub);
        destroySet(gamma);
        i++;
    }while (true);
    
}

BinMatrix *SupercodeDecoding(BinMatrix G, BinMatrix H, BinMatrix b, int n, int k, int e, int y, int b_param, int d){

    BinMatrix *H_t = transpose(H);

    // The decoded vector, to return at the end
    BinMatrix* decoded=zeroVector(n);

    // The best distance so far
    int best_dist=HammingDistance(b,*decoded);

    int e2=compute_e2(n,k,y,e,b_param,0.11002786443835955);

    printf("--------------------------\n"
           "Best distance starts at %d\n",best_dist);

    // Compute the syndrome
    BinMatrix* synd = product(b,*H_t);

    printf("Found syndrome: \n");
    printMatrix(*synd);

    // Build the set {1,...,n}
    int i;
    int* set_n_array = (int*) malloc(sizeof(int)*n);

    for (i=0;i<n;++i)
        set_n_array[i]=i;

    Set* n_set = buildSet(set_n_array,n);
    Set* gamma=NULL;

    // Iterate Ln(k,e) times
    for (i=0;i<computeLn(n,k,e);++i){

        printf("--------------------------\nITERATION: %d\n", i);
        // generate the information set
        gamma=getInformationSet(H,n_set,k, i);
        printf("Generated subset gamma: ");
        printSet(gamma);

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
        printf("Generated complementary subset of gamma:   {");
        for (int j=0; j<n-k-1; j++)
            printf("%d, ", compl_gamma_indexes[j]);
        printf("%d}\n", compl_gamma_indexes[n-k-1]);

        BinMatrix* I_n_k=identityMatrix(n-k);
        BinMatrix* A = sampleFromMatrix(gamma_indexes,n-k,H,MATRIX_SAMPLE_COLUMNS);
        //BinMatrix* A_compl = sampleFromMatrix(compl_gamma_indexes,n-k,H,MATRIX_SAMPLE_COLUMNS);
        BinMatrix* H_prime = concat(*A,*I_n_k,0);
        //H_prime = concat(*A, *A_compl, 0);

        //printf("Found H'\n");
        //printMatrix(*H_prime);

        int s = ceil( (n-k)*1.0/y );
        int* indexes = (int*) malloc(sizeof(int)*y);

        // An array of VectorLists that will contain the results of SplitSyndrome
        VectorList* K = (VectorList*)malloc(sizeof(VectorList)*s);

        // Iterate s times
        for (int j=0;j<s;++j){

            printf("\nITERATION %d/%d\n",j+1,s);

            // We are splitting for the last time
            bool cond = j==s-1 && (n-k)%y != 0;
            int samplingLen = cond ? (n-k)%y : y;

            BinMatrix* Iy=identityMatrix(samplingLen);
            
            for(int z=0; z<samplingLen; ++z){
                indexes[z]=z+samplingLen*j;
            }

            // Build Hi = [Ai|Iy]
            BinMatrix* A_i = sampleFromMatrix(indexes, samplingLen,*A,MATRIX_SAMPLE_ROWS);
            BinMatrix* H_i=concat(*A_i,*Iy,0);
            destroyMatrix(A_i);
            destroyMatrix(Iy);

            // Build s_i
            BinMatrix* s_i = sampleFromMatrix(indexes,samplingLen,*synd,MATRIX_SAMPLE_COLUMNS);
            //BinMatrix* s_i_transp = sampleFromMatrix(indexes,y,*synd,MATRIX_SAMPLE_COLUMNS);
            //BinMatrix* s_i = transpose(*s_i_transp);
            //destroyMatrix(s_i_transp);

            VectorList u_1=NULL;
            VectorList u_2=NULL;

            printf("Applying SplitSyndrome on s_i:\n");
            printMatrix(*s_i);
            printf("and H_i\n");
            printMatrix(*H_i);
            SplitSyndrome(*H_i,*s_i,3, &u_1,&u_2,e,e2);

            /*puts("INIT---U_1");
            VectorList_print(u_1);
            puts("INIT----U_2");
            VectorList_print(u_2);
            BinMatrix *con = concat(*(vectorList_pop(&u_1)->v),*(vectorList_pop(&u_2)->v),0);
            puts("INIT----CONCAT");
            printMatrix(*con);
            puts("SYNDROME");
            printMatrix(*product(*con, *H_t));
            scanf("%d", &k);*/
            
            K[j]=u_1;

        }

        // The list of vectors KGamma
        VectorList KGamma=NULL;
        printf("Building K(gamma), containing %d lists\n",s);
        buildKGamma(K,s,b_param,&KGamma);
        //printf("Gamma successfully built\n");

        puts("SOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOS");
        VectorList_print(KGamma);
        if (KGamma == NULL)printf("WARNING: iteration %d\n", i);
        puts("ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ");

        // Inspect the list KGamma
        VectorList temp=KGamma;
        
        while (temp != NULL)
        {
            BinMatrix* ul = temp->v;
            BinMatrix* ur=zeroVector(gamma->length-ul->cols);
            BinMatrix* u_full=concat(*ul,*ur,0);
            BinMatrix* b_Gamma=sampleFromMatrix(gamma->data,gamma->length,b,MATRIX_SAMPLE_COLUMNS);
            //BinMatrix* b_Gamma_t=transpose(*b_Gamma);
            //printMatrix(*b_Gamma);
            BinMatrix* m = vectorSum(*b_Gamma,*u_full);

            printf("Found m:\n");
            printMatrix(*m);

            // Iterate over the vectors c' whose projection on Gamma equals m
            loopOverProjectionOnGamma(gamma,*m,H,b,&decoded,&best_dist);

            // Cleaning
            destroyMatrix(ul);
            destroyMatrix(b_Gamma);
            destroyMatrix(m);

            temp=temp->next;
        }
        
        free(K);

    }

    /*printf("FINAL GUESS: ");
    printMatrix(*decoded);
    printMatrix(b);*/
    //BinMatrix *pr = product(*decoded,*H_t);
    //printMatrix(*pr);

    return decoded;

}