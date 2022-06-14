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
 * @param lists Array of the lists K1,K2,...,Kn
 * @param len The number of lists
 * @param l In how many lists the target vectors should appear
 * @param KGamma At the end of the algorithm it will contain all the
 * desired vectors
 */
void buildKGamma(VectorList* lists, int len, int l, VectorList* KGamma){

    int i,counter;
    VectorList examined=NULL;


    for (i=0;i<len;++i){

        
        // Retrieve K_i
        VectorList vector;
        
        // iterate over all the vectors in K_i
        while (1)
        {   
            counter=0;
            vector=vectorList_pop(&lists[i]);

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

            VectorList_destroy(&vector);
        }
    }

    VectorList_destroy(&examined);
}

/**
 * @brief Interprets the row vector v as a binary number and increments it
 * 
 * @param v Ptr to the vector to increment
 */
void incrementVector(BinMatrix* v){

    char k=1;
    int len= isRowVector( (*v) )? v->cols:v->rows;
    int i;

    if ( isRowVector( (*v) ) ){
        for (i=len-1;i>=0;--i){
            char curr=getElement(*v,0,i);
            char new_digit=curr+k;
            k=new_digit/2;
            new_digit=new_digit%2;
            putElement(v,0,i,new_digit);
        }
        return;
    }

    if (isColumnVector( (*v) )){
        for (i=len-1;i>=0;--i){
            char curr=getElement(*v,i,0);
            char new_digit=curr+k;
            k=new_digit/2;
            new_digit=new_digit%2;
            putElement(v,i,0,new_digit);
        }
    }
    

}

/**
 * @brief This function takes a vector m and generates all the codewords
 * belonging to the target code whos eprojection on the set Gamma equals m.
 * For each of these codewords, it checks if their distance from the codeword to
 * decode is less than the previous best guess, and possibly updates the new best guess.
 * 
 * @param G Generator matrix
 * @param H Parity-check matrix
 * @param m The vector to match
 * @param b Vector to decode
 * @param gamma Set of k positions to consider
 * @param guess The current best guess for the decoded codeword
 * @param guess_dist The current best distance
 */
void generateAllProjections(BinMatrix G, BinMatrix H, BinMatrix m, BinMatrix b, Set* gamma,BinMatrix** guess, int* guess_dist){
    int n=G.cols;
    int k=G.rows;

    BinMatrix* old;

    /* 
    IDEA: We'll apply Gauss-Jornad elimination to solve the system G(gamma)*x=m^t.
    Once we have found all the column vectors x that satisfy this, we can simply find 
    all the codewods with the dsired propery by computing x^t * G.
    NOTE: x^t means the transpose of x.
    */

    // First sample from G using indexes in gamma to find G(gamma)
    BinMatrix* G_gamma=(sampleFromMatrix(gamma->data,k,G,MATRIX_SAMPLE_COLUMNS));
    old=G_gamma;
    G_gamma=transpose(*G_gamma);
    destroyMatrix(old);
    
    // Construct the augmented matrix A=[G(gamma)|m]
    BinMatrix *m_transpose = transpose(m);
    BinMatrix *Gelim_input=concat(*G_gamma,*m_transpose,0);
    destroyMatrix(m_transpose);

    // Apply elimination
    BinMatrix* reduced=GaussElimination(*Gelim_input);

    //Let's check the effective number of vectors that we'll need to try
    int eff_dim=0,i=0,j=0;

    // At the end of this loop, eff_dim will contain the rank
    // of the G(gamma) matrix
    do{

        if (getElement(*reduced,i,i)!=1)
            break;
        
        eff_dim++;
        i++;

    }while (i<k);

    // On our machines, trying more than 2^18 becomes impractical.
    if ((k-eff_dim)>=18) 
        return;

    /*
    If this condition holds, then the original matix was invertible. Thus
    you can simply invert gamma G(gamma) and find the only solution x.
    */
    if (eff_dim==k){
        
        begin=clock();
        // Final is the solution, a.k.a. x
        BinMatrix *G_gamma_inv = inverse(*G_gamma);
        BinMatrix *m_transpose = transpose(m);
        BinMatrix* final= product(*G_gamma_inv,*m_transpose);
        end=clock();
        PRINTF("Inversion took %ld\n",end-begin);
        destroyMatrix(G_gamma_inv);
        destroyMatrix(m_transpose);
        
        // Transpose x and compute the new guess
        BinMatrix *final_transpose = transpose(*final);
        BinMatrix *new_guess=product(*final_transpose,G);
        destroyMatrix(final_transpose);

        // A double-check to ensure that the math is actually correct. Can be removed later.
        BinMatrix *test=sampleFromMatrix(gamma->data,k,*new_guess,MATRIX_SAMPLE_COLUMNS);

        if ( compareVectors(*test,m)!=0){
            puts("ALARM");
            printMatrix(*transpose(*final));
            printMatrix(m);
            printMatrix(*test);
            exit(-1);
        }
        

        int new_dist = HammingDistance(b,*new_guess);

        /* 
        If you improved the distance, update your guess.
        Since this is the only guess you can have here, cleanup and return
        */
        if (new_dist<*guess_dist){
            printf("\nImproved: new distance %d\nnew_vector: ",new_dist);
            printMatrix(*new_guess);
            printf("\n");
            BinMatrix* old=*guess;
            *guess=new_guess;
            *guess_dist=new_dist;

            destroyMatrix(old);
            destroyMatrix(final);
            destroyMatrix(test);
            destroyMatrix(reduced);
            destroyMatrix(G_gamma);
            destroyMatrix(Gelim_input);
            return;
        }
        else{
            destroyMatrix(new_guess);
            destroyMatrix(final);
            destroyMatrix(test);
            destroyMatrix(reduced);
            destroyMatrix(G_gamma);
            destroyMatrix(Gelim_input);
            return;
        }
    }

    /*
     If did not return before, than it means there is 
     more than one result to check. We will need to try
     2^(k-eff_dim) vectors.
    */

    // Initialize curr_trial to k zeros
    BinMatrix *zero = zeroVector(k);
    BinMatrix* curr_trial=transpose(*zero);
    destroyMatrix(zero);

    // The indexes from 0 to eff_dim-1
    int* indexes= (int*) (malloc(sizeof(int)*(eff_dim)));

    // The indexes fro eff_dim to k
    int* non_eff_indexes= (int*) (malloc(sizeof(int)*(k-eff_dim)));

    for (int i=0;i<eff_dim;++i)
        indexes[i]=i;

    for (int i=0;i<k-eff_dim;++i)
        non_eff_indexes[i]=eff_dim+i;
    
    // Invert the matrix corresponding to the first eff_dim columns of the reduced matrix
    BinMatrix* invertible=sampleFromMatrix(indexes,eff_dim,*reduced,MATRIX_SAMPLE_COLUMNS);
    old=invertible;
    invertible=sampleFromMatrix(indexes,eff_dim,*invertible,MATRIX_SAMPLE_ROWS);
    destroyMatrix(old);

    begin=clock();
    BinMatrix* inv=inverse(*invertible);
    end=clock();
    PRINTF("Inversion took: %d\n",end-begin);

    int curr_as_int = 0;
    int final=(1<< (k-eff_dim));

    int last_index = k;

    // Retrieve the target value to find after the GaussJordan elimination
    BinMatrix* new_m = sampleFromMatrix(&last_index,1,*reduced,MATRIX_SAMPLE_COLUMNS);

    // The indexes from 0 to k-1
    int* first_indexs=(int*) malloc(sizeof(int)*k);

    for(int i=0;i<k;++i)
        first_indexs[i]=i;

    // The all-bu-last columns of the reduced matrix, i.e. we exclude the column that represents m
    BinMatrix* reduced_left=sampleFromMatrix(first_indexs,k,*reduced,MATRIX_SAMPLE_COLUMNS);

    while (curr_as_int < final)
    {

        /*
        At this point the linear system has the form (using code variables' names)
            reduced_left*final = new_m

        Curr_trial always has the form 
            [   0...0      | cur_as_int ]
              eff_dim bits    k-eff_dim bits
        
        Thus the product reduced_left*curr_trial is equivalent to multiplying the right part
        of curr_trial and the last k-eff_dim columns of reduced_left.
        We then add this result to new_m.
        */
        BinMatrix* x = product(*reduced_left,*curr_trial);
        BinMatrix* t = vectorSum(*new_m,*x);
        destroyMatrix(x);

        /*
        We sample the last k-eff_dim entries of t. If they are zero, it means
        we are able to generate a proper solution with this trial. Otherwise, just move 
        to the next one.
        */
        BinMatrix *matrix_sampled = sampleFromMatrix(non_eff_indexes,k-eff_dim,*t,MATRIX_SAMPLE_ROWS);
        BinMatrix *last_target = transpose(*matrix_sampled);
        BinMatrix *ref=zeroVector(k-eff_dim);
        destroyMatrix(matrix_sampled);

        if (compareVectors(*last_target,*ref)!=0){
            ++curr_as_int;
            incrementVector(curr_trial);
            destroyMatrix(last_target);
            destroyMatrix(ref);
            destroyMatrix(t);
            continue;
        }

        // Let's consider only the first eff_dim entries of t now
        x=sampleFromMatrix(indexes,eff_dim,*t,MATRIX_SAMPLE_ROWS);

        // Multiplying by inv will return the left part of final, which is the solution of the system we were looking for
        BinMatrix* x_left = product(*inv,*x );
        destroyMatrix(x);
        
        BinMatrix* non_effective= sampleFromMatrix(non_eff_indexes,k-eff_dim,*curr_trial,MATRIX_SAMPLE_ROWS);

        // Concatenating with the last part of curr_trial yields the full solutions

        // Now generate the guess corresponding to full
        BinMatrix *final_transpose = transpose(*final);
        BinMatrix* new_guess=product(*final_transpose,G);
        destroyMatrix(final_transpose);
 
        BinMatrix* test=sampleFromMatrix(gamma->data,k,*new_guess,MATRIX_SAMPLE_COLUMNS);

        int new_dist = HammingDistance(b,*new_guess);

        // If you improved the distance, update your guess
        if (new_dist<*guess_dist){
            printf("\nImproved: new distance %d\nnew_vector: ",new_dist);
            printMatrix(*new_guess);
            printf("\n");
            BinMatrix* old=*guess;
            *guess=new_guess;
            *guess_dist=new_dist;
            destroyMatrix(old);
        }
        else
            destroyMatrix(new_guess);

        incrementVector(curr_trial);
        ++curr_as_int;

        destroyMatrix(t);
        destroyMatrix(x_left);
        destroyMatrix(non_effective);
        destroyMatrix(final);
        destroyMatrix(test);
        destroyMatrix(last_target);
        destroyMatrix(ref);
        //destroyMatrix(x_left);
        
    }


 final_cleanup:
    free(indexes);
    free(non_eff_indexes);
    free(first_indexs);
    destroyMatrix(invertible);
    destroyMatrix(inv);
    destroyMatrix(curr_trial);
    destroyMatrix(reduced);
    destroyMatrix(reduced_left);
    destroyMatrix(new_m);
    destroyMatrix(Gelim_input);
    destroyMatrix(G_gamma);

}


/**
 * @brief Generates a random subset of the set {1,...,n} of length k.
 * In spite of the function's name, this set is not always guaranteed to be an information
 * set.
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

    gamma=getRandomElements(n_set,k,clock()*(iteration+7) + time(NULL)*(iteration+3));
    //H_sub=sampleFromMatrix(gamma->data,k,H,MATRIX_SAMPLE_COLUMNS);

    quicksort(gamma->data,0,k-1);
    //destroyMatrix(H_sub);
    return gamma;
}

BinMatrix *SupercodeDecoding(BinMatrix G, BinMatrix H, BinMatrix b, int n, int k, int e, int y, int b_param, int d){

    BinMatrix *H_t = transpose(H);

    // The decoded vector, to return at the end
    BinMatrix* decoded=zeroVector(n);

    // The best distance so far
    int best_dist=HammingDistance(b,*decoded);

    // Compute the parameter e2 and the number of iterations
    int e2=compute_e2(n,k,y,e,b_param,DELTA_0);
    unsigned long int tot_iterations=computeLn(n,k,e);

    printf("--------------------------\n"
           "Best distance starts at %d\n",best_dist);
    printf("Parameters:\n");
    printf("\t-y:%d\n",y);
    printf("\t-e1:%d\n",e);
    printf("\t-e2:%d\n",e2);
    printf("\t-b:%d\n",b_param);
    printf("Total iterations:%ld\n",tot_iterations);

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
    free (set_n_array);
    Set* gamma=NULL;

    // Iterate Ln(k,e) times
    for (i=0;i<tot_iterations;++i){
        
        printf("\r--------------------------ITERATION: %d", i);
        fflush(stdout);
        
        // generate the information set
        begin=clock();
        gamma=getInformationSet(H,n_set,k, i);
        end=clock();
;
        PRINTF("Generated subset gamma: ");
        PRINTSET(gamma);

        // Find H'= [A|I_n-k]
        int* gamma_indexes=gamma->data;

        BinMatrix* I_n_k = identityMatrix(n-k);
        BinMatrix* A = sampleFromMatrix(gamma_indexes,n-k,H,MATRIX_SAMPLE_COLUMNS);
        BinMatrix* H_prime = concat(*A,*I_n_k,0);

        int s = ceil( (n-k)*1.0/y );
        int* indexes = (int*) malloc(sizeof(int)*y);

        // An array of VectorLists that will contain the results of SplitSyndrome
        VectorList* K = (VectorList*)malloc(sizeof(VectorList)*s);

        // Iterate s times
        for (int j=0;j<s;++j){

            PRINTF("\nMICRO-ITERATION %d/%d\n",j+1,s);

            /* 
            Corrections needed to check if we are splitting for the last time.
            If that's the case and s does not divide n-k, the last split must
            take that into account.
            */
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

            // Build s_i by samplying from the syndrome
            BinMatrix* s_i = sampleFromMatrix(indexes,samplingLen,*synd,MATRIX_SAMPLE_COLUMNS);

            VectorList u_1=NULL;
            VectorList u_2=NULL;

            PRINTF("length of H_i: %d\n",H_i->cols);
            /*
            Apply the SplitSyndrome to find a list of vectors K[i]
            such that K[i]={u_1 | u=(u_1 | u_2) such that H_i * u =s_i}
            */
            SplitSyndrome(*H_i,*s_i,e+e2, &u_1,&u_2,e,e2,k,samplingLen);
            
            K[j]=u_1;

            destroyMatrix(H_i);
            destroyMatrix(s_i);
            VectorList_destroy(&u_2);

        }

        free(indexes);
        destroyMatrix(I_n_k);
        destroyMatrix(A);
        destroyMatrix(H_prime);

        // The list of vectors KGamma
        VectorList KGamma=NULL;
        PRINTF("Building K(gamma), processing s=%d lists\n",s);

        begin=clock();
        /*
        Build the list KGamma of those vectors that appear in at least
        b_param lists K[i]
        */
        buildKGamma(K,s,b_param,&KGamma);
        end=clock();

        PRINTF("K(gamma) built in %ld ticks\n",end-begin);

        if (KGamma == NULL)
            printf("WARNING: null K(gamma) iteration %d\n", i);


        // Inspect the list KGamma
        VectorList temp=KGamma;
        BinMatrix* b_Gamma=sampleFromMatrix(gamma->data,gamma->length,b,MATRIX_SAMPLE_COLUMNS);

        // Iterate over all the vectors ul in KGamma
        while (temp != NULL)
        {
            BinMatrix* ul = temp->v;

            // Get m=b_Gamma+ul
            BinMatrix* m = vectorSum(*b_Gamma,*ul);

            /* 
            Iterate over the vectors c' whose projection on Gamma equals m,
            checl if some of them improve you HammingDistance from the codeword to decode
            */
            begin=clock();
            generateAllProjections(G,H,*m,b,gamma,&decoded,&best_dist);
            end=clock();
            PRINTF("Looping over Projection on Gamma took %ld ticks\n",end-begin);

            // Cleaning
            destroyMatrix(m);

            temp=temp->next;
        }
        
        destroyMatrix(b_Gamma);
        destroySet(gamma);
        VectorList_destroy(&KGamma);
        free(K);
    }

    destroyMatrix(H_t);
    destroyMatrix(synd);
    destroySet(n_set);

    return decoded;

}