#ifndef SPLIT_SYNDROME_H
#define SPLIT_SYNDROME_H

#include "../Tree/BST.h"
#include "PairSet.h"
#include <math.h>
#include <time.h>

unsigned int binCoefficients[100][100];
int valid[100][100];

/**
 * @brief compute the (n,k) binomial coefficient
 * 
 * @param n 
 * @param k 
 * @return int 
 */
unsigned int binomialCoeff(unsigned int n, unsigned int k){
    // Base Cases
    if (k > n){
        valid[n][k]=1;
        binCoefficients[n][k]=0;
        return 0;
    }
    if (k == 0 || k == n)
    {
        valid[n][k]=1;
        binCoefficients[n][k]=1;
        return 1;
    }
    if (k==1)
    {
        valid[n][k]=1;
        binCoefficients[n][k]= n;
        return n;
    }

    if (valid[n][k] != -1)
        return binCoefficients[n][k];
 
    // Recur
    unsigned int t1,t2;

    if (valid[n-1][k-1] == -1){
        t1=binomialCoeff(n - 1, k - 1);
        binCoefficients[n-1][k-1]=t1;
        valid[n-1][k-1]=1;
    }
    else
        t1=binCoefficients[n-1][k-1];

    if (valid[n-1][k] == -1){
        t2=binomialCoeff(n - 1, k);
        binCoefficients[n-1][k]=t2;
        valid[n-1][k]=1;
    }
    else
        t1=binCoefficients[n-1][k];

    return t1+t2;
}

/**
 * @brief Precomputes and stores all the
 * coeeficients of the type binCoeff(n-m,t)
 * for m from 1 to n-1 and t from 1 to w
 * 
 * @param n 
 * @param t 
 */
void precomputeBinCoefficients(unsigned int n, unsigned int w){

    unsigned int m,t;

    for (m = 1; m <= n; m++) binCoefficients[0][m] = 0;
    for (t = 0; t <= w; t++) binCoefficients[t][0] = 1;

    for (m = 1; m <= n; m++)
        for (t = 1; t <= 100; t++)
            binCoefficients[m][t] = binCoefficients[m-1][t-1] + binCoefficients[m-1][t];


}

/**
 * @brief Finds the pairs (u,m) such that
 * the tables Xl and Xr have a size whose difference
 * is below threshold
 * 
 * @param n The vectors length to be considered
 * @param t The maximum number of errors
 * @param threshold The maximum difference between the sizes (ratio big_table over small_table)
 * @param E Pointer to the list the found entries
 */
void findEqualSize_u_m(int n, int t, float threshold, PairSet* E){

    int u,m;
    int size_l, size_r;

    for (m=n-1; m>0;--m){
        /*
            (n-m)!/t!(n-m-t)!= (n-m+1)!/t!(n-m-t+1)! * (n-m-t+1)/(n-m+1)
        */

        for (u=0; u<=m && u<=t && t-u<=n-m;++u){

            size_l=binCoefficients[m][u];
            size_r=binCoefficients[n-m][t-u];

            float ratio = size_r >= size_l ? size_r *1.0/size_l : size_l*1.0/size_r;
            if ( ratio < threshold ){
                //printf("Found pair: u=%d, m=%d\nSize(Xl)=%d\nSize(Xl)=%d\n",u,m,size_l,size_r);
                PairSet_addHead(u,m,E);
            }
        }
    }

}

bool updateIndexes(int* indexes, int* moduli, int u){

    int i;
    int carry=1;
    bool clean=false;

    /*puts("Before");
    for(i=0;i<u;++i)
        printf("%d",indexes[i]);
    printf("\n");*/

    for (i=u-1; i >=0; --i){

        indexes[i] = (indexes[i]+carry) % moduli[i];
        carry = indexes[i] == 0 ? 1:0;
    }

    int count;
    for (i=0; i<u; ++i){

        if (indexes[i]!=0 && indexes[i+1]==0 && !clean && i!=u-1  ){
            clean = true;
            count=indexes[i];
            continue;
        }
        
        if (clean){
            indexes[i]=++count;
        }

    }

    /*for(i=0;i<u;++i)
        printf("%d",indexes[i]);
    printf("\n");*/

    // Returns a bool flag saying if you should stop
    return !(indexes[0]==0 && indexes[u-1]==0);
}

/**
 * @brief Enumerates all the m-vectors with a number of 1s equals to u
 * 
 * @param m 
 * @param u 
 */
void iterateOverM_Vectors(int m, int u, BinMatrix H_l_r, BST* X){

    int* indexes = malloc(sizeof(int)*u);
    int* moduli = malloc(sizeof(int)*u);
    int* array = malloc(sizeof(int)*m);

    int i,j,count=0;

    // Case when u=0 or u=m, only one vector is considered
    if (u==0 || u==m){
        int elem_to_set = (u==0? 0:1);
        for(i=0;i<m;++i){
            array[i]=elem_to_set;
        }
        BinMatrix* e = buildMatrix(array,1,m);
        BinMatrix* s = product(*e,H_l_r);
        addNode((void*)s, (void*) e, X,BST_COMPARISON_BINMATRIX);
        free(indexes);
        free(moduli);
        free(array);
        return;
    }

    for (i=0;i<u;++i){
        indexes[i]=i;
        moduli[i]=m-(u-i)+1;
    }

    do
    {
        for (i=0,j=0;i<m;++i){
            
            if (j < u && i==indexes[j]){
                array[i]=1;
                j++;
            }
            else
                array[i]=0;
                
        }

        BinMatrix* e = buildMatrix(array,1,m);
        BinMatrix* s = product(*e,H_l_r);
        addNode((void*)s, (void*) e, X,BST_COMPARISON_BINMATRIX);
        count++;
    }while (updateIndexes(indexes,moduli,u));
    

    //printf("Tested %d vectors\n",count);
    free(indexes);
    free(moduli);
    free(array);

}

/**
 * @brief Converts an integer to binary vector
 * 
 * @param x The integer to convert
 * @param v The vector to populate
 * @param len The length of v
 * @returns The number of 1s in the vector
 */
int intToBinVector(int x, int* v, int len){
    int count =0;
    int res=0;

    for (int i=0; i<len; ++i){
        res = x % 2;
        v[len-i-1]= res;
        count += res;
        x = x >> 1;
    }

    return count;
}

/**
 * @brief Obtains the left syndrome from the right syndrome and the full syndrome
 * 
 * @param s 
 * @param sr 
 * @return BinMatrix* 
 */
BinMatrix* getLeftSyndrome(BinMatrix s, BinMatrix sr){

    return vectorSum(s,sr);
}

/**
 * @brief Implements step 4 of the split-syndrome algorithm
 * 
 * @param Xr The right table
 * @param Xl The left table
 * @param s The full syndrome
 * @param el Contains the left error at the end of the algortithm
 * @param er Contains the right error at the end of the algortithm
 * @return true 
 * @return false 
 */
bool inspectTables(BST Xr, BST Xl, BinMatrix s, VectorList* el, VectorList* er){

    BSTNode* node;
    bool found=false;

    if (Xr->l != NULL){
        found = inspectTables(Xr->l,Xl,s, el,er);

        if (found)
            return found;
    }

    BinMatrix* left_syndrome = getLeftSyndrome(s, *(BinMatrix*) Xr->key);
    node = searchNode((void*)left_syndrome,Xl,BST_COMPARISON_BINMATRIX);

    if (node != NULL){

        VectorList temp = (VectorList)Xr->data;
        while(temp!=NULL)
           VectorList_addHead(vectorList_pop(&temp)->v, er);
        
        temp = (VectorList)node->data;
        while(temp!=NULL)
           VectorList_addHead(vectorList_pop(&temp)->v, el);

        /*
        *er = (VectorList) Xr->data;
        *el = (VectorList) node->data;
        //printf("FOUND!\n");
        */
        return true;
    }
    else{
        destroyMatrix(left_syndrome);
    }

    if (Xr->r != NULL){
        found = inspectTables(Xr->r,Xl,s,el,er);

        return found;
    }
    else
        return false;
}

/**
 * @brief Builds the left table Xl
 * 
 * @param m The length of the left side
 * @param u The number of errors considered
 * @param H The parity check matrix
 * @param Xl The left table to populate
 */
void buildLeftTable(int m,int u, BinMatrix H, BST* Xl){

    int i;
    int* indexes= malloc(sizeof(int)*m);
    int* error_as_vector= malloc(sizeof(int)*m);

    for (i=0; i<m;++i){
        indexes[i]=i;
        error_as_vector[i]=0;
    }

    BinMatrix* samples = sampleFromMatrix(indexes,m,H,MATRIX_SAMPLE_COLUMNS);
    BinMatrix* Hl = transpose(*samples);

    /*for (i=0; i<pow(2,m);++i){

        if ( intToBinVector(i,error_as_vector,m) != u )
            continue;
        BinMatrix* el = buildMatrix(error_as_vector,1,m);
        BinMatrix* sl = product(*el,*Hl);
        addNode((void*)sl, (void*) el, Xl,BST_COMPARISON_BINMATRIX);
    }*/

    iterateOverM_Vectors(m,u,*Hl,Xl);

    destroyMatrix(samples);
    destroyMatrix(Hl);
    free(indexes);
    free(error_as_vector);
}

/**
 * @brief Builds the right table Xr
 * 
 * @param m The length of the left part
 * @param t_minus_u The number of errors in the right part, i.e. t-u
 * @param len The total syndrome length
 * @param H Parity check matrix
 * @param Xr The table to populate
 */
void buildRightTable(int m,int t_minus_u, int len, BinMatrix H, BST* Xr){

    int i;
    int* indexes= malloc(sizeof(int)*(len-m));
    int* error_as_vector= malloc(sizeof(int)*(len-m));

    for (i=0; i<len-m;++i){
        indexes[i]=i+m;
        error_as_vector[i]=0;
    }

    BinMatrix* samples = sampleFromMatrix(indexes,len-m,H,MATRIX_SAMPLE_COLUMNS);
    BinMatrix* Hr = transpose(*samples);

    /*for (i=0; i<pow(2,len-m);++i){

        if ( intToBinVector(i,error_as_vector,len-m) != t_minus_u)
            continue;
        
        BinMatrix* er = buildMatrix(error_as_vector,1,len-m);
        BinMatrix* sr = product(*er,*Hr);
        addNode((void*)sr, (void*) er, Xr,BST_COMPARISON_BINMATRIX);
    }*/

    iterateOverM_Vectors(len-m,t_minus_u,*Hr,Xr);

    destroyMatrix(samples);
    destroyMatrix(Hr);
    free(indexes);
    free(error_as_vector);
}

/**
 * @brief Implements the split-syndrome algorithm
 * 
 * @param H
 * @param s 
 * @param d 
 */
void SplitSyndrome(BinMatrix H, BinMatrix s, int d, VectorList* left,VectorList* right){

    int t,u,m;
    BST Xl=NULL;
    BST Xr=NULL;
    /*VectorList el;
    VectorList er;*/
    PairSet* tables=malloc(sizeof(PairSet)*(d+1));

    // First do the precomputation step
    //printf("Starting the precomputation stage\n");

    // Compute some bin coefficients (already computed in main)
    // precomputeBinCoefficients(100,100);

    // Build the tables E(1),...,E(d)
    for (t=1; t<=d;++t){
        tables[t]=NULL;
        findEqualSize_u_m(2*s.cols,t,3,&tables[t]);
        //printf("Found E(%d)\n",t);
    }

    //printf("Precomputation stage complete!\n");

    
    for (t=1; t<=d;++t){

        PairSet E = tables[t];
        //printf("Trying t=%d\n",t);

        while (E != NULL)
        {   
            // Get a (u,m) pair from E(t)
            Pair* pair=PairSet_pop(&E);
            u=pair->x;
            m=pair->y;

            PairSet_destroy(pair);

            //printf("Trying the values u=%d,m=%d,t=%d\n",u,m,t);
            // Build the two tables
            //printf("Building the table X_left...");
            buildLeftTable(m,u,H,&Xl);
            //printf("DONE!\nBuilding the table X_right...");
            buildRightTable(m,t-u,2*s.cols,H,&Xr);
            
            //printf("DONE!\n");

            if ( inspectTables(Xr,Xl,s,left,right) ){
                
                /* do nothing */;

                // If you found the two errors, build the full error and return it
                //printf("Found code!\n");
                //VectorList_print(el);
                //VectorList_print(er);
                
                /*
                *left=el;
                *right=er;
                return;
                */
            }

            destroyTree(&Xl,BST_COMPARISON_BINMATRIX);
            destroyTree(&Xr,BST_COMPARISON_BINMATRIX);
            Xl=NULL;
            Xr=NULL;
        }
        
    }

    free(tables);
}

#endif