#ifndef SPLIT_SYNDROME_H
#define SPLIT_SYNDROME_H

#include "../Tree/BST.h"
#include "PairSet.h"
#include <math.h>
#include <time.h>
#include "../Utilities/debug.h"

double binCoefficients[200][200];
int valid[200][200];

/**
 * @brief compute the (n,k) binomial coefficient
 * 
 * @param n 
 * @param k 
 * @return double 
 */
double binomialCoeff(unsigned int n, unsigned int k){
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
    double t1,t2;

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
 * is below threshold.
 * 
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

/**
 * @brief This function is called by iterateOverM_vectors and is 
 * used to algoritmically update an array of indexes. Each position
 * in the array represents a position in the error vector that contains
 * a 1. The algorithm updates this array in a way that allows us to generate
 * all the error patterns with weight=u.
 * 
 * @param indexes The array of indexes
 * @param moduli An array of moduli, used to update the indexes
 * @param u The desired error's weight
 * @return true If you have tried all the possible error patterns
 * @return false Otherwise
 */
bool updateIndexes(int* indexes, int* moduli, int u){

    int i;
    int carry=1;
    bool clean=false;

    for (i=u-1; i >=0; --i){

        indexes[i] = (indexes[i]+carry) % moduli[i];
        carry = indexes[i] == 0 ? 1:0;
    }

    int count;
    for (i=0; i<u; ++i){

        if (indexes[i]!=0 && (i==(u-1) || indexes[i+1]==0) && !clean && i!=u-1  ){
            clean = true;
            count=indexes[i];
            continue;
        }
        
        if (clean){
            indexes[i]=++count;
        }

    }

    // Returns a bool flag saying if you should stop
    return !(indexes[0]==0 && indexes[u-1]==0);
}

/**
 * @brief Enumerates all the m-vectors with a number of 1s equals to u,
 * computes their syndrome and stores both in a table
 * 
 * @param m The vector's length
 * @param u The desired error weight
 * @param H_l_r The parity-check matrix to use to compute the syndrome
 * @param X The table
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
        addNode((void*)s, (void*) e, X, BST_COMPARISON_BINMATRIX);
        destroyMatrix(e);
        free(indexes);
        free(moduli);
        free(array);
        return;
    }

    /*
    Initialize the arrays. 
    Indexes <- {0,1,...,u-1}
    Moduli <- {m-u+1,m-u+2,...,m}
    */
    for (i=0;i<u;++i){
        indexes[i]=i;
        moduli[i]=m-(u-i)+1;
    }

    /*
    Iteratively generate all the error patterns with u errors,
    compute the syndrome and store the pair in the target table
    */
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
        PRINTMATRIX_PTR(e);
        addNode((void*)s, (void*) e, X, BST_COMPARISON_BINMATRIX);
        destroyMatrix(e);
        count++;
    }while (updateIndexes(indexes,moduli,u));
    
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
 * @param s Full syndrome
 * @param sr Right syndrome
 * @return BinMatrix* Left syndrome
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
bool inspectTables(BST Xr, BST Xl, BinMatrix s, VectorList* el, VectorList* er, int e1, int e2){

    BSTNode* node;
    bool found=false;

    if (Xr->l != NULL){
        found = inspectTables(Xr->l,Xl,s, el,er,e1,e2);

        if (found)
            return found;
    }

    BinMatrix* left_syndrome = getLeftSyndrome(s, *(BinMatrix*) Xr->key);
    node = searchNode((void*)left_syndrome,Xl,BST_COMPARISON_BINMATRIX);

    if (node != NULL){
        
        while(node->data!=NULL){

            VectorList curr=vectorList_pop( (VectorList *)(&(node->data)) );
            if ( HammingWeight( *(curr->v) ) <= e1 )
                VectorList_addHead(curr->v, el);            
            VectorList_destroy(&curr);

        }

        while(Xr->data!=NULL){

            VectorList curr=vectorList_pop((VectorList *)(&(Xr->data)));

            if ( HammingWeight( *(curr->v) ) <= e2 )
                VectorList_addHead(curr->v, er);
            VectorList_destroy(&curr);
            
        }
        
        destroyMatrix(left_syndrome);
        return true;
    }
    else{
        destroyMatrix(left_syndrome);
    }

    if (Xr->r != NULL){
        found = inspectTables(Xr->r,Xl,s,el,er,e1,e2);

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
void SplitSyndrome(BinMatrix H, BinMatrix s, int d, VectorList* left,VectorList* right, int e1, int e2,int k, int y){

    int t,u,m;
    BST Xl=NULL;
    BST Xr=NULL;

    /* 
    We can skip the precomputation step since we already know
    the how to split the syndrome and the error.
    */

    for (t=0; t<=d;++t){

        int i=0;
        while (i<=k && i<=t && t-i<=y)
        {   
            /*
            We always split a syndrome of length k+y into
            s_l of length k and s_r of length y.
            */
            u=i;
            m=k;

            buildLeftTable(m,u,H,&Xl);
            buildRightTable(m,t-u,y+k,H,&Xr);

            if ( inspectTables(Xr,Xl,s,left,right,e1,e2) ){


                // If you found the two errors, build the full error and return it
                PRINTF("Found code!\n");
                PRINTVLIST(*left);
                destroyTree(&Xl,BST_COMPARISON_BINMATRIX);
                destroyTree(&Xr,BST_COMPARISON_BINMATRIX);
                return;

            }

            destroyTree(&Xl,BST_COMPARISON_BINMATRIX);
            destroyTree(&Xr,BST_COMPARISON_BINMATRIX);
            Xl=NULL;
            Xr=NULL;

            ++i;
        }
    }
}

#endif