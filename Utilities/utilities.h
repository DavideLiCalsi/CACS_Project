#ifndef UTILITIES_HEADER
#define UTILITIES_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../Matrix/BinaryMatrix.h"
#include "../SplitSyndrome/SupercodeSplitSyndrome.h"
#include "dataStructures.h"
#include "randomSelector.h"


/**
 * @brief compute the Gilbert-Vashamov distance (d_0)
 * 
 * @param n code length
 * @param k dimension subspace of vector space F^n_q
 * @param q elements of finite field F
 * @return int 
 */
int gilbertVashamovDistance(int n, int k, int q){

    double y = pow(q, n-k);
    int distance = 1;

    double partial_sum = 0;
    while(partial_sum <= y){
        distance++;
        partial_sum += binomialCoeff(n,distance-1)*pow(q-1,distance-1);
    }
    return distance;
}


/**
 * @brief compute the L_n(k, e) parameter
 * 
 * @param n 
 * @param k 
 * @param e 
 * @return unsigned long int 
 */
unsigned long int computeLn(int n, int k, int e){

    //int GV_distance = 12;
    int GV_distance = gilbertVashamovDistance(n, k, 2);
    
    double result_d = (n*log(n)/log(2)) * binomialCoeff(n,GV_distance) /  (binomialCoeff(k,e) * binomialCoeff(n-k,GV_distance-e) );
    unsigned long int result = round(result_d);
    return round(result);

}



/**
 * @brief transform generator matrix to standard form
 * 
 * @param generator generator matrix
 * @return BinMatrix* generator matrix in standard (systematic) form
 */
BinMatrix* standardizeGeneratorMatrix(BinMatrix *generator){

    // generator is a (k x n) matrix
    int k = generator->rows;
    

    // goal: make the first kxk matrix an Identity (generator = g)
    for (int j = 0; j < k; j++){

        /** 
        * Step 1) If g_jj != 0, go to Step 2. If g_jj==0, and if for some i>j, g_ij != 0, then interchange row_j and row_i.
        *         If g_jj==0 and g_ij==0 for all i>j, then choose h such that g_jh != 0 and interchange col_j and col_h.
        *         (last case is impossibile because the k columns should be linear indipendent (belonging to the Information Set))
        */
        if (getElement(*generator,j,j) == 0){
            bool success = false;
            for (int i=j+1; i<k && !success; i++){
                if (getElement(*generator,i,j) != 0){
                    swapRows(generator,j,i);
                    success = true;
                }
            }
        }
        /**
         * Step 2) We now have g_jj != 0. Multiply row_j by g_jj(^-1) (useless working with binary code)
         */

        /**
         * Step 3) We now have g_jj==1. For each of i>=0 and i<k (i!=j), replace row_i by row_i - g_ij*row_j
         */
        for (int i=0; i<k; i++)
            if (getElement(*generator,i,j) != 0 && i!=j)
                addRows(generator, i, j);
    }
    return generator;
}



/**
 * @brief transform parity check matrix to standard form
 * 
 * @param parityMatrix parity check matrix
 * @return BinMatrix* parity check matrix in standard (systematic) form
 */
BinMatrix* standardizeParityMatrix(BinMatrix *parityMatrix){

    // parityMatrix is a (n-k) x n matrix
    int rows = parityMatrix->rows;
    int cols = parityMatrix->cols;
    

    // goal: make the first (n-k)x(n-k) matrix an Identity (parityMatrix = g)
    for (int z=0, j = rows; z < rows; z++, j++){

        /** 
        * Step 1) If g_zj != 0, go to Step 2. If g_zj==0, and if for some i>z, g_ij != 0, then interchange row_z and row_i.
        *         If g_zj==0 and g_ij==0 for all i>z, then choose h such that g_zh != 0 and interchange col_j and col_h.
        *         (last case is impossibile because the n-k columns should be linear indipendent (belonging to the Information Set))
        */
        if (getElement(*parityMatrix,z,j) == 0){
            bool success = false;
            for (int i=z+1; i<rows && !success; i++){
                if (getElement(*parityMatrix,i,j) != 0){
                    swapRows(parityMatrix,z,i);
                    success = true;
                }
            }
        }
        /**
         * Step 2) We now have g_zj != 0. Multiply row_z by g_zj(^-1) (useless working with binary code)
         */

        /**
         * Step 3) We now have g_zj==1. For each of i>=0 and i<k (i!=z), replace row_i by row_i - g_ij*row_z
         */
        for (int i=0; i<rows; i++)
            if (getElement(*parityMatrix,i,j) != 0 && i!=z)
                addRows(parityMatrix, i, z);
    }
    return parityMatrix;
}


/**
 * @brief generate an admissable codeword using generator matrix G
 * 
 * @param G generator matrix
 * @param seed random seed
 * @return BinMatrix* generated codeword
 */
BinMatrix* generateCodeword(BinMatrix *G, int seed){

    srand(seed);
    int k = G->rows;
    int n = G->cols;

    // take randomly some of the generator matrix's rows and combine them to obtain a codeword 
    int *rows_array = (int *)malloc(sizeof(int)*k);
    for (int i=0; i<k; i++)
        rows_array[i] = i;
    Set* G_rows = buildSet(rows_array, k);
    free(rows_array);
    int num_G_rows = 1 + rand() % (k-1); // number of the generator's rows to combine
    Set* G_rows_index = getDistinctRandomNumbers(G_rows, num_G_rows, seed); // indeces of the G's rows to select
    destroySet(G_rows);
    
    // generate the codeword
    BinMatrix *codeword = getRow(*G, G_rows_index->data[0]);
    
    BinMatrix *tmp, *tmp2;
    for (int i=1; i<num_G_rows; i++){
        tmp=codeword;
        tmp2 = getRow(*G, G_rows_index->data[i]);
        codeword = vectorSum(*tmp, *tmp2);
        destroyMatrix(tmp); destroyMatrix(tmp2);
    }
    destroySet(G_rows_index);

    return codeword;
}


/**
 * @brief generate a random error due to a noisy channel
 * 
 * @param n codeword's length
 * @param w minimum Hamming weight
 * @param seed random seed
 * @return BinMatrix* generated random error
 */
BinMatrix *generateError(int n, int w, int seed){

    srand(seed);

    // select number of 1s of the error
    int ones = 1 + rand() % (w-1);

    // select the position to set to 1
    int *array = (int *)malloc(sizeof(int)*n);
    for (int i=0; i<n; i++)
        array[i] = i;
    Set* indeces = buildSet(array, n);
    free(array);
    Set *selected_indeces = getDistinctRandomNumbers(indeces, ones, seed);
    destroySet(indeces);
    
    // generate the error
    int *error_array = (int *)malloc(sizeof(int)*n);
    for (int i=0; i<n; i++)
        error_array[i] = 0;
    for (int i=0; i<selected_indeces->length; i++)
        error_array[selected_indeces->data[i]] = 1;

    destroySet(selected_indeces);

    BinMatrix *error = buildMatrix(error_array, 1, n);
    free(error_array);
    return error;
}

typedef struct _allBinaryString{
    int **binaryStrings;
    int index;
}AllBinaryStrings;


AllBinaryStrings *buildAllBinaryStrings(int lenStrings){
    AllBinaryStrings *res = (AllBinaryStrings *)malloc(sizeof(AllBinaryStrings));

    res->binaryStrings = (int **)malloc(sizeof(int *)*pow(2,lenStrings));
    for(int i=0; i<pow(2,lenStrings); i++)
        res->binaryStrings[i] = (int *)malloc(sizeof(int)*lenStrings);
    res->index = 0;

}


void generateAllBinaryStrings(int n, int arr[], int i, AllBinaryStrings *res){
    if (i == n) {
        int *copy = (int *)malloc(sizeof(int)*n);
        for (int j=0; j<n; j++)
            copy[j] = arr[j];
        res->binaryStrings[res->index] = copy;
        res->index += 1;
        return;
    }
 
    // First assign "0" at ith position
    // and try for all other permutations
    // for remaining positions
    arr[i] = 0;
    generateAllBinaryStrings(n, arr, i + 1, res);
 
    // And then assign "1" at ith position
    // and try for all other permutations
    // for remaining positions
    arr[i] = 1;
    generateAllBinaryStrings(n, arr, i + 1, res);
}


/**
 * @brief generate all the admissible codeword
 * 
 * @param G Generator matrix
 * @return BinMatrix** Array of all the admissible codewords
 */
BinMatrix **generateAllCodeword(BinMatrix *G){

    int k = G->rows;
    BinMatrix **allCodewords = (BinMatrix **)malloc(sizeof(BinMatrix *));

    int *array = (int *)malloc(sizeof(int)*k);
    for(int i=0; i<k; i++)
        array[i] = i;
    Set *set = buildSet(array, k);
    free(array);

    AllBinaryStrings *allBinaryStrings = buildAllBinaryStrings(k);
    
    generateAllBinaryStrings(k, set->data, 0, allBinaryStrings);

    for (int i=0; i<pow(2,k); i++){
        for (int j=0; j<k; j++)
            printf("%d ", allBinaryStrings->binaryStrings[i][j]);
        printf("\n");
    }

    return NULL;
}





#endif