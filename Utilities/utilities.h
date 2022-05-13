#ifndef UTILITIES_HEADER
#define UTILITIES_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../Matrix/BinaryMatrix.h"

/**
 * @brief compute the (n,k) binomial coefficient
 * 
 * @param n 
 * @param k 
 * @return int 
 */
int binomialCoefficients(int n, int k){
    if (k == 0 || k == n)
        return 1;
    return binomialCoefficients(n - 1, k - 1) + binomialCoefficients(n - 1, k);
}


/**
 * @brief compute the Gilbert-Vashamov distance (d_0)
 * 
 * @param n code length
 * @param k dimension subspace of vector space F^n_q
 * @param q elements of finite field F
 * @return int 
 */
int gilbertVashamovDistance(int n, int k, int q){

    int y = pow(q, n-k);
    int distance = 1;

    int partial_sum = 0;
    while(partial_sum < y){
        distance++;
        partial_sum += binomialCoefficients(n-1,distance-2)*pow(q-1,distance-2);        
    }
    return distance;
    
    /*
    implementation based on a different formula 
    (see "General Methods of Decoding of Linear Codes")
    distance = 1;
    partial_sum = 0;
    while(partial_sum < y){
        distance++;
        partial_sum += binomialCoefficients(n,distance-1)*pow(q-1,distance-1);
    }
    */
}


/**
 * @brief compute the L_n(k, e) parameter
 * 
 * @param n 
 * @param k 
 * @param e 
 * @return int 
 */
int computeLn(int n, int k, int e){

    int GV_distance = gilbertVashamovDistance(n, k, 2);
    return (n*log(n)/log(2)) * (n,GV_distance) / binomialCoefficients(k,e) * binomialCoefficients(n-k,GV_distance-e);

}



/**
 * @brief transform generator matrix to standard form
 * 
 * @param generator generator matrix
 * @return BinMatrix* generator matrix in standard (systematic) form
 */
BinMatrix* standardizeGeneratorMatrix(BinMatrix *generator){

    int k = generator->rows;    // generator is a (k x n) matrix
    

    // goal: make the first kxk matrix an Identity (generator = g)
    for (int j = 0; j < k; j++){

        /** 
        * Step 1) If g_jj != 0, go to Step 2. If g_jj==0, and if for some i>j, g_ij != 0, then interchange row_j and row_i.
        *         If g_jj==0 and g_ij==0 for all i>j, then choose h such that g_jh != 0 and interchange col_j and col_h.
        *         (last case is impossibile because the k columns should be linear indipendet (belonging to the Information Set))
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

    // parityMatrix is a (k x n) matrix
    int k = parityMatrix->rows;
    int n = parityMatrix->cols;
    

    // goal: make the first kxk matrix an Identity (parityMatrix = g)
    for (int z=0, j = n-k; z < k; z++, j++){

        /** 
        * Step 1) If g_zj != 0, go to Step 2. If g_zj==0, and if for some i>z, g_ij != 0, then interchange row_z and row_i.
        *         If g_zj==0 and g_ij==0 for all i>z, then choose h such that g_zh != 0 and interchange col_j and col_h.
        *         (last case is impossibile because the n-k columns should be linear indipendet (belonging to the Information Set))
        */
        if (getElement(*parityMatrix,z,j) == 0){
            bool success = false;
            for (int i=z+1; i<k && !success; i++){
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
        for (int i=0; i<k; i++)
            if (getElement(*parityMatrix,i,j) != 0 && i!=z)
                addRows(parityMatrix, i, z);
    }
    return parityMatrix;
}



#endif