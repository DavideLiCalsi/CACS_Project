#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

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