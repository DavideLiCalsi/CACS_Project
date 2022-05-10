#include "Matrix.h"

/**
 * @brief Computes the syndrome of vector y
 * 
 * @param parityCheck The parity check matrix
 * @param y The row vector whose syndome you should compute
 * @return Matrix* A column  vector representing the syndome, NULL for error
 */
Matrix* computeSyndrome(Matrix parityCheck, Matrix y){

    if (y.cols != parityCheck.cols || y.rows != 1){
        printf("Error: input vector of size (%d,%d) is not a columns vector\n", y.rows,y.cols);
        return NULL;
    }

    Matrix* res = product(parityCheck, *transpose(y));

    return res;

}

/**
 * @brief Computes the parity check matrix from the
 * generator matrix. The generator matrix should be in standard form
 * G = [I_k | A]
 * 
 * @param A
 * @param identity_size The size of the identity matrix in the generator's
 * standard form
 * @return Matrix* Pointer to the parity check matrix, null for error
 */
Matrix* computeParityCheck(Matrix A, int identity_size){

    if (identity_size != A.rows){
        printf("Invalid inputs. Matrix A and I_k have %d and %d rows\n",A.rows,identity_size);
        return NULL;
    }

    Matrix* A_t = transpose(A);
    Matrix* I_k = identityMatrix(identity_size);
    return(concat(*A_t,*I_k,0));
}

