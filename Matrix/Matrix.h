#include <stdio.h>
#include <stdlib.h>

#define MATRIX_SUCCESS 0
#define MATRIX_FAILURE 1

struct Matrix{
    int rows;
    int cols;
    int* data;
    __uint128_t* data;
};

typedef struct Matrix Matrix;

/**
 * @brief Get the Element of indexes (i,j)
 * 
 * @param m The matrix
 * @param i Row index
 * @param j Column index
 * @return int The element m[i][j]
 */
int getElement(Matrix m, int i, int j){

    if ( 0<=i && i<= m.rows-1 && 0<=j && j<=m.cols -1 )
        return m.data[i*m.cols + j];
    else{
        printf("Invalid indexes (%d,%d) for a matrix of size (%d,%d)!\n", i,j,m.rows,m.cols);
        return __INT_MAX__;
    }
}

/**
 * @brief Change m[i][j] to val
 * 
 * @param m The matrix
 * @param i The row index
 * @param j The column index
 * @param val The new values
 * @return int 0 on success, 1 otherwise
 */
int putElement(Matrix* m, int i, int j, int val){

    if ( 0<=i && i<= m->rows-1 && 0<=j && j<=m->cols -1 ){
        m->data[i*m->cols + j]=val;
        return MATRIX_SUCCESS;
    }
    else{
        printf("Invalid indexes (%d,%d) for a matrix of size (%d,%d)!\n", i,j,m->rows,m->cols);
        return MATRIX_FAILURE;
    }
}

/**
 * @brief Get the i-th row of the matrix
 * 
 * @param m The matrix
 * @param i The row index
 * @return int* An array representing the matrix row
 */
int* getRow(Matrix m, int i){

    if ( 0>i || i>=m.rows){
        printf("Invalid row index %d for a matrix having %d rows!\n", i, m.rows);
    }

    int* row = (int*)( malloc(sizeof(int) * m.cols) );

    for(int j=0; j<m.cols;++j){
        row[j]=getElement(m,i,j);
    }

    return row;
}

/**
 * @brief Get the j-th column of the matrix
 * 
 * @param m The matrix
 * @param j The column index
 * @return int* An array representing the matrix columns
 */
int* getColumn(Matrix m, int j){

    if ( 0>j || j>=m.cols){
        printf("Invalid column index %d for a matrix having %d colums!\n", j, m.cols);
    }

    int* column = (int*)( malloc(sizeof(int) * m.rows) );

    for(int i=0; i<m.rows;++i){
        column[i]=getElement(m,i,j);
    }

    return column;
}

/**
 * @brief Computes the transpose matrix of m
 * 
 * @param m Matrix that you wish to transpose
 * @return Matrix* pointer to the transpose matrix, NULL on failure
 */
Matrix* transpose(Matrix m){

    Matrix* m_t = (Matrix*) malloc(sizeof(Matrix));
    m_t->cols=m.rows;
    m_t->rows=m.cols;
    m_t->data=(int*) malloc(sizeof(int)*m.rows*m.cols);

    for(int i=0; i<m.rows;++i){
        for(int j=0; j<m.cols;++j){

            if ( putElement(m_t,j,i, getElement(m,i,j)) != MATRIX_SUCCESS)
                return NULL;
        }
    }

    return m_t;
}

/**
 * @brief Computes the product of the row vector v1 and the
 * column vector v2
 * 
 * @param v1 The row vector
 * @param v2 The column vector
 * @param size The size of the arrays
 * @return int The product
 */
int elementwiseProduct(int* v1, int* v2, int size){

    int res=0;

    for(int i=0; i<size;++i)
        res+=v1[i]*v2[i];
    
    return res;
}

/**
 * @brief Compute the matrix product between m1 and m2
 * 
 * @param m1 
 * @param m2 
 * @return Matrix* The product between m1 and m2, NULL if fails
 */
Matrix* product(Matrix m1, Matrix m2){

    if (m1.cols != m2.rows){
        printf("Cannot multiply matrices with %d columns and %d rows\n",m1.cols,m2.rows);
        return NULL;
    }
        
    Matrix* res = (Matrix*) (malloc(sizeof(Matrix)));
    res->rows=m1.rows;
    res->cols=m2.cols;
    res->data = (int*)(malloc(sizeof(int) * res->rows * res->cols));

    for(int i=0; i<res->rows; ++i){
        for(int j=0;j<res->cols; ++j){
            
            int val = elementwiseProduct( getRow(m1,i),getColumn(m2,j), m1.cols );
            if( putElement(res,i,j,val) != MATRIX_SUCCESS)
                return NULL;
        }
    }

    return res;
}

/**
 * @brief Pretty prints the input matrix
 * 
 * @param m The input matrix
 */
void printMatrix(Matrix m){

    printf("Matrix(%d rows, %d cols)\n",m.rows,m.cols);

    for(int i=0;i<m.rows;++i){
        for(int j=0;j<m.cols;++j){
            printf("%d ",getElement(m,i,j));
        }
        printf("\n");
    }

}

/**
 * @brief Build a matrix from an array
 * @param array The array to convert
 * @param rows The rows of the matrix
 * @param cols The columns of the matrix
 * @return Matrix* The ptr to the new matrix object, NULL for failure
 */
Matrix* buildMatrix(int* array, int rows, int cols){

    Matrix* m = (Matrix*)(malloc(sizeof(Matrix)));

    m->rows=rows;
    m->cols=cols;
    m->data=(int*) malloc(sizeof(int)*rows*cols);

    for(int i=0; i<rows*cols; ++i){
        
        int row_index, col_index;
        row_index=i/cols;
        col_index= (i%cols);
        if ( putElement(m, row_index, col_index,array[i]) != MATRIX_SUCCESS )
            return NULL;
    }

    return m;
}

/**
 * @brief Destroys the input matrix, i.e. it frees all the
 * allocated memory
 * 
 * @param m Pointer to the matrix to free
 */
void destroyMatrix(Matrix* m){

    if (m == NULL || m->data == NULL){
        printf("Error: cannot free null pointer\n");
        return;
    }

    free(m->data);
    free(m);

}

/**
 * @brief Returns a k x k identity matrix
 * 
 * @param k The size of the identity matriz
 * @return Matrix* The identity matrix, null for error
 */
Matrix* identityMatrix(int k){

    if (k<0){
        printf("Please specify a valid size, %d is invalid\n", k);
        return NULL;
    }

    Matrix* res = (Matrix*) malloc(sizeof(Matrix));
    res->rows=k;
    res->cols=k;
    res->data=(int*) (malloc(sizeof(int)*k*k));

    for (int i=0; i<k; ++i){
        for (int j=0; j<k; ++j){

            if (i==j)
                putElement(res,i,j,1);
            else
                putElement(res,i,j,0);
        }
    }

    return res;

}

/**
 * @brief Concatenates two matrices along the specified axis
 * 
 * @param m1 First matrix
 * @param m2 Second matrix
 * @param axis The axis along which you want to concat, either 0 or 1
 * @return Matrix* result, Null for errors
 */
Matrix* concat(Matrix m1, Matrix m2, int axis){

    Matrix* res = (Matrix*) malloc(sizeof(Matrix));

    switch (axis)
    {
    case 0:
        
        if (m1.rows != m2.rows){
            printf("Cannot concatenate two matrices having %d and %d rows\n", m1.rows, m2.rows);
            free(res);
            return NULL;
        }

        res->rows=m1.rows;
        res->cols=m1.cols+m2.cols;
        res->data=(int*) malloc(sizeof(int)*res->cols*res->rows);

        for (int i=0; i<res->rows; ++i){
            for (int j=0; j<res->cols; ++j){

                if (j<m1.cols)
                    putElement(res,i,j, getElement(m1,i,j));
                else
                    putElement(res,i,j,getElement(m2,i,j-m1.cols));
            }
        }

        return res;
        break;

    case 1:
        
        if (m1.cols != m2.cols){
            printf("Cannot concatenate two matrices having %d and %d colums\n", m1.cols, m2.cols);
            free(res);
            return NULL;
        }
        res->rows=m1.rows + m2.rows;
        res->cols=m1.cols;
        res->data=(int*) malloc(sizeof(int)*res->cols*res->rows);

        for (int i=0; i<res->rows; ++i){
            for (int j=0; j<res->cols; ++j){

                if (i<m1.rows)
                    putElement(res,i,j, getElement(m1,i,j));
                else
                    putElement(res,i,j,getElement(m2,i-m1.rows,j));
            }
        }

        return res;
        break;
    
    default:
        printf("Invalid axis %d\n",axis);
        return NULL;
        break;
    }
}
