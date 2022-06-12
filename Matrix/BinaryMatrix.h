#ifndef BINMATRIX_H
#define BINMATRIX_H
#define N 50

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#define MAX(a,b) (a>b?a:b)

#define MATRIX_SUCCESS 0
#define MATRIX_FAILURE 1
#define MATRIX_INVALID_ELEMENT 2

#define MATRIX_SAMPLE_ROWS 0
#define MATRIX_SAMPLE_COLUMNS 1

#define MATRIX_INVALID_WEIGHT -1

/**
 * @brief Macros to manage Matrix checks
 *
 */
#define isRowVector(x) (x.rows == 1)
#define isColumnVector(x) (x.cols == 1)
#define indexWithinBounds(m,i,j) (i<m.rows && i>=0 && j<m.cols && j>=0)
#define indexWithinBoundsPtr(m,i,j) (i<m->rows && i>=0 && j<m->cols && j>=0)
#define indexOutOfBounds(m,i,j) (i>=m.rows || i<0 && j>=m.cols && j<0)
#define indexOutOfBoundsPtr(m,i,j) (i>=m->rows && i<0 && j>=m->cols && j<0)
#define isSquareMatrix(m) (m.cols==m.rows)
#define isSquareMatrixPtr(m) (m->cols==m->rows)

struct BinMatrix{
    int rows;
    int cols;
    unsigned long* data;
};

typedef struct BinMatrix BinMatrix;



/**
 * @brief Get the Element of indexes (i,j)
 *
 * @param m The matrix
 * @param i Row index
 * @param j Column index
 * @return int The element m[i][j]
 */
char getElement(BinMatrix m, int i, int j){

    if ( indexWithinBounds(m,i,j) ){

        int array_index= (m.cols * i + j) / (8*sizeof(unsigned long));
        unsigned long array_selector = (m.cols * i + j) % (8*sizeof(unsigned long));
        array_selector = 8*sizeof(unsigned long)- array_selector - 1;
        // printf("Array index %d\nArray selector %d\n",array_index,array_selector);
        unsigned long row = m.data[array_index];
        // printf("%ld\n",row);
        unsigned long selector = 1UL << array_selector;
        return ( row & selector) != 0UL ? 1 : 0;

    }
    else{
        printf("Invalid indexes (%d,%d) for a matrix of size (%d,%d)!\n", i,j,m.rows,m.cols);
        return 2;
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
int putElement(BinMatrix* m, int i, int j, int val){

    if ( indexWithinBoundsPtr(m,i,j) ){

        int array_index = (m->cols * i + j) / (8*sizeof(unsigned long));
        unsigned long array_selector = (m->cols * i + j) % (8*sizeof(unsigned long));
        array_selector = 8*sizeof(unsigned long) - array_selector - 1;

        switch (val)
        {
        case 1:
            m->data[array_index] |= (1UL << array_selector);
            return MATRIX_SUCCESS;
            break;

        case 0:
            m->data[array_index] &= ~(1UL << array_selector);
            return MATRIX_SUCCESS;
            break;

        default:
            printf("Bit must be 0 or 1. %d is not a valid bit\n", val);
            return MATRIX_FAILURE;
            break;
        }
    }
    else{
        printf("Invalid indexes (%d,%d) for a matrix of size (%d,%d)!\n", i,j,m->rows,m->cols);
        return MATRIX_FAILURE;
    }
}



/**
 * @brief Pretty prints the input matrix
 *
 * @param m The input matrix
 */
void printMatrix(BinMatrix m){

    printf("Matrix(%d rows, %d cols)\n",m.rows,m.cols);

    for(int i=0;i<m.rows;++i){
        for(int j=0;j<m.cols;++j){
            printf("%d ",getElement(m,i,j));
        }
        printf("\n");
    }

}



/**
 * @brief Compare two row vectors. Comparison is done
 * by interpreting the two vectors as binary integers
 *
 * @param v1
 * @param v2
 * @return int 0 if the content is equal, 1 if v1>v2, -1 if v1<v2, 2 for error
 */
int compareVectors(BinMatrix v1, BinMatrix v2){

    if (!isRowVector(v1) || !isRowVector(v2)){
        printf("I can only compare two row vectors\n");
        return MATRIX_INVALID_ELEMENT;
    }

    if (v1.cols != v2.cols){
        //printf("Cannot compare two vectors of size %d and %d\n", v1.cols,v2.cols);
        return MATRIX_INVALID_ELEMENT;
    }

    int ulong_needed = ceil( v1.cols*1.0 / (8*sizeof(unsigned long)) );
    int excess = v1.cols % (8*sizeof(unsigned long));
    unsigned long w1,w2;

    for (int i=0; i<ulong_needed;++i){

        w1 = v1.data[i];
        w2 = v2.data[i];

        if (i<ulong_needed-1){

            if (w1 > w2)
                return 1;

            if (w1<w2)
                return -1;
        }
        else{

            w1 = v1.data[i] >> (8*sizeof(unsigned long)-excess);
            w2 = v2.data[i] >> (8*sizeof(unsigned long)-excess);

            if (w1 > w2)
                return 1;

            if (w1<w2)
                return -1;

            if (w1==w2)
                return 0;
        }
    }
}

/**
 * @brief Compare the two matrices
 *
 * @param m1 First matrix
 * @param m2 Second matrix
 * @return true The matrices are equal
 * @return false They are different
 */
bool compareMatrices(BinMatrix m1, BinMatrix m2){

    if (m1.cols!=m2.cols || m1.rows != m2.rows)
        return false;

    int needed_ulongs = ceil( (m1.rows * m1.cols*1.0) / (8*sizeof(unsigned long)) );

    for (int i=0; i<needed_ulongs;++i){
        if (m1.data[i] != m2.data[i])
            return false;
    }

    return true;
}



/**
 * @brief Get the i-th row of the matrix
 *
 * @param m The matrix
 * @param i The row index
 * @return BinMatrix* An array representing the matrix row
 */
BinMatrix* getRow(BinMatrix m, int i){

    if ( 0>i || i>=m.rows){
        printf("Invalid row index %d for a matrix having %d rows!\n", i, m.rows);
    }

    BinMatrix* row = (BinMatrix*)( malloc(sizeof(BinMatrix)) );
    row->rows=1;
    row->cols=m.cols;

    int needed_u_long = ceil ( row->rows*row->cols*1.0 / (8*sizeof(unsigned long)) );
    row->data= (unsigned long*) malloc(sizeof(unsigned long) * needed_u_long);
    memset(row->data,0,sizeof(unsigned long) * needed_u_long);

    for(int j=0; j<m.cols;++j){
        unsigned long new = (unsigned long) getElement(m,i,j);
        int array_offset = j / (8*sizeof(unsigned long) );
        int bin_offset = j % (8*sizeof(unsigned long) );
        row->data[array_offset] |= new << (63-bin_offset);
    }

    return row;
}

/**
 * @brief Get the j-th column of the matrix
 *
 * @param m The matrix
 * @param j The column index
 * @return BinMatrix* An array representing the matrix column
 */
BinMatrix* getColumn(BinMatrix m, int j){

    if ( 0>j || j>=m.cols){
        printf("Invalid column index %d for a matrix having %d columns!\n", j, m.cols);
    }

    BinMatrix* col = (BinMatrix*)( malloc(sizeof(BinMatrix)) );
    col->rows=m.rows;
    col->cols=1;

    int needed_u_long = ceil( (col->rows*col->cols)*1.0 / (8*sizeof(unsigned long)) );
    col->data= (unsigned long*) malloc(sizeof(unsigned long) * needed_u_long);
    memset(col->data,0,sizeof(unsigned long) * needed_u_long);
    

    for(int i=0; i<m.rows;++i){
        unsigned long new = (unsigned long) getElement(m,i,j);
        int array_offset = i/ (8*sizeof(unsigned long) );
        int bin_offset = i % (8*sizeof(unsigned long) );
        col->data[array_offset] |= new << (63- bin_offset);
    }

    return col;
}

/**
 * @brief Computes the transpose matrix of m
 *
 * @param m BinMatrix that you wish to transpose
 * @return BinMatrix* pointer to the transpose matrix, NULL on failure
 */
BinMatrix* transpose(BinMatrix m){

    BinMatrix* m_t = (BinMatrix*) malloc(sizeof(BinMatrix));
    m_t->cols=m.rows;
    m_t->rows=m.cols;
    m_t->data=(unsigned long*) malloc(sizeof(unsigned long)*m.rows*m.cols);

    for(int i=0; i<m.rows;++i){
        for(int j=0; j<m.cols;++j){

            if ( putElement(m_t,j,i, getElement(m,i,j)) != MATRIX_SUCCESS)
                return NULL;
        }
    }

    return m_t;
}


/**
 * @brief Build a matrix from an array
 * @param array The array to convert
 * @param rows The rows of the matrix
 * @param cols The columns of the matrix
 * @return BinMatrix* The ptr to the new matrix object, NULL for failure
 */
BinMatrix* buildMatrix(int* array, int rows, int cols){

    BinMatrix* m = (BinMatrix*)(malloc(sizeof(BinMatrix)));
    unsigned long ulong_needed = ceil ( rows*cols*1.0/ (8*sizeof(unsigned long)) );
    m->rows=rows;
    m->cols=cols;
    m->data=(unsigned long*) malloc(sizeof(unsigned long)*ulong_needed);

    for(int i=0; i<rows*cols; ++i){

        int row_index, col_index;
        row_index=i/cols;
        col_index= (i%cols);
        int val = array[i];

        if ( val != 0 && val != 1){
            printf("Invalid entry %d in position %d for binary matrix\n",val,i);
            return NULL;
        }

        if ( putElement(m, row_index, col_index,val) != MATRIX_SUCCESS )
            return NULL;
    }

    return m;
}


/**
 * @brief return a copy of the input matrix
 * 
 * @param matrix input matrix
 * @return BinMatrix* copy of the input matrix
 */
BinMatrix* copyMatrix(BinMatrix *matrix){

    int rows = matrix->rows;
    int cols = matrix->cols;

    BinMatrix* m = (BinMatrix*)(malloc(sizeof(BinMatrix)));
    unsigned long ulong_needed = ceil( matrix->rows*cols*1.0/ (8*sizeof(unsigned long)) );
    m->rows=rows;
    m->cols=cols;
    m->data=(unsigned long*) malloc(sizeof(unsigned long)*ulong_needed);

    for(int i=0; i<rows; ++i)
        for (int j=0; j<cols; j++)
        if (putElement(m, i, j,getElement(*matrix,i,j)) != MATRIX_SUCCESS )
            return NULL;

    return m;

}

/**
 * @brief Destroys the input matrix, i.e. it frees all the
 * allocated memory
 *
 * @param m Pointer to the matrix to free
 */
void destroyMatrix(BinMatrix* m){

    if (m == NULL || m->data == NULL){
        printf("Error: cannot free null pointer\n");
        return;
    }

    free(m->data);
    free(m);

}

/**
 * @brief Returns a row vector of length k whose
 * entries are all 1
 *
 * @param k
 * @return BinMatrix*
 */
BinMatrix* oneVector(int k){

    BinMatrix* res = malloc(sizeof(BinMatrix));
    res->rows=1;
    res->cols=k;

    int ulong_needed = ceil( k*1.0 / (8*sizeof(unsigned long)) );
    res->data = malloc(sizeof(unsigned long)*ulong_needed);

    for (int i=0;i<ulong_needed;++i)
        res->data[i]=0xffffffffffffffffUL;

    return res;
}

/**
 * @brief Returns a row vector of length k whose
 * entries are all 0
 *
 * @param k
 * @return BinMatrix*
 */
BinMatrix* zeroVector(int k){

    BinMatrix* res = malloc(sizeof(BinMatrix));
    res->rows=1;
    res->cols=k;

    int ulong_needed = ceil( k*1.0 / (8*sizeof(unsigned long)) );
    res->data = malloc(sizeof(unsigned long)*ulong_needed);

    for (int i=0;i<ulong_needed;++i)
        res->data[i]=0;

    return res;
}

/**
 * @brief Returns a k x k identity matrix
 *
 * @param k The size of the identity matriz
 * @return Matrix* The identity matrix, null for error
 */
BinMatrix* identityMatrix(int k){

    if (k<0){
        printf("Please specify a valid size, %d is invalid\n", k);
        return NULL;
    }

    int ulong_needed = ceil ( (k*k) *1.0/ (8*sizeof(unsigned long)) );

    //printf("%d\n",ulong_needed);
    BinMatrix* res = (BinMatrix*) malloc(sizeof(BinMatrix));
    res->rows=k;
    res->cols=k;
    res->data=(unsigned long*) (malloc(sizeof(unsigned long)*ulong_needed));

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
BinMatrix* concat(BinMatrix m1, BinMatrix m2, int axis){

    BinMatrix* res = (BinMatrix*) malloc(sizeof(BinMatrix));
    int ulong_needed;

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
        ulong_needed = ceil( (res->cols*res->rows*1.0) /(8*sizeof(unsigned long)) );
        res->data=(unsigned long*) malloc(sizeof(unsigned long)*ulong_needed);

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
        ulong_needed = ceil( (res->cols*res->rows*1.0) /(8*sizeof(unsigned long)) );
        res->data=(unsigned long*) malloc(sizeof(unsigned long)*ulong_needed);

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

/**
 * @brief Sums two vectors
 *
 * @param v1 First operand
 * @param v2 Second operand
 * @return BinMatrix* Result, NULL for error
 */
BinMatrix* vectorSum(BinMatrix v1, BinMatrix v2){

    if (v1.rows != v2.rows || v1.cols != v2.cols){
        printf("Cannot sum two vectors of size (%d,%d) and (%d,%d)\n",v1.rows,v1.cols,v2.rows,v2.cols);
        return NULL;
    }

    int ulong_needed = ceil( (v1.cols*v1.rows*1.0) / (8*sizeof(unsigned long)) );

    BinMatrix* z = (BinMatrix*) malloc(sizeof(BinMatrix));
    z->rows=v1.rows;
    z->cols=v1.cols;
    z->data=(unsigned long*) malloc(sizeof(unsigned long)*ulong_needed);

    for (int i=0; i<ulong_needed;++i){

        z->data[i] = v1.data[i] ^ v2.data[i];
    }

    return z;
}

/**
 * @brief Computes the inner (scalar) product of the two vectors
 *
 */
char vectorProduct(BinMatrix v1, BinMatrix v2){

    if ( !( isRowVector(v1) && isColumnVector(v2) ) ){
        printf("Error. Scalar product can be computed on a row vector and a column vector only");
        return MATRIX_INVALID_ELEMENT;
    }

    if ( v1.cols != v2.rows ){
        printf("Error. Incompatible sizes (%d,%d) and (%d,%d)\n", v1.rows,v1.cols,v2.rows,v2.cols);
        return MATRIX_INVALID_ELEMENT;
    }

    int ulong_needed = ceil( (v1.cols*v1.rows*1.0) / ( 8*(sizeof(unsigned long)) ) );
    int j;
    char res = 0;

    for(int i=0; i<ulong_needed; ++i){

        unsigned long and = v1.data[i] & v2.data[i];

        bool condition = (i==(ulong_needed-1)) && (v1.cols%(8*sizeof(unsigned long))!= 0);
        int stop = condition ? v1.cols%(8*sizeof(unsigned long)) : 8*sizeof(unsigned long);
        for (j=0; j<stop; ++j){

            res ^= (and >> (8*sizeof(unsigned long)-j-1) ) & 1;
        }
    }
    
    
    return res;
}

/**
 * @brief Compute the matrix product between m1 and m2
 *
 * @param m1
 * @param m2
 * @return BinMatrix* The product between m1 and m2, NULL if fails
 */
BinMatrix* product(BinMatrix m1, BinMatrix m2){

    if (m1.cols != m2.rows){
        printf("Cannot multiply matrices with %d columns and %d rows\n",m1.cols,m2.rows);
        return NULL;
    }

    BinMatrix* res = (BinMatrix*) (malloc(sizeof(BinMatrix)));
    res->rows=m1.rows;
    res->cols=m2.cols;
    int ulong_needed = ceil( (res->rows * res->cols*1.0) / (8*sizeof(unsigned long)) );
    res->data = (unsigned long*)(malloc(sizeof(unsigned long) * ulong_needed));

    for(int i=0; i<res->rows; ++i){
        for(int j=0;j<res->cols; ++j){

            BinMatrix* row=getRow(m1,i);
            BinMatrix* column=getColumn(m2,j);
            int val = vectorProduct(*row,*column);
            destroyMatrix(row);
            destroyMatrix(column);
            if( putElement(res,i,j,val) != MATRIX_SUCCESS)
                return NULL;
        }
    }

    return res;
}

/**
 * @brief Swaps rows r1 and r2 in the input matrix
 *
 * @param m The input matrix
 * @param r1 The 1st row
 * @param r2 The 2nd row
 * @return int error code
 */
int swapRows(BinMatrix* m, int r1, int r2){

    if (r1 <0 || r1 > m->rows || r2<0 || r2>m->rows){
        printf("Invalid rows to swap %d and %d for matrix of size (%d,%d)\n", r1, r2, m->rows,m->cols);
        return MATRIX_INVALID_ELEMENT;
    }

    int row_len=m->cols;

    unsigned long bitmask=0, bitmask1, bitmask2;

    for (int i=0;i<row_len;++i){
        bitmask |= (1<<i);
    }

    int r1_array_index = row_len*(r1) / (8*sizeof(unsigned long));
    int r2_array_index = row_len*(r2) / (8*sizeof(unsigned long));
    int r1_bit_index = (row_len*r1) % (8*sizeof(unsigned long));
    int r2_bit_index = (row_len*r2) % (8*sizeof(unsigned long));

    bitmask1 = bitmask << ( (8*sizeof(unsigned long)) - row_len - r1_bit_index );
    bitmask2 = bitmask << ( (8*sizeof(unsigned long)) - row_len - r2_bit_index );

    // Extract the two rows
    unsigned long exctracted_row1 = m->data[r1_array_index] & bitmask1;
    unsigned long exctracted_row2 = m->data[r2_array_index] & bitmask2;

    if (r1 < r2 )
    {
        exctracted_row2 = exctracted_row2 << row_len*(r2-r1);
        exctracted_row1 = exctracted_row1 >> row_len*(r2-r1);
    }

    if (r2_bit_index < r1_bit_index)
    {
        exctracted_row2 = exctracted_row2 >> row_len*(r1-r2);
        exctracted_row1 = exctracted_row1 << row_len*(r1-r2);
    }

    m->data[r1_array_index] = ( m->data[r1_array_index] & (-1UL ^ bitmask1) );
    m->data[r1_array_index] |= exctracted_row2;
    m->data[r2_array_index] = ( m->data[r2_array_index] & (-1UL ^ bitmask2) ) | exctracted_row1;

    return 0;
}

/**
 * @brief Adds row r2 to row r1 in matrix m.
 * In this context, adding is binary (i.e. xor)
 *
 * @param m Pointer to the matrix to manipulate
 * @param r1 Row to increment
 * @param r2 Row to add
 * @return int
 */
int addRows(BinMatrix* m, int r1, int r2){

    // TODO add check of row indexes
    int row_len=m->cols;

    unsigned long bitmask=0, bitmask1, bitmask2;

    for (int i=0;i<row_len;++i){
        bitmask |= (1<<i);
    }

    int r1_array_index = row_len*(r1) / (8*sizeof(unsigned long));
    int r2_array_index = row_len*(r2) / (8*sizeof(unsigned long));
    int r1_bit_index = (row_len*r1) % (8*sizeof(unsigned long));
    int r2_bit_index = (row_len*r2) % (8*sizeof(unsigned long));

    bitmask1 = bitmask << ( (8*sizeof(unsigned long)) - row_len - r1_bit_index );
    bitmask2 = bitmask << ( (8*sizeof(unsigned long)) - row_len - r2_bit_index );

    // Extract the two rows
    unsigned long exctracted_row1 = m->data[r1_array_index] & bitmask1;
    unsigned long exctracted_row2 = m->data[r2_array_index] & bitmask2;

    if (r1 < r2 )
    {
        exctracted_row2 = exctracted_row2 << row_len*(r2-r1);
    }

    if (r2_bit_index < r1_bit_index)
    {
        exctracted_row2 = exctracted_row2 >> row_len*(r1-r2);
        exctracted_row1 = exctracted_row1 << row_len*(r1-r2);
    }

    m->data[r1_array_index] ^= exctracted_row2;

}

/**
 * @brief Computes the determinant but inefficiently
 *
 * @param m
 * @return char
 */
char addRowsSlow(BinMatrix* m,int i, int j){

    int col_index=0;
    int row_len=m->cols;

    while (col_index<row_len)
    {   
        char new = ( getElement(*m,i,col_index) ^ getElement(*m,j,col_index) );
        putElement(m,i,col_index,new);
        col_index++;
    }

}

void swapRowsSlow(BinMatrix* m, int i, int j){
    int col_index=0;
    int row_len=m->cols;

    while (col_index<row_len)
    {
        char i_val= getElement(*m,i,col_index);
        char j_val= getElement(*m,j,col_index);
        putElement(m,i,col_index,j_val);
        putElement(m,j,col_index,i_val);
        col_index++;
    }
}

BinMatrix* GaussElimination(BinMatrix m){
    int i,j,k;

   // puts("BEGIN DET");
    for (j=0; j<m.cols;++j){

        k=j;

        /*
        With this loop, you turn the j-th column in the form [1,1,...,1,0,0,...,0]
        */
        for (i=j; i<m.rows;++i){

            if ( getElement(m,i,j) == 1 ){
                swapRowsSlow(&m,k,i);
                k++;
            }
        }

        /*
        Sum rows to bring the matrix in inferior triangular form
        */
        for (i=j+1; i<m.rows;++i ){

            if ( getElement(m,i,j) == 1 ){

                addRowsSlow(&m,i,j);
            }
        }
    }

    return copyMatrix(&m);
}

/**
 * @brief Computes the determinant of matrix m
 *
 * @param m
 * @return char
 */
char determinant(BinMatrix m){

    int i,j,k;

    if (!isSquareMatrix(m)){
        printf("Matrix of size (%d,%d) is not a square matrix\n", m.rows,m.cols);
        return MATRIX_INVALID_ELEMENT;
    }
   // puts("BEGIN DET");
    for (j=0; j<m.cols;++j){

        k=j;

        /*
        With this loop, you turn the j-th column in the form [1,1,...,1,0,0,...,0]
        */
        for (i=j; i<m.rows;++i){

            if ( getElement(m,i,j) == 1 ){
                swapRowsSlow(&m,k,i);
                k++;
            }
        }

     //   printf("Fixed column %d\n",j);

        /*
        Sum rows to bring the matrix in inferior triangular form
        */
        for (i=j+1; i<m.rows;++i ){

            if ( getElement(m,i,j) == 1 ){
              //  printf("Adding %d to %d\n",j,i);
                addRowsSlow(&m,i,j);
            }
        }
    }

    //printf("DONE DET\n");
    //printMatrix(m);
    /*
    If there is a 0 on the main diagonal, the determinant is 0
    */
    for (i=0;i<m.rows;++i){
        if (getElement(m,i,i) == 0)
            return 0;
    }

    return 1;
}

/**
 * @brief Computes the inverse of matrix m
 * through the Gauss-Jordan elimination
 *
 * @param m The matrix to invert
 * @return BinMatrix* The inverse matrix
 */
BinMatrix* inverse(BinMatrix m){

    int i,j,k;

    if (!isSquareMatrix(m)){
        printf("Matrix of size (%d,%d) is not a square matrix\n", m.rows,m.cols);
        return NULL;
    }

    // Build the augmented matrix by concatenating the matrix to invert and the Identity
    BinMatrix *identity = identityMatrix(m.rows);
    BinMatrix *augmented = concat(m,*identity,0);
    destroyMatrix(identity);

    for (j=0; j<m.cols;++j){

        k=j;

        /*
        With this loop, you turn the j-th column in the form [1,1,...,1,0,0,...,0]
        */
        for (i=j; i<m.rows;++i){

            if ( getElement(*augmented,i,j) == 1 && k !=i){
                swapRowsSlow(augmented,k,i);
                k++;
            }
        }

        /*
        Sum rows to bring the matrix in inferior triangular form
        */
        for (i=j+1; i<m.rows;++i ){

            if ( i!=j && getElement(*augmented,i,j) == 1 )
                addRowsSlow(augmented,i,j);
        }
    }

    //Now obtain the identity

    for (i=1;i<m.rows;++i){

        for(j=0;j<i;++j){

            if (i!=j && getElement(*augmented,j,i)==1){
                addRowsSlow(augmented,j,i);
            }

        }
    }

    //printMatrix(*augmented);

    BinMatrix* inv = (BinMatrix*) ( malloc(sizeof(BinMatrix)) );
    inv->rows=m.rows;
    inv->cols=m.cols;
    inv->data = (unsigned long*) ( malloc( sizeof(unsigned long) * inv->rows * inv->cols ) );

    /*
    This loop copies bits one by one from the augmented matrix
    to the matrix to return.
    NOTE: this is highly inefficient, should be optimized later
    */
    for (i=0; i<m.rows;++i){

        for (j=m.cols; j<augmented->cols; ++j){

            putElement(inv,i,j-m.cols, getElement(*augmented,i,j));
        }
    }

    /*BinMatrix* check=product(*inv,m);
    if (!compareMatrices(*check,*identityMatrix(m.rows)) ){
        puts("ERROR");
        printMatrix(*check);
        printMatrix(*inv);
        printMatrix(*augmented);
        exit(0);
    }*/

    // free memory
    destroyMatrix(augmented);

    return inv;
}

/**
 * @brief Subsample some rows or columns from a matrix
 *
 * @param indexes Array of rows(columns) to sample
 * @param len Length of the array
 * @param m Matrix to sample
 * @param mode Specify what to sample, rows or columns
 * @return BinMatrix* Matrix of the sampled rows(columns), NULL for error
 *
 * TODO: this implementation is very naive and slow. Optimize later
 */
BinMatrix* sampleFromMatrix(int* indexes, int len, BinMatrix m, int mode){

    BinMatrix* res = (BinMatrix*) (malloc(sizeof(BinMatrix)));
    int i,j;
    int needed_ulong;

    switch (mode)
    {

    case MATRIX_SAMPLE_ROWS:
        if ( len > m.rows ){
            printf("Cannot sample %d rows: matrix only has %d rows\n", len, m.rows);
            free(res);
            return NULL;
        }

        res->rows=len;
        res->cols=m.cols;

        needed_ulong = ceil( ( res->rows*res->cols *1.0) / (8*sizeof(unsigned long)) );
        res->data=(unsigned long*) malloc(sizeof(unsigned long)*needed_ulong);

        const int row_len = m.cols;

        // Iterate over all rows to extract
        for (i=0; i<len; ++i){

            for (j=0;j<row_len;++j){

                putElement(res,i,j, getElement(m,indexes[i],j));
            }
        }
        break;

    case MATRIX_SAMPLE_COLUMNS:
        if ( len > m.cols ){
            printf("Cannot sample %d columns: matrix only has %d columns\n", len, m.cols);
            free(res);
            return NULL;
        }

        res->rows=m.rows;
        res->cols=len;

        needed_ulong = ceil( ( res->rows*res->cols*1.0 ) / (8*sizeof(unsigned long)) );
        res->data=(unsigned long*) malloc(sizeof(unsigned long)*needed_ulong);

        const int col_len = m.rows;

        // Iterate over all columns to extract
        for (j=0; j<len; ++j){

            for (i=0;i<col_len;++i){

                putElement(res,i,j, getElement(m,i,indexes[j]));
            }
        }
        break;

    default:
        printf("Invalid extraction mode\n");
        return NULL;
        break;
    }

    return res;
}

/**
 * @brief Computes the hamming distance between two row vectors
 *
 * @param m1
 * @param m2
 * @return int The hamming distance
 */
int HammingDistance(BinMatrix m1, BinMatrix m2){

    if (m1.rows >1 || m2.rows >1){
        puts("Hamming distance can be computed between row vectors!!!");
        return -1;
    }

    int dist=0;
    for(int i=0;i< m1.cols;++i){
        if (getElement(m1,0,i)!=getElement(m2,0,i))
            dist++;
    }
    return dist;
}


/**
 * @brief Computes the hamming weight of a row vector
 *
 * @param m1
 * @param m2
 * @return int The hamming weight
 */
int HammingWeight(BinMatrix m1){

    if (m1.rows >1){
        puts("Hamming distance can be computed between row vectors!!!");
        return -1;
    }

    int weight=0;

    for(int i=0; i< m1.cols; ++i)
        if (getElement(m1,0,i) == 1)
            weight++;
    
    return weight;
}

/**
 * @brief Return the weight of a code
 *
 * @param v The code as a row or column vector
 * @return int The weight, -1 for error
 */
int codeWeight(BinMatrix v){

    int w=0;

    if (isRowVector(v)){

        for (int j=0; j<v.cols; ++j){

            if (getElement(v,0,j)==1)
                w++;
        }

        return w;
    }

    if (isColumnVector(v)){

        for (int i=0; i<v.rows; ++i){

            if (getElement(v,i,0)==1)
                w++;
        }
        return w;
    }

    printf("Input matrix is not a vector!\n");
    return MATRIX_INVALID_WEIGHT;

}
#endif
