#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MATRIX_SUCCESS 0
#define MATRIX_FAILURE 1
#define MATRIX_INVALID_ELEMENT 2

#define isRowVector(x) (x.rows == 1)
#define isColumnVector(x) (x.cols == 1)
#define MAX(a,b) a>=b?a:b

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

    if ( 0<=i && i<= m.rows-1 && 0<=j && j<=m.cols -1 ){

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

    if ( 0<=i && i<= m->rows-1 && 0<=j && j<=m->cols -1 ){

        int array_index= (m->cols * i + j) / (8*sizeof(unsigned long));
        unsigned long array_selector = (m->cols * i + j) % (8*sizeof(unsigned long));
        array_selector = 8*sizeof(unsigned long)- array_selector - 1;
        unsigned long selector = 1UL << array_selector;

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

    int needed_u_long = 1 + (row->rows*row->cols) / (8*sizeof(unsigned long));
    row->data= (unsigned long*) malloc(sizeof(unsigned long) * needed_u_long);
    memset(row->data,0,sizeof(unsigned long) * needed_u_long);

    for(int j=0; j<m.cols;++j){
        unsigned long new = (unsigned long) getElement(m,i,j);
        int array_offset = ceil( j / (8*sizeof(unsigned long) ) );
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

    int needed_u_long =  (col->rows*col->cols) / (8*sizeof(unsigned long));
    col->data= (unsigned long*) malloc(sizeof(unsigned long) * needed_u_long);
    memset(col->data,0,sizeof(unsigned long) * needed_u_long);

    for(int i=0; i<m.rows;++i){
        unsigned long new = (unsigned long) getElement(m,i,j);
        int array_offset = ceil( i / (8*sizeof(unsigned long) ) );
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
 * @brief Build a matrix from an array
 * @param array The array to convert
 * @param rows The rows of the matrix
 * @param cols The columns of the matrix
 * @return BinMatrix* The ptr to the new matrix object, NULL for failure
 */
BinMatrix* buildMatrix(int* array, int rows, int cols){

    BinMatrix* m = (BinMatrix*)(malloc(sizeof(BinMatrix)));

    m->rows=rows;
    m->cols=cols;
    m->data=(unsigned long*) malloc(sizeof(unsigned long)*rows*cols);

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

    BinMatrix* res = (BinMatrix*) malloc(sizeof(BinMatrix));
    res->rows=k;
    res->cols=k;
    res->data=(unsigned long*) (malloc(sizeof(unsigned long)*k*k));

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
        res->data=(unsigned long*) malloc(sizeof(unsigned long)*res->cols*res->rows);

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
        res->data=(unsigned long*) malloc(sizeof(unsigned long)*res->cols*res->rows);

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
 * @brief Computes the inner productof the two vectors
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

    int ulong_needed = 1 + (v1.cols*v1.rows) / ( 8*(sizeof(unsigned long)) );
    int j;
    char res;

    for(int i=0; i<ulong_needed;++i){

        unsigned long and = v1.data[i] & v2.data[i];

        for (j=0, res=0; j< 8*sizeof(unsigned long); ++j){

            res ^= (and >> j) & 1;
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
    res->data = (unsigned long*)(malloc(sizeof(unsigned long) * res->rows * res->cols));

    for(int i=0; i<res->rows; ++i){
        for(int j=0;j<res->cols; ++j){
            
            int val = vectorProduct( *getRow(m1,i),*getColumn(m2,j));
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
    printf("%lu\n",exctracted_row1);
    m->data[r1_array_index] = ( m->data[r1_array_index] & (-1UL ^ bitmask1) );
    m->data[r1_array_index] |= exctracted_row2;
    m->data[r2_array_index] = ( m->data[r2_array_index] & (-1UL ^ bitmask2) ) | exctracted_row1;

    return 0;
}