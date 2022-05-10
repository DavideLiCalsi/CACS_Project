#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MATRIX_SUCCESS 0
#define MATRIX_FAILURE 1

#define isRowVector(x) (x.rows == 1)
#define isColumnVector(x) (x.cols == 1)

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
        printf("Array index %d\nArray selector %d\n",array_index,array_selector);
        unsigned long row = m.data[array_index];
        printf("%ld\n",row);
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

    int needed_u_long =  (row->rows*row->cols) / (8*sizeof(unsigned long));
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