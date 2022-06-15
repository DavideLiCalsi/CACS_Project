#ifndef DATA_READER_HEADER
#define DATA_READER_HEADER

#include <stdlib.h>
#include "../Matrix/BinaryMatrix.h"

typedef struct _Info{
    int n;          // code length
    int seed;       // seed
    int w;          // minimum weight
    BinMatrix *H_t; // transpose of the parity check matrix
    BinMatrix *s;   // syndrome
}Info;


/**
 * @brief Read from file the needed information
 * 
 * @param path path to the file to read
 * @return Info* structure containing all the read information
 */
Info *readData(char *path){

    FILE *file_ptr = fopen(path, "r");
    char comment[100];
    char ch;
    int n, seed, w;

    if (NULL == file_ptr) 
        printf("file can't be opened \n");

    // read "# n"
    fgets(comment, 100, file_ptr);
    fscanf(file_ptr, "%d\n", &n);

    // read "# seed"
    fgets(comment, 100, file_ptr);
    fscanf(file_ptr, "%d\n", &seed);

    // read "# w"
    fgets(comment, 100, file_ptr);
    fscanf(file_ptr, "%d\n", &w);

    // read "# H^transpose (each line corresponds to column of H, the identity part is omitted)"
    fgets(comment, 100, file_ptr);
    
    // the dimension of the matrix is (n-k) x (n-k) where n=2k --> dimension of the array: k^2 = (n/2)^2
    unsigned long *tmp_array = (unsigned long *)malloc(sizeof(unsigned long)*n/2);
    int *H_array = (int *)malloc(sizeof(int)*(n*n/4));
    int i,j;
    for (i=0; i<n/2; i++){
        for (j=0; j<n/2; j++){
            fscanf(file_ptr, "%c", &ch);
            H_array[(n/2)*i+j] = (int)(ch - '0');
        }
        fgets(comment, 100, file_ptr);    //read "\n"
    }

    // read "# s^transpose"
    fgets(comment, 100, file_ptr);

    int *s_array = (int *)malloc(sizeof(int)*n/2);
    for (int i=0; i<n/2; i++){
        fscanf(file_ptr, "%c", &ch);
        s_array[i] = (int)(ch - '0');
    }

    Info *info = (Info *)malloc(sizeof(Info));
    info->n = n;
    info->seed = seed;
    info->w = w;
    info->H_t = buildMatrix(H_array, n/2, n/2);
    info->s = buildMatrix(s_array, 1, n/2);

    free(tmp_array);
    free(H_array);
    free(s_array);
    fclose(file_ptr);


    return info;
}

/**
 * @brief Destroys the input Info, i.e. it frees all the allocated memory
 * 
 * @param info Pointer to the Info to free
 */
void destroyInfo(Info *info){
    destroyMatrix(info->H_t);
    destroyMatrix(info->s);
    free(info);
}

#endif