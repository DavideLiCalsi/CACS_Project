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


Info *readData(char *path){

    FILE *file_ptr = fopen(path, "r");
    char ch[1000];
    int n, seed, w;

    if (NULL == file_ptr) 
        printf("file can't be opened \n");

    // read "# n"
    fgets(ch, 1000, file_ptr);
    fscanf(file_ptr, "%d\n", &n);

    // read "# seed"
    fgets(ch, 1000, file_ptr);
    fscanf(file_ptr, "%d\n", &seed);

    // read "# w"
    fgets(ch, 1000, file_ptr);
    fscanf(file_ptr, "%d\n", &w);

    // read "# H^transpose (each line corresponds to column of H, the identity part is omitted)"
    fgets(ch, 1000, file_ptr);
    
    // the dimension of the matrix is (n-k) x (n-k) where n=2k --> dimension of the array: k^2 = (n/2)^2
    int *tmp_array = (int *)malloc(sizeof(int)*n/2);
    int *H_array = (int *)malloc(sizeof(int)*(n*n/4));
    int i,j;
    for (i=0; i<n/2; i++){
        fscanf(file_ptr, "%d", &tmp_array[i]);
        for (j=n/2-1; j>=0; j--){
            H_array[(n/2)*i+j] = tmp_array[i]%10;
            tmp_array[i] /= 10;
        }
    }
    // read "\n"
    fgets(ch, 1000, file_ptr);

    // read "# s^transpose"
    fgets(ch, 1000, file_ptr);
    printf("%s", ch);
    int tmp;
    fscanf(file_ptr, "%d", &tmp);
    int *s_array = (int *)malloc(sizeof(int)*n/2);
    for (int i=n/2-1; i>=0; i--){
        s_array[i] = tmp%10;
        tmp /= 10;
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


void destroyInfo(Info *info){
    destroyMatrix(info->H_t);
    destroyMatrix(info->s);
    free(info);
}

#endif