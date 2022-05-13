#ifndef DATA_STRUCTURES_HEADER
#define DATA_STRUCTURES_HEADER

#include <stdio.h>

/**********************************************************************************************
 *                                                                                            *
 *                                          SET                                               *
 *                                                                                            *
 *********************************************************************************************/

typedef struct _Set{
    int length;
    int* data;
}Set;


/**
 * @brief convert the array into a set (no duplicates are allowed)
 * 
 * @param array array to convernt into a set
 * @param length length of the array (number of elements)
 * @return Set* converted array into set
 */
Set* buildSet(int* array, int length){

    if (length<0){
        printf("Negative length is not allowed!!!\n");
        return NULL;
    }
    
    int setLength = 0;
    int *arraySet = (int *)malloc(sizeof(int)*length);

    for (int i=0; i<length; i++){
        bool alreadyPresent = false;
        for (int j=0; j<i && !alreadyPresent; j++)
            if (arraySet[j]==array[i])
                alreadyPresent=true;
        if (!alreadyPresent){
            arraySet[setLength] = array[i];
            setLength++;
        }
    }

    arraySet = (int*)realloc(arraySet, sizeof(int)*setLength);
    Set *set = (Set *)malloc(sizeof(Set));
    set->length = setLength;
    set->data = arraySet;
    return set;
}


/**
 * @brief Destroys the input set, i.e. it frees all the allocated memory
 * 
 * @param set Pointer to the set to free
 */
void destroySet(Set* set){

    if (set == NULL || set->data == NULL){
        printf("Error: cannot free null pointer\n");
        return;
    }

    free(set->data);
    free(set);

}



void addSetElem (Set *set, int elem){
    bool alreadyPresent = false;
    for (int i=0; i<set->length && !alreadyPresent;i++)
        if ((set->data)[i] == elem)
            alreadyPresent = true;
    if (!alreadyPresent){
        set->length += 1;
        set->data = (int *)realloc(set->data, sizeof(int)*set->length);
        (set->data)[set->length-1] = elem;
    }
}


/**
 * @brief Pretty print the input set
 * 
 * @param set input set
 */
void printSet(Set *set){

    if (set->length < 0){
        printf("A set cannot have a negative number of elements!!!\n");
        return;
    }
    printf("Set of %d elements: {", set->length);
    for (int i=0; i<set->length-1; i++)
        printf("%d, ", (set->data)[i]);
    if (set->length == 0)
        printf("}\n");
    else
        printf("%d}\n", set->data[set->length-1]);
}

/**
 * @brief quicksort algorithm implementation
 * 
 * @param array array of integers
 * @param first index of the array from where to start sorting
 * @param last index of the array where to stop sorting 
 */
void quicksort(int *array,int first,int last){

    int i, j, pivot, temp;

    if(first<last){
        pivot=first;
        i=first;
        j=last;

        while(i<j){
            while(array[i]<=array[pivot]&&i<last)
                i++;
            while(array[j]>array[pivot])
                j--;
            if(i<j){
                temp=array[i];
                array[i]=array[j];
                array[j]=temp;
            }
        }

        temp=array[pivot];
        array[pivot]=array[j];
        array[j]=temp;
        quicksort(array,first,j-1);
        quicksort(array,j+1,last);

    }
}



/**********************************************************************************************
 *                                                                                            *
 *                                  SET_LINKED_LIST                                           *
 *                                                                                            *
 *********************************************************************************************/


#endif