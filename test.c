#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Matrix/BinaryMatrix.h"
#include "Utilities/utilities.h"
#include "Utilities/randomSelector.h"


#define length 4


int main(){
    
    /*BinMatrix* G;
    
    int G_vect[28]=
    {
        0, 1, 1, 0, 0, 0, 1,
        1, 1, 1, 1, 1, 1, 1,
        1, 0, 0, 0, 1, 0, 1,
        1, 1, 0, 0, 0, 1, 0
    };
    G = buildMatrix(G_vect,4,7);

    //printMatrix(*G);
    //printMatrix(*standardizeGeneratorMatrix(G));


    BinMatrix* H;
    
    int H_vect[28]=
    {
        1, 1, 1, 0, 0, 1, 0,
        0, 1, 0, 0, 1, 1, 1,
        1, 0, 0, 1, 0, 1, 1
    };
    H = buildMatrix(H_vect,3,7);

    printMatrix(*H);
    printMatrix(*standardizeParityMatrix(H));*/


    // SET + RANDOM SUBSET SELECTION testing

    /*
    int *data = (int *)malloc(sizeof(int)*length);
    for (int i=0; i<length; i++)
        (data)[i] = (i+2)%5;

    Set *set = buildSet(data, length);
    printSet(set);

    quicksort(set->data,0,set->length-1);
    printSet(set);

    addSetElem(set, 5);
    addSetElem(set, 10);
    addSetElem(set, 0);
    addSetElem(set, -1);
    addSetElem(set, -2);
    printSet(set);

    quicksort(set->data,0,set->length-1);
    printSet(set);

    Set *subset = getRandomElements(set, 7);
    printSet(subset);
    quicksort(subset->data, 0, subset->length-1);
    printSet(subset);*/


    // SET LINKED LIST testing
    int arrayset[10];
    for(int i=0; i <10; i++)
        arrayset[i] = i+1;
    Set *set = buildSet(arrayset, 10);
    SetLinkedList *list = NULL;

    for(int i=0; i<10; i++){
        while (!sortedInsert_SetLinkedList(&list, getRandomElements(set, 4, time(NULL) + 2*i)))
            /* nothing */ ;
    }
    printSetLinkedList(list);


    destroy_SetLinkedList(list);
    destroySet(set);


    return 0;
}