#include <stdio.h>
#include <stdlib.h>
#include "Matrix/BinaryMatrix.h"
#include "Utilities/utilities.h"
#include "Utilities/randomSelector.h"


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

    int length = 10;
    int *data = (int *)malloc(sizeof(int)*length);
    for (int i=0; i<length; i++)
        (data)[i] = (i+2)%5;

    Set *set = buildSet(data, length);
    printSet(set);

    quicksort(set->data,0,set->length-1);
    printSet(set);

    shuffle(set->data, set->length);
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
    printSet(subset);

    return 0;
}