#ifndef RANDOM_SELECTOR_HEADER
#define RANDOM_SELECTOR_HEADER

#include <time.h>
#include <stdlib.h>
#include "dataStructures.h"


/**
 * @brief shuffle the input array
 * 
 * @param array input array of integers
 * @param length length of the input array
 */
void shuffle(int *array, int length){

    srand(time(NULL));
    int val, j;
    for(int i=0; i<length; i++){
        j = rand() % (length - 0) + 0;  // i <= j < n
        int val = array[i];
        array[i] = array[j];
        array[j] = val;
    }
}


/**
 * @brief Get "amount" distinct random numbers from inputSet
 * 
 * @param inputset input set from where to take elements randomly
 * @param amount total number of elements to select
 * @return Set* set containing "amount" distinct random numbers belonging to inputSet
 */

Set* getDistinctRandomNumbers(Set *inputset, int amount){

    if (amount > inputset->length)
        amount = inputset->length;
    Set *set = (Set*)malloc(sizeof(Set));
    set->length=0;
    srand(time(NULL));
    while (set->length != amount){
        int randomIndex = rand() % inputset->length;
        addSetElem(set, (inputset->data)[randomIndex]);
    }
    return set;
}

/**
 * @brief Get a random subset of length "amount" from set inputSet 
 * 
 * @param inputSet input set
 * @param amount number of elements of the subset
 * @return Set* subset of the inputset
 */
Set * getRandomElements(Set *inputSet, int amount){
    /**
     * idea: if amount is near the total number of elements in the inputSet, it may be inefficient to pick randomly
     *       elements of the inputSet according to function getDistinctRandomNumbers. Thus, in case "amount" is not
     *       "small", we shuffle the copy of the elements of the inputSet and take the first "amount" elements.
     *       Empirically, it has been shown that  0.63212055882 * lengthInputSet is a good estimate whether to
     *       consider "amount" small or not.
     * Credit to: https://www.codementor.io/@alexanderswilliams/how-to-efficiently-generate-a-random-subset-150hbz3na4
     */
    int inputSetLen = inputSet->length;
    if (amount > inputSetLen)
        amount = inputSetLen;
    
    Set *result;
    if (amount < 0.63212055882 * inputSetLen)
        result = getDistinctRandomNumbers(inputSet, amount);
    else{
        int *setNumbers = (int *)malloc(sizeof(int)*inputSet->length);
        for (int i=0; i<inputSet->length; i++)
            setNumbers[i] = (inputSet->data)[i];
        shuffle(setNumbers, inputSet->length);
        result = (Set *)malloc(sizeof(Set));
        result->length = 0;
        for (int i=0; i<amount; i++)
            addSetElem(result, setNumbers[i]);
        free (setNumbers);
    }
    return result;
}
    
    

#endif