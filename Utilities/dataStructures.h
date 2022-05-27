#ifndef DATA_STRUCTURES_HEADER
#define DATA_STRUCTURES_HEADER

#include <stdio.h>
#include <stdbool.h>

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

    if (set->length == 0){
        set->data = (int *)malloc(sizeof(int));
        set->data[0] = elem;
        set->length += 1;
        return;
    }

    bool alreadyPresent = false;
    for (int i=0; i<set->length && !alreadyPresent;i++)
        if ((set->data)[i] == elem)
            alreadyPresent = true;
    if (!alreadyPresent){
        set->length += 1;
        set->data = (int *)realloc(set->data, sizeof(int)*set->length);
        if (set->data == NULL)
            printf("SUUUUUUUUUUUUUUS\n");
        (set->data)[set->length-1] = elem;
    }
}

/**
 * @brief compute the first position where the input sets of equal length differ
 * 
 * @param set1 first input set
 * @param set2 second input set
 * @return int first position where the input sets differ (if they are the same return -1)
 */
int FindDiffBetweenSets(Set *set1, Set *set2){

    if (set1->length != set2->length){
        printf("The sets should have the same length!!!\n");
        return -2;
    }
    for (int i=0; i<set1->length; i++)
        if ((set1->data)[i] != (set2->data)[i])
            return i;
    
    return -1;
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

typedef struct _SetLinkedList {
   Set *set;
   struct _SetLinkedList *next;
}SetLinkedList;


/**
 * @brief add a new set (if not already present) to the linked list of sets in a sorted way
 * 
 * @param head the head of the linked list of set
 * @param newSet the new set to add
 * @return true if the set has been added
 * @return false if the set was already present
 */
bool sortedInsert_SetLinkedList(SetLinkedList **head, Set* newSet){

    // sort the input set to add
    quicksort(newSet->data, 0, newSet->length-1);

    // create the new node
    SetLinkedList *newNode = (SetLinkedList *)malloc(sizeof(SetLinkedList));
    newNode->set = newSet;

    // if the linked list is empty we can add the new node
    if (*head==NULL){
        newNode->next = NULL;
        *head = newNode;
        return true;
    }

    // find the right position where to add the new node (if not alredy present), then add it
    bool alreadyPresent = false;
    bool success = false;
    SetLinkedList *current = *head, *previous = NULL;
    while (current!=NULL && !alreadyPresent && !success){
        //find the first index where the input sets differ
        int indexDiff = FindDiffBetweenSets(current->set, newSet);
        // if the sets are equal, do not add the element and destroy it;
        if (indexDiff == -1){
            alreadyPresent = true;
            free(newNode);
            destroySet(newSet);
            continue;
        }
        // if the element at index indexDiff of the set of the node pointed by current is greater than the one of newSet, 
        // then current points to node NEXT to the new node containing newSet
        if (current->set->data[indexDiff] > newSet->data[indexDiff]){
            if (previous == NULL){
                /* current == head */
                newNode->next = current;
                *head = newNode;}
            else{
                previous->next = newNode;
                newNode->next = current;}
            success = true;
            continue;
        }
        else{ /* current->set[indexDiff] < newSet->data[indexDiff] */
            if (current->next == NULL){
                current->next = newNode;
                newNode->next = NULL;
                success=true;
            }
            else if (current->next->set->data[indexDiff] > newSet->data[indexDiff]){
                previous = current;
                newNode->next = current->next;
                previous->next = newNode;
                success = true;
            }
        }
        previous = current;
        current = current->next;
    }
    return success;
}


/**
 * @brief pretty prints of the set linked list
 * 
 * @param head head pointing to the linked list
 */
void printSetLinkedList(SetLinkedList *head){
    SetLinkedList *current = head;
    printf("BEGIN Set Linked list:\n");
    while (current != NULL){
        printSet(current->set);
        current = current->next;
    }
    printf("END\n");
}


/**
 * @brief Destroys the input set linked list, i.e. it frees all the allocated memory
 * 
 * @param head Pointer to the set linked list to free
 */
void destroy_SetLinkedList(SetLinkedList *head){
    if (head == NULL)
        return;
    SetLinkedList *current = head, *previous;
    while (current != NULL){
        previous = current;
        current = current->next;
        destroySet(previous->set);
        free(previous);
    }
}




#endif