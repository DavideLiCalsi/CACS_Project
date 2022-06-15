/**
 * @brief This file implements data structures to handle
 * pairs of integers and lists of pairsof integers.
 */
#include <stdio.h>
#include <stdlib.h>

struct Pair{

    int x;
    int y;

    struct  Pair* next;
    
};

typedef struct  Pair Pair;
typedef Pair* PairSet;

/**
 * @brief Pretty prints a list of pairs.
 * 
 * @param p 
 */
void PairSet_print(PairSet p){

    PairSet temp=p;

    while (temp!= NULL)
    {
        printf("Pair (%d,%d)\n",temp->x,temp->y);
        temp=temp->next;
    }
    
}

/**
 * @brief Add a pair to the head of the list
 * 
 * @param x 1st element of the pair
 * @param y 2nd elements of the pair
 * @param p The pairlist
 */
void PairSet_addHead(int x, int y, PairSet* p){

    PairSet newPair = (PairSet) malloc(sizeof(Pair));
    newPair->x=x;
    newPair->y=y;

    if (*p==NULL){
        
        newPair->next=NULL;
        *p=newPair;

    }
    else{

        newPair->next=*p;
        *p=newPair;
    }
}

/**
 * @brief Pops a pair from the head of the list
 * 
 * @param p The list
 * @return Pair* The former list head
 */
Pair* PairSet_pop(PairSet* p){

    if (*p == NULL){
        printf("End of pair reached\n");
        return NULL;
    }

    Pair* ret;
    ret = *p;
    (*p)=(*p)->next;
    return ret;
}

/**
 * @brief Destroys a pair
 * 
 * @param p The pair to destroy
 */
void PairSet_destroy(PairSet p){
    free(p);
}