#include <stdio.h>
#include <stdlib.h>

struct Pair{

    int x;
    int y;

    struct  Pair* next;
    
};

typedef struct  Pair Pair;
typedef Pair* PairSet;

void PairSet_print(PairSet p){

    PairSet temp=p;

    while (temp!= NULL)
    {
        printf("Pair (%d,%d)\n",temp->x,temp->y);
        temp=temp->next;
    }
    
}

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

void PairSet_destroy(PairSet p){
    free(p);
}