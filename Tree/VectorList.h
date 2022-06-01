#ifndef VECTOR_LIST_HEADER
#define VECTOR_LIST_HEADER

#include "../Matrix/BinaryMatrix.h"

struct VectorNode{
    BinMatrix* v;
    struct VectorNode* next;
};

typedef struct VectorNode* VectorList;

/**
 * @brief Prints a VectorList
 * 
 * @param p The VectorList to print
 */
void VectorList_print(VectorList p){

    VectorList temp=p;

    while (temp!= NULL)
    {
        printMatrix(*(p->v));
        temp=temp->next;
    }
    
}

/**
 * @brief Adds a new vector to the list
 * 
 * @param m The vector to add
 * @param p Pointer to the list
 */
void VectorList_addHead(BinMatrix* m, VectorList* p){

    VectorList newNode = (VectorList) malloc(sizeof(struct VectorNode));
    newNode->v=m;

    if (*p==NULL){
        
        newNode->next=NULL;
        *p=newNode;

    }
    else{

        newNode->next=*p;
        *p=newNode;
    }
}

/**
 * @brief Destroys a vector list
 * 
 * @param l Pointer to the list to destroy
 */
void VectorList_destroy(VectorList* l){

    if (l==NULL )
        return;

    VectorList temp = *l;
    VectorList old;

    while (temp != NULL)
    {   
        destroyMatrix(temp->v);
        old=temp;
        temp=temp->next;
        free(old);
    }
    
}

/**
 * @brief Pops an item from the list head
 * 
 * @param p 
 * @return VectorList 
 */
VectorList VectorList_pop(VectorList* p){

    if (*p == NULL){
        printf("End of vectorlist reached\n");
        return NULL;
    }

    VectorList ret;
    ret = *p;
    (*p)=(*p)->next;

    VectorList_print(ret);
    return ret;
}

bool VectorList_search(BinMatrix v, VectorList l){

    VectorList temp=l;

    while (temp!=NULL)
    {
        if (compareMatrices(v,*temp->v)==0)
            return true;
    }

    return false;
    
}

#endif