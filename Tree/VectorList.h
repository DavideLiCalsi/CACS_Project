#ifndef VECTOR_LIST_HEADER
#define VECTOR_LIST_HEADER

#include "../Matrix/BinaryMatrix.h"

typedef struct _VectorNode{
    BinMatrix* v;
    struct _VectorNode* next;
}VectorNode;

typedef VectorNode* VectorList;

/**
 * @brief Prints a VectorList
 * 
 * @param p The VectorList to print
 */
void VectorList_print(VectorList p){

    VectorList curr=p;

    while (curr!= NULL)
    {
        printMatrix(*(curr->v));
        curr=curr->next;
    }
    
}

/**
 * @brief Adds a new vector to the list
 * 
 * @param m The vector to add
 * @param p Pointer to the pointer to the list
 */
void VectorList_addHead(BinMatrix* m, VectorList* p){

    VectorList newNode = (VectorList) malloc(sizeof(VectorNode));
    newNode->v = m;

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
 * @brief Destroy the node of a VectorList
 * 
 * @param l Pointer to the node to destroy
 */
void VectorList_nodeDestroy(VectorList v){

    destroyMatrix(v->v);
    free(v);
}

/**
 * @brief Destroys a whole vector list
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
VectorList vectorList_pop(VectorList* p){

    if (*p == NULL){
        printf("End of vectorlist reached\n");
        return NULL;
    }

    VectorList ret = (VectorList)malloc(sizeof(VectorNode));
    ret->v = copyMatrix((*p)->v);
    ret->next = NULL;
    VectorList old = *p;
    (*p)=(*p)->next;
    //free(old);

    //VectorList_print(ret);
    return ret;
}


/**
 * @brief get the element at position "index"
 * 
 * @param p input list
 * @param index position of the element to take
 * @return VectorList the element at position "index"
 */
VectorList vectorList_get(VectorList* p, int index){

    if (*p == NULL){
        //printf("The list is empty!\n");
        return NULL;
    }

    VectorList curr = *p;
    int count = 0;
    while (curr!=NULL){
        if (count == index){
            VectorList ret = (VectorList)malloc(sizeof(VectorNode));
            ret->v = copyMatrix(curr->v);
            ret->next = NULL;
            return ret;
        }
        count++;
        curr = curr->next;
    }

    //printf("The index %d is larger than the list's length!\n", index);
    return NULL;

}

/**
 * @brief Searches the binary vector v in a list
 * 
 * @param v The vector to search
 * @param l The list to inspect
 * @return true Vector found
 * @return false Not found
 */
bool VectorList_search(BinMatrix v, VectorList l){

    VectorList temp=l;

    while (temp!=NULL)
    {
        if (compareMatrices(v,*temp->v)==true)
            return true;

        temp=temp->next;
    }

    return false;
    
}

#endif