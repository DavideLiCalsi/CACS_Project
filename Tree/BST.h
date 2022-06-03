#include <stdio.h>
#include <stdlib.h>
#include "../Matrix/BinaryMatrix.h"
#include "VectorList.h"

#define BST_RES_ERROR_NULL_PTR 1
#define BST_RES_SUCCESS 0

/**
 * @brief Constants defining the type of key to compare
 * 
 */
#define BST_COMPARISON_INT 0
#define BST_COMPARISON_BINMATRIX 1

typedef struct _BSTNode
{
    void* key;
    void* data;
    struct _BSTNode* father;
    struct _BSTNode* l;
    struct _BSTNode* r;
    
}BSTNode;

typedef BSTNode* BST;

/**
 * @brief Compare two key stored in the node of a BST
 * 
 * @param d1 First value to compare
 * @param d2 Second value to compare
 * @param type Type of comparison
 * @return int 1 if d1>d2, -1 if d1<d2, 0 if d1=d2
 */
int compareData(void* d1, void* d2, int type){
    int* i1 = (int*) d1;
    int* i2 = (int*) d2;

    BinMatrix m1 = *(BinMatrix*) d1;
    BinMatrix m2 = *(BinMatrix*) d2;

    switch (type)
    {
    case BST_COMPARISON_INT:

        if ( *i1 > *i2)
            return 1;
        
        if ( *i1 < *i2)
            return -1;
        
        return 0;
        
        break;
    
    case BST_COMPARISON_BINMATRIX:

        return compareVectors(m1,m2);
        break;
    
    default:
        return -2;
        break;
    }

}

void destroyTree(BST* tree, int type){

    BST temp = *tree;
    BST left=temp->l;
    BST right=temp->r;

    if (left != NULL)
        destroyTree(&left, type);

    if (right != NULL)
        destroyTree(&right, type);

    switch (type)
    {
    case BST_COMPARISON_BINMATRIX:
        destroyMatrix((BinMatrix*) temp->key);
        VectorList_destroy((VectorList*) &(temp->data) );
        free(temp);
        //printf("Freed node!\n");
        break;
    
    default:
        break;
    }
}

/**
 * @brief Adds a node to a BST
 * 
 * @param node The node to add
 * @param tree The tree to whom you add the node
 * @param type The type of key
 * @return int A error or success code
 */
int addNode(void* key,void* data, BST* tree, int type){

    int* new_key;
    BinMatrix* new_key_vect;

    // First create the new node
    BST new_node = (BST) malloc(sizeof(BSTNode));

    if (type == BST_COMPARISON_INT){
        new_key = (int*) malloc(sizeof(int));
        *new_key=*(int*)key;
        new_node->key=(void*)new_key;
    }

    if (type == BST_COMPARISON_BINMATRIX){
        new_key_vect = (BinMatrix*) malloc(sizeof(BinMatrix));
        *new_key_vect = *(BinMatrix*)key;
        new_node->key=(void*)new_key_vect;
    }
    
    new_node->l=NULL;
    new_node->r=NULL;

    /*
    If the tree is empty, just set the pointer to point to
    the new node.
    */
    if (*tree == NULL){
        new_node->father=NULL;
        *tree = new_node;
        VectorList errors = NULL;
        VectorList_addHead((BinMatrix*)data,&errors);
        new_node->data=errors;

        return BST_RES_SUCCESS;
    }

    /*
    If the tree is not empty, find the right leaf node
    */
    BST tmp = *tree;

    while (1)
    {
        /* code */
        int res = compareData(key,tmp->key,type);

        switch (res)
        {
        case 1:

            if (tmp->r != NULL)
                tmp=tmp->r;
            else{
                tmp->r = new_node;
                new_node->father=tmp;
                VectorList errors = NULL;
                VectorList_addHead((BinMatrix*)data,&errors);
                new_node->data=errors;
                
                return BST_RES_SUCCESS;
            }
            break;
        
        case -1:
            if (tmp->l != NULL)
                tmp=tmp->l;
            else{
                tmp->l = new_node;
                new_node->father=tmp;
                VectorList errors = NULL;
                VectorList_addHead((BinMatrix*)data,&errors);
                new_node->data=errors;
                return BST_RES_SUCCESS;
            }
            break;
        
        default:
            
            if (type==BST_COMPARISON_BINMATRIX){
                VectorList_addHead((BinMatrix*)data, (VectorList*)&(tmp->data) );
                free(new_node);
            }
            return BST_RES_SUCCESS;
            break;
        }
    }

    return BST_RES_ERROR_NULL_PTR;
}

void printNode(BSTNode n, int type){

    puts("######## NODE ########");

    switch (type)
    {
    case BST_COMPARISON_INT:
        printf("KEY: %d\t",*(int*)(n.key));
        printf("DATA: %d\n", *(int*)(n.data));
        break;

    case BST_COMPARISON_BINMATRIX:
        printf("KEY:\n");
        printMatrix(*(BinMatrix*)(n.key));
        printf("DATA:\n");
        VectorList_print((VectorList)n.data);
        //printMatrix(*(BinMatrix*)(n.data));
    
    default:
        break;
    }

    puts("######## END ########\n");
}

/**
 * @brief Prints a BST
 * 
 * @param tree The tree to print
 */
void printTree(BST tree, int type){

    BST tmp = tree;

    if (tree->father==NULL)
        puts("######## TREE ########\n");

    if (tree->l != NULL)
        printTree(tree->l,type);

    printNode(*tree,type);
    
    if (tree->r != NULL)
        printTree(tree->r,type);
}

/**
 * @brief Search a node with a given key in the BST
 * 
 * @param key The key to search for
 * @param tree The tree to inspect
 * @return BSTNode* Pointer to node with that key if present, null if not present
 */
BSTNode* searchNode(void* key, BST tree, int type){

    if (tree == NULL)
        return NULL;

    BST tmp = tree;

    while (1)
    {
        /* code */
        int res =compareData(key,tmp->key,type);

        switch (res)
        {
        case 1:

            if (tmp->r != NULL)
                tmp=tmp->r;
            else{
                return NULL;
            }
            break;
        
        case -1:
            if (tmp->l != NULL)
                tmp=tmp->l;
            else{
                return NULL;
            }
        break;

        case 0:
            return tmp;
        
        default:
            return NULL;
            break;
        }
    }
}