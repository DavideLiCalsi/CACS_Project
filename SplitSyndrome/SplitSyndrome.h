#include "../Tree/BST.h"
#include "PairSet.h"
#include <math.h>


/**
 * @brief compute the (n,k) binomial coefficient
 * 
 * @param n 
 * @param k 
 * @return int 
 */
int binomialCoeff(int n, int k){
    // Base Cases
    if (k > n)
        return 0;
    if (k == 0 || k == n)
        return 1;
 
    // Recur
    return binomialCoeff(n - 1, k - 1)
           + binomialCoeff(n - 1, k);
}

/**
 * @brief Finds the pairs (u,m) such that
 * the tables Xl and Xr have a size whose difference
 * is below threshold
 * 
 * @param n The vectors length to be considered
 * @param t The maximum number of errors
 * @param threshold The maximum difference between the sizes
 * @param E Pointer to the list the found entries
 */
void findEqualSize_u_m(int n, int t, int threshold, PairSet* E){

    int u,m;
    int size_l, size_r;

    for (m=1; m<n;++m){

        size_l = 1;
        size_r = binomialCoeff(n-m,t-0);

        for (u=0; u<=m && u<t && t-u<=n-m;++u){

            if ( abs(size_r-size_l) < threshold ){
                //printf("Found pair: u=%d, m=%d\nSize(Xl)=%d\nSize(Xl)=%d\n",u,m,size_l,size_r);
                PairSet_addHead(u,m,E);
            }
            
            size_l = size_l * (m-u)/(u+1);
            size_r = size_r * (t-u)/(n-m-t+u+1);
        }
    }

}

/**
 * @brief Converts an integer to binary vector
 * 
 * @param x The integer to convert
 * @param v The vector to populate
 * @param len The length of v
 * @returns The number of 1s in the vector
 */
int intToBinVector(int x, int* v, int len){
    int count =0;
    int res=0;

    for (int i=0; i<len; ++i){
        res = x % 2;
        v[len-i-1]= res;
        count += res;
        x = x >> 1;
    }

    return count;
}

/**
 * @brief Obtains the left syndrome from the right syndrome
 * and the full syndrome
 * 
 * @param s 
 * @param sr 
 * @return BinMatrix* 
 */
BinMatrix* getLeftSyndrome(BinMatrix s, BinMatrix sr){

    return vectorSum(s,sr);
}

/**
 * @brief Implements step 4 of the split-syndrome algorithm
 * 
 * @param Xr The right table
 * @param Xl The left table
 * @param s The full syndrome
 * @param el Contains the left error at the end of the algortithm
 * @param er Contains the right error at the end of the algortithm
 * @return true 
 * @return false 
 */
bool inspectTables(BST Xr, BST Xl, BinMatrix s, BinMatrix** el, BinMatrix** er){

    BSTNode* node;
    bool found=false;

    if (Xr->l != NULL){
        found = inspectTables(Xr->l,Xl,s, el,er);

        if (found)
            return found;
    }

    BinMatrix* left_syndrome = getLeftSyndrome(s, *(BinMatrix*) Xr->key);
    node = searchNode((void*)left_syndrome,Xl,BST_COMPARISON_BINMATRIX);

    if (node != NULL){
        *er = (BinMatrix*) Xr->data;
        *el = (BinMatrix*) node->data;
        return true;
    }

    if (Xr->r != NULL){
        found = inspectTables(Xr->r,Xl,s,el,er);

        return found;
    }
    else
        return false;
}

/**
 * @brief Builds the left table Xl
 * 
 * @param m The length of the left side
 * @param u The number of errors considered
 * @param H The parity check matrix
 * @param Xl The left table to populate
 */
void buildLeftTable(int m,int u, BinMatrix H, BST* Xl){

    int i;
    int indexes[m];
    int error_as_vector[m];

    for (i=0; i<m;++i){
        indexes[i]=i;
        error_as_vector[i]=0;
    }

    BinMatrix* Hl = transpose(*sampleFromMatrix(indexes,m,H,MATRIX_SAMPLE_COLUMNS));

    for (i=0; i<pow(2,m);++i){

        if ( intToBinVector(i,error_as_vector,m) != u )
            continue;
        BinMatrix* el = buildMatrix(error_as_vector,1,m);
        BinMatrix* sl = product(*el,*Hl);
        addNode((void*)sl, (void*) el, Xl,BST_COMPARISON_BINMATRIX);
    }

}

/**
 * @brief Builds the right table Xr
 * 
 * @param m The length of the left part
 * @param t_minus_u The number of errors in the right part, i.e. t-u
 * @param len The total syndrome length
 * @param H Parity check matrix
 * @param Xr The table to populate
 */
void buildRightTable(int m,int t_minus_u, int len, BinMatrix H, BST* Xr){

    int i;
    int indexes[len-m];
    int error_as_vector[len-m];

    for (i=0; i<len-m;++i){
        indexes[i]=i+m;
        error_as_vector[i]=0;
    }

    BinMatrix* Hr = transpose(*sampleFromMatrix(indexes,len-m,H,MATRIX_SAMPLE_COLUMNS));

    for (i=0; i<pow(2,len-m);++i){

        if ( intToBinVector(i,error_as_vector,len-m) != t_minus_u)
            continue;
        
        BinMatrix* er = buildMatrix(error_as_vector,1,len-m);
        BinMatrix* sr = product(*er,*Hr);
        addNode((void*)sr, (void*) er, Xr,BST_COMPARISON_BINMATRIX);
    }

}

/**
 * @brief Implements the split-syndrome algorithm
 * 
 * @param H
 * @param s 
 * @param d 
 */
void SplitSyndrome(BinMatrix H, BinMatrix s, int d, BinMatrix** e){

    int t,u,m;
    BST Xl=NULL;
    BST Xr=NULL;
    BinMatrix* el;
    BinMatrix* er;
    PairSet tables[d];

    // First do the precomputation step
    printf("Starting the precomputation stage\n");

    for (t=0; t<d;++t){
        tables[t]=NULL;
        findEqualSize_u_m(s.cols,t,15,&tables[t]);
        printf("Found E(%d)\n",t);
    }

    printf("Precomputation stage complete!\n");

    
    for (t=1; t<d;++t){

        PairSet E = tables[t];

        while (E != NULL)
        {   
            Pair* pair=PairSet_pop(&E);
            u=pair->x;
            m=pair->y;

            printf("Trying the values u=%d,m=%d,t=%d\n",u,m,t);
            // Build the two tables
            printf("Building the table X_left...");
            buildLeftTable(m,u,H,&Xl);
            printf("DONE!\nBuilding the table X_right...");
            buildRightTable(m,t-u,s.cols,H,&Xr);
            printf("DONE!\n");

            if ( inspectTables(Xr,Xl,s,&el,&er) ){

                // If you found the two errors, build the full error and return it
                printf("Found code!\n");
                *e = concat(*el,*er,0);
                printMatrix(**e);
                return;
            }
        }
        
    }
}