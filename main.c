//#include "Matrix.h"
//#include "Matrix/BinaryMatrix.h"
//#include "Tree/BST.h"
#include "Utilities/randomSelector.h"
#include <stdbool.h>

int main(){
    
    /*BinMatrix* H;
    int H_vect[128]=
    {
        0,0,1,1,1,0,1,0,0,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,
        1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1
    };

    H = buildMatrix(H_vect,2,64);
    printMatrix(*H);

    addRows(H,0,1);

    printMatrix(*H);
    int rows[3]={0,1,3};
    BinMatrix* code = sampleFromMatrix(rows,1,*H, MATRIX_SAMPLE_ROWS);
    printf("%d\n",codeWeight(*code));*/

    int arr[100];

    for (int i=0; i<100;++i)
        arr[i]=i;
    
    Set* s = buildSet(arr,100); 
    printSet(s);
    
    return 0;
}