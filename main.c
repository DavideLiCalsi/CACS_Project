//#include "Matrix.h"
#include "BinaryMatrix.h"

int main(){
    
    BinMatrix* m = (BinMatrix*) malloc(sizeof(BinMatrix));
    m->rows=8;
    m->cols=16;
    m->data=(unsigned long*) malloc(sizeof(unsigned long)*2);
    m->data[0]=0x8000000000000000;
    m->data[1]=0xffffffffffffffff;

   long a = 0x8000ffffL;
    printf("%d\n",getElement(*m,7,15));

    putElement(m,7,15,1);

    printf("%d\n",getElement(*m,7,15));
    printf("%ld\n",getRow(*m,1)->data[0]);

    return 0;
}