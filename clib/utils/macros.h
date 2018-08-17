
#ifndef  SAFEDEL
#define SAFEDEL(x) {if (x != NULL) {delete x; x=NULL;}}
#endif