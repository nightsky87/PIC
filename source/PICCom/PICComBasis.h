#ifndef PIC_COM_BASIS_H
#define PIC_COM_BASIS_H

#include "ComDef.h"

void basisSearch(cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 qp);
void basisInverse(cuStruct *cu, const u8 sWidth, const u8 sHeight);

#endif