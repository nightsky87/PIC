#ifndef PIC_COM_BASIS_H
#define PIC_COM_BASIS_H

#include "ComDef.h"

void basisSearch(cuStruct *cu, const u8 sWidth, const u8 sHeight);
void basisInverse(cuStruct *cu, const u8 sWidth, const u8 sHeight);
void hadamard(const pel *img, pel *coeff, const u8 sWidth, const u8 sHeight);

#endif