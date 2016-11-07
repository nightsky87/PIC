#ifndef PIC_COM_BASIS_H
#define PIC_COM_BASIS_H

#include "ComDef.h"

void basisSearch(cuStruct *cu, const u8 sWidth, const u8 sHeight);
void basisInverse(cuStruct *cu, const u8 sWidth, const u8 sHeight);
u32 hadamardMetric(pel *tb, const u8 stride);

#endif