#ifndef PIC_COM_QUANT_H
#define PIC_COM_QUANT_H

#include "ComDef.h"

void quantConst(s32 *tb, u16 qp);
void quantConst(cpStruct *cp, u16 qp);

void dequantConst(s32 *tb, u16 qp);
void dequantConst(cpStruct *cp, u16 qp);

s16 coeffCast(double val);

#endif
