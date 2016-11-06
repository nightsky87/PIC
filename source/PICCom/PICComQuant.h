#ifndef PIC_COM_QUANT_H
#define PIC_COM_QUANT_H

#include "ComDef.h"

s32 quantVal(const s32 val, const u8 qp);
s32 dequantVal(const s32 val, const u8 qp);

#endif
