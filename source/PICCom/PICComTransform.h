#ifndef PIC_COM_TRANSFORM_H
#define PIC_COM_TRANSFORM_H

#include "ComDef.h"

void dct(pel *tb, const u8 shift, const u8 qp);
void dct(const cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 shift, const u8 qp);
void dct(pel *tb, const u8 width, const u8 height, const u8 shift, const u8 qp);

void idct(pel *tb, const u8 shift, const u8 qp);
void idct(const cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 shift, const u8 qp);
void idct(pel *tb, const u8 sWidth, const u8 sHeight, const u8 qp, const u8 shift);

#endif
