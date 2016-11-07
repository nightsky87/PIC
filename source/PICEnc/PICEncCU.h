#ifndef PIC_ENC_CU_H
#define PIC_ENC_CU_H

#include "ComDef.h"

void PICEncCU(s16 *img, u16 width, u16 height, paramStruct param);
void reduce(const cuStruct *cu, pel *pLuma, pel *pChroma1, pel *pChroma2, const u8 sWidth, const u8 sHeight);
void filter(cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 preserve);

#endif