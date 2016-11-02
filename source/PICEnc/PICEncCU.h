#ifndef PIC_ENC_CU_H
#define PIC_ENC_CU_H

#include "ComDef.h"

void PICEncCU(s16 *img, u16 width, u16 height, paramStruct param);

void generateCP(s16 *img, u16 width, u16 height, cpStruct *cpLargest);
void sliceFilt(cpStruct *cp, u8 preserve);

#endif