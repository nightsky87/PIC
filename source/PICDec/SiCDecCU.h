#ifndef SIC_DEC_CU_H
#define SIC_DEC_CU_H

#include "ComDef.h"
#include "SiCComTU.h"
#include "SiCComQuant.h"
#include "SiCComPred.h"
#include "SiCComBAC.h"
#include <cstring>

void SiCDecCU(s16 *img, u16 width, u16 height, paramStruct param);

void copyFromCU(s16 *img, u16 width, u16 height, cuStruct &cu);
void copyFromSCU(cuStruct &cu, scuStruct &hscu, scuStruct &vscu);

#endif
