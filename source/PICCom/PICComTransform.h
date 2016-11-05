#ifndef PIC_COM_TRANSFORM_H
#define PIC_COM_TRANSFORM_H

#include "ComDef.h"
#include <stdlib.h>

void dct(cpStruct *cp);
void dct(s32 *tb);
void dctHorz(s32 *tb, u8 width, u8 height, u8 bitshift);
void dctVert(s32 *tb, u8 width, u8 height, u8 bitshift);

void idct(cpStruct *cp);
void idct(s32 *tb);
void idctHorz(s32 *tb, u8 width, u8 height, u8 bitshift);
void idctVert(s32 *tb, u8 width, u8 height, u8 bitshift);

#endif
