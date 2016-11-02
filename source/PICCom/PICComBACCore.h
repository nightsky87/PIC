#ifndef PIC_COM_BAC_CORE_H
#define PIC_COM_BAC_CORE_H

#include "ComDef.h"
#include "PICComBACContexts.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>

u16 rangeTabLPS(u8 pStateIdx, u8 qCodIRangeIdx);
u8 transIdxLPS(u8 pStateIdx);

void InitializeEncoder(u32 numBytes);
void EncodeDecision(u8 ctxIdx, u8 binVal);
void EncodeBypass(u8 binVal);
void EncodeTerminate();
void PutBit(u8 B);
void WriteBitstream(char *fName, u16 width, u16 height, u8 channels, paramStruct param);

void InitializeDecoder();
void ReadBitstream(char *fName, u16 *width, u16 *height, u8 *channels, paramStruct *param);
u8 DecodeDecision(u8 ctxIdx);
u8 DecodeBypass();
u8 DecodeTerminate();
u32 ReadBits(u8 numBits);

#endif
