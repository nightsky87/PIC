#ifndef COMDEF_H
#define COMDEF_H

#define CU_SIZE 64
#define PRED_MODE_NN 0

typedef unsigned char	u8;
typedef unsigned short	u16;
typedef unsigned long	u32;
typedef signed char		s8;
typedef signed short	s16;
typedef signed long		s32;

enum Component
{
	COMPONENT_LUMA,
	COMPONENT_CHROMA
};

enum ChromaSub
{
	CHROMA_444,
	CHROMA_420,
	CHROMA_400
};

enum ScanDir
{
	SCAN_VERT,
	SCAN_HORZ,
	SCAN_DIAG
};

enum PredType
{
	TYPE_UNI,
	TYPE_LINEAR,
	TYPE_BILINEAR
};

struct paramStruct
{
	u16 qp;
	ChromaSub chromaSub;
};

struct cpStruct
{
	s32 *pLuma;
	s32 *pChroma1;
	s32 *pChroma2;

	u8 *predMode;

	u8 width;
	u8 height;
	u8 bitshift;
	ChromaSub chromaSub;

	cpStruct *larger;
	cpStruct *smaller;
};

#endif