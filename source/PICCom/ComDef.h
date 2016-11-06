#ifndef COMDEF_H
#define COMDEF_H

#define CU_SIZE 64

#if CU_SIZE == 16
#define NUM_BASIS_SELECTORS 5 
#elif CU_SIZE == 32
#define NUM_BASIS_SELECTORS 21 
#elif CU_SIZE == 64
#define NUM_BASIS_SELECTORS 85 
#elif CU_SIZE == 128
#define NUM_BASIS_SELECTORS 341 
#endif


typedef unsigned char	u8;
typedef unsigned short	u16;
typedef unsigned long	u32;
typedef signed char		s8;
typedef signed short	s16;
typedef signed long		s32;

typedef s32 pel;

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

struct paramStruct
{
	u16 qp;
	ChromaSub chromaSub;
};

struct cuStruct
{
	pel *pLuma;
	pel *pChroma1;
	pel *pChroma2;

	u8 *basisSelect;
	u16 basisIndex;

	u8 width;
	u8 height;
	ChromaSub chromaSub;
};

#endif