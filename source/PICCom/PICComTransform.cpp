#include "PICComTransform.h"
#include "PICComQuant.h"
#include <cstring>

void dct(pel *tb, const u8 shift, const u8 qp)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < 4; y++)
	{
		v0 = (s32)tb[CU_SIZE * y + 0] + (s32)tb[CU_SIZE * y + 3];
		v1 = (s32)tb[CU_SIZE * y + 1] + (s32)tb[CU_SIZE * y + 2];
		v2 = (s32)tb[CU_SIZE * y + 2] - (s32)tb[CU_SIZE * y + 1];
		v3 = (s32)tb[CU_SIZE * y + 3] - (s32)tb[CU_SIZE * y + 0];
		tb[CU_SIZE * y + 0] = (pel)((64 * (v0 + v1)) >> (shift + 1));
		tb[CU_SIZE * y + 1] = (pel)((-36 * v2 - 83 * v3) >> (shift + 1));
		tb[CU_SIZE * y + 2] = (pel)((64 * (v0 - v1)) >> (shift + 1));
		tb[CU_SIZE * y + 3] = (pel)((83 * v2 - 36 * v3) >> (shift + 1));
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 x = 0; x < 4; x++)
	{
		v0 = (s32)tb[CU_SIZE * 0 + x] + (s32)tb[CU_SIZE * 3 + x];
		v1 = (s32)tb[CU_SIZE * 1 + x] + (s32)tb[CU_SIZE * 2 + x];
		v2 = (s32)tb[CU_SIZE * 2 + x] - (s32)tb[CU_SIZE * 1 + x];
		v3 = (s32)tb[CU_SIZE * 3 + x] - (s32)tb[CU_SIZE * 0 + x];
		tb[CU_SIZE * 0 + x] = (pel)quantVal((64 * (v0 + v1)) >> 8, qp);
		tb[CU_SIZE * 1 + x] = (pel)quantVal((-36 * v2 - 83 * v3) >> 8, qp);
		tb[CU_SIZE * 2 + x] = (pel)quantVal((64 * (v0 - v1)) >> 8, qp);
		tb[CU_SIZE * 3 + x] = (pel)quantVal((83 * v2 - 36 * v3) >> 8, qp);
	}
}

void dct(const cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 shift, const u8 qp)
{
	// Define the half dimensions
	const u8 hWidth = sWidth / 2;
	const u8 hHeight = sHeight / 2;

	// Transform the luma component
	dct(&cu->pLuma[hWidth], hWidth, hHeight, shift, qp);
	dct(&cu->pLuma[CU_SIZE * hHeight], hWidth, hHeight, shift, qp);
	dct(&cu->pLuma[CU_SIZE * hHeight + hWidth], hWidth, hHeight, shift, qp);

	// Transform the chroma components
	if (cu->chromaSub != CHROMA_400)
	{
		// Transform the first chroma component
		dct(&cu->pChroma1[hWidth], hWidth, hHeight, shift, qp);
		dct(&cu->pChroma1[CU_SIZE * hHeight], hWidth, hHeight, shift, qp);
		dct(&cu->pChroma1[CU_SIZE * hHeight + hWidth], hWidth, hHeight, shift, qp);

		// Transform the second chroma component
		dct(&cu->pChroma2[hWidth], hWidth, hHeight, shift, qp);
		dct(&cu->pChroma2[CU_SIZE * hHeight], hWidth, hHeight, shift, qp);
		dct(&cu->pChroma2[CU_SIZE * hHeight + hWidth], hWidth, hHeight, shift, qp);
	}
}

void dct(pel *tb, const u8 width, const u8 height, const u8 shift, const u8 qp)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < height; y++)
	{
		for (u8 x = 0; x < width; x += 4)
		{
			v0 = (s32)tb[CU_SIZE * y + x + 0] + (s32)tb[CU_SIZE * y + x + 3];
			v1 = (s32)tb[CU_SIZE * y + x + 1] + (s32)tb[CU_SIZE * y + x + 2];
			v2 = (s32)tb[CU_SIZE * y + x + 2] - (s32)tb[CU_SIZE * y + x + 1];
			v3 = (s32)tb[CU_SIZE * y + x + 3] - (s32)tb[CU_SIZE * y + x + 0];
			tb[CU_SIZE * y + x + 0] = (pel)((64 * (v0 + v1)) >> (shift + 1));
			tb[CU_SIZE * y + x + 1] = (pel)((-36 * v2 - 83 * v3) >> (shift + 1));
			tb[CU_SIZE * y + x + 2] = (pel)((64 * (v0 - v1)) >> (shift + 1));
			tb[CU_SIZE * y + x + 3] = (pel)((83 * v2 - 36 * v3) >> (shift + 1));
		}
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 x = 0; x < width; x++)
	{
		for (u8 y = 0; y < height; y += 4)
		{
			v0 = (s32)tb[CU_SIZE * (y + 0) + x] + (s32)tb[CU_SIZE * (y + 3) + x];
			v1 = (s32)tb[CU_SIZE * (y + 1) + x] + (s32)tb[CU_SIZE * (y + 2) + x];
			v2 = (s32)tb[CU_SIZE * (y + 2) + x] - (s32)tb[CU_SIZE * (y + 1) + x];
			v3 = (s32)tb[CU_SIZE * (y + 3) + x] - (s32)tb[CU_SIZE * (y + 0) + x];
			tb[CU_SIZE * (y + 0) + x] = (pel)quantVal((64 * (v0 + v1)) >> 8, qp);
			tb[CU_SIZE * (y + 1) + x] = (pel)quantVal((-36 * v2 - 83 * v3) >> 8, qp);
			tb[CU_SIZE * (y + 2) + x] = (pel)quantVal((64 * (v0 - v1)) >> 8, qp);
			tb[CU_SIZE * (y + 3) + x] = (pel)quantVal((83 * v2 - 36 * v3) >> 8, qp);
		}
	}
}

void idct(pel *tb, const u8 shift, const u8 qp)
{
	s32 t0, t1, t2, t3;
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < 4; y++)
	{
		t0 = dequantVal(tb[CU_SIZE * y + 0], qp);
		t1 = dequantVal(tb[CU_SIZE * y + 1], qp);
		t2 = dequantVal(tb[CU_SIZE * y + 2], qp);
		t3 = dequantVal(tb[CU_SIZE * y + 3], qp);
		v0 = 64 * (t0 + t2);
		v1 = 64 * (t0 - t2);
		v2 = -36 * t1 + 83 * t3;
		v3 = -83 * t1 - 36 * t3;
		tb[CU_SIZE * y + 0] = (pel)((v0 - v3) >> 7);
		tb[CU_SIZE * y + 1] = (pel)((v1 - v2) >> 7);
		tb[CU_SIZE * y + 2] = (pel)((v1 + v2) >> 7);
		tb[CU_SIZE * y + 3] = (pel)((v0 + v3) >> 7);
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 x = 0; x < 4; x++)
	{
		v0 = 64 * ((s32)tb[CU_SIZE * 0 + x] + (s32)tb[CU_SIZE * 2 + x]);
		v1 = 64 * ((s32)tb[CU_SIZE * 0 + x] - (s32)tb[CU_SIZE * 2 + x]);
		v2 = -36 * (s32)tb[CU_SIZE * 1 + x] + 83 * (s32)tb[CU_SIZE * 3 + x];
		v3 = -83 * (s32)tb[CU_SIZE * 1 + x] - 36 * (s32)tb[CU_SIZE * 3 + x];
		tb[CU_SIZE * 0 + x] = (pel)((v0 - v3) >> (12 - shift));
		tb[CU_SIZE * 1 + x] = (pel)((v1 - v2) >> (12 - shift));
		tb[CU_SIZE * 2 + x] = (pel)((v1 + v2) >> (12 - shift));
		tb[CU_SIZE * 3 + x] = (pel)((v0 + v3) >> (12 - shift));
	}
}

void idct(const cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 shift, const u8 qp)
{
	// Define the half dimensions
	const u8 hWidth = sWidth / 2;
	const u8 hHeight = sHeight / 2;

	// Transform the luma component
	idct(&cu->pLuma[hWidth], hWidth, hHeight, shift, qp);
	idct(&cu->pLuma[CU_SIZE * hHeight], hWidth, hHeight, shift, qp);
	idct(&cu->pLuma[CU_SIZE * hHeight + hWidth], hWidth, hHeight, shift, qp);

	// Transform the chroma components
	if (cu->chromaSub != CHROMA_400)
	{
		// Transform the first chroma component
		idct(&cu->pChroma1[hWidth], hWidth, hHeight, shift, qp);
		idct(&cu->pChroma1[CU_SIZE * hHeight], hWidth, hHeight, shift, qp);
		idct(&cu->pChroma1[CU_SIZE * hHeight + hWidth], hWidth, hHeight, shift, qp);

		// Transform the second chroma component
		idct(&cu->pChroma2[hWidth], hWidth, hHeight, shift, qp);
		idct(&cu->pChroma2[CU_SIZE * hHeight], hWidth, hHeight, shift, qp);
		idct(&cu->pChroma2[CU_SIZE * hHeight + hWidth], hWidth, hHeight, shift, qp);
	}
}

void idct(pel *tb, const u8 width, const u8 height, const u8 shift, const u8 qp)
{
	s32 t0, t1, t2, t3;
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < height; y++)
	{
		for (u8 x = 0; x < width; x += 4)
		{
			t0 = dequantVal(tb[CU_SIZE * y + x + 0], qp);
			t1 = dequantVal(tb[CU_SIZE * y + x + 1], qp);
			t2 = dequantVal(tb[CU_SIZE * y + x + 2], qp);
			t3 = dequantVal(tb[CU_SIZE * y + x + 3], qp);
			v0 = 64 * (t0 + t2);
			v1 = 64 * (t0 - t2);
			v2 = -36 * t1 + 83 * t3;
			v3 = -83 * t1 - 36 * t3;
			tb[CU_SIZE * y + x + 0] = (pel)((v0 - v3) >> 7);
			tb[CU_SIZE * y + x + 1] = (pel)((v1 - v2) >> 7);
			tb[CU_SIZE * y + x + 2] = (pel)((v1 + v2) >> 7);
			tb[CU_SIZE * y + x + 3] = (pel)((v0 + v3) >> 7);
		}
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 x = 0; x < width; x++)
	{
		for (u8 y = 0; y < height; y += 4)
		{
			v0 = 64 * ((pel)tb[CU_SIZE * (y + 0) + x] + (pel)tb[CU_SIZE * (y + 2) + x]);
			v1 = 64 * ((pel)tb[CU_SIZE * (y + 0) + x] - (pel)tb[CU_SIZE * (y + 2) + x]);
			v2 = -36 * (pel)tb[CU_SIZE * (y + 1) + x] + 83 * (pel)tb[CU_SIZE * (y + 3) + x];
			v3 = -83 * (pel)tb[CU_SIZE * (y + 1) + x] - 36 * (pel)tb[CU_SIZE * (y + 3) + x];
			tb[CU_SIZE * (y + 0) + x] = (pel)((v0 - v3) >> (12 - shift));
			tb[CU_SIZE * (y + 1) + x] = (pel)((v1 - v2) >> (12 - shift));
			tb[CU_SIZE * (y + 2) + x] = (pel)((v1 + v2) >> (12 - shift));
			tb[CU_SIZE * (y + 3) + x] = (pel)((v0 + v3) >> (12 - shift));
		}
	}
}
