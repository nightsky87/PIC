#include "PICComTransform.h"
#include "PICComQuant.h"
#include <cstring>

void dct(pel *tb, const u8 qp)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < 4; y++)
	{
		v0 = (s32)tb[CU_SIZE * y + 0] + (s32)tb[CU_SIZE * y + 3];
		v1 = (s32)tb[CU_SIZE * y + 1] + (s32)tb[CU_SIZE * y + 2];
		v2 = (s32)tb[CU_SIZE * y + 2] - (s32)tb[CU_SIZE * y + 1];
		v3 = (s32)tb[CU_SIZE * y + 3] - (s32)tb[CU_SIZE * y + 0];
		tb[CU_SIZE * y + 0] = (pel)((64 * (v0 + v1)) >> 1);
		tb[CU_SIZE * y + 1] = (pel)((-36 * v2 - 83 * v3) >> 1);
		tb[CU_SIZE * y + 2] = (pel)((64 * (v0 - v1)) >> 1);
		tb[CU_SIZE * y + 3] = (pel)((83 * v2 - 36 * v3) >> 1);
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

void dct(const cuStruct *cu, const u8 sWidth, const u8 shift)
{
	// Assign the height of all CU scales to be equal to the width
	const u8 sHeight = sWidth;

	// Allocate space for the transform coefficients
	static pel *buf = new pel[CU_SIZE * CU_SIZE];

	// Transform the luma component
	dct(cu->pLuma, buf, sWidth, sHeight, shift);
	dct(&cu->pLuma[1], &buf[sWidth / 2], sWidth, sHeight, shift);
	dct(&cu->pLuma[CU_SIZE], &buf[CU_SIZE * sHeight / 2], sWidth, sHeight, shift);
	dct(&cu->pLuma[CU_SIZE + 1], &buf[CU_SIZE * sHeight / 2 + sWidth / 2], sWidth, sHeight, shift);

	// Copy the resulting coefficients
	for (u8 y = 0; y < sHeight; y++)
		memcpy(&cu->pLuma[CU_SIZE * y], &buf[CU_SIZE * y], sWidth * sizeof(pel));

	// Transform the chroma components
	if (cu->chromaSub != CHROMA_400)
	{
		// Transform the first chroma component
		dct(cu->pChroma1, buf, sWidth, sHeight, shift);
		dct(&cu->pChroma1[1], &buf[sWidth / 2], sWidth, sHeight, shift);
		dct(&cu->pChroma1[CU_SIZE], &buf[CU_SIZE * sHeight / 2], sWidth, sHeight, shift);
		dct(&cu->pChroma1[CU_SIZE + 1], &buf[CU_SIZE * sHeight / 2 + sWidth / 2], sWidth, sHeight, shift);

		// Copy the resulting coefficients
		for (u8 y = 0; y < sHeight; y++)
			memcpy(&cu->pChroma1[CU_SIZE * y], &buf[CU_SIZE * y], sWidth * sizeof(pel));

		// Transform the second chroma component
		dct(cu->pChroma2, buf, sWidth, sHeight, shift);
		dct(&cu->pChroma2[1], &buf[sWidth / 2], sWidth, sHeight, shift);
		dct(&cu->pChroma2[CU_SIZE], &buf[CU_SIZE * sHeight / 2], sWidth, sHeight, shift);
		dct(&cu->pChroma2[CU_SIZE + 1], &buf[CU_SIZE * sHeight / 2 + sWidth / 2], sWidth, sHeight, shift);

		// Copy the resulting coefficients
		for (u8 y = 0; y < sHeight; y++)
			memcpy(&cu->pChroma2[CU_SIZE * y], &buf[CU_SIZE * y], sWidth * sizeof(pel));
	}
}

void dct(const pel *tb, pel *buf, const u8 sWidth, const u8 sHeight, const u8 shift)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < sHeight; y += 2)
	{
		for (u8 x = 0; x < sWidth; x += 8)
		{
			v0 = (s32)tb[CU_SIZE * y + x + 0] + (s32)tb[CU_SIZE * y + x + 6];
			v1 = (s32)tb[CU_SIZE * y + x + 2] + (s32)tb[CU_SIZE * y + x + 4];
			v2 = (s32)tb[CU_SIZE * y + x + 4] - (s32)tb[CU_SIZE * y + x + 2];
			v3 = (s32)tb[CU_SIZE * y + x + 6] - (s32)tb[CU_SIZE * y + x + 0];
			buf[CU_SIZE * y / 2 + x / 2 + 0] = (pel)((64 * (v0 + v1)) >> (shift + 1));
			buf[CU_SIZE * y / 2 + x / 2 + 1] = (pel)((-36 * v2 - 83 * v3) >> (shift + 1));
			buf[CU_SIZE * y / 2 + x / 2 + 2] = (pel)((64 * (v0 - v1)) >> (shift + 1));
			buf[CU_SIZE * y / 2 + x / 2 + 3] = (pel)((83 * v2 - 36 * v3) >> (shift + 1));
		}
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 x = 0; x < sWidth / 2; x++)
	{
		for (u8 y = 0; y < sHeight / 2; y += 4)
		{
			v0 = (s32)buf[CU_SIZE * (y + 0) + x] + (s32)buf[CU_SIZE * (y + 3) + x];
			v1 = (s32)buf[CU_SIZE * (y + 1) + x] + (s32)buf[CU_SIZE * (y + 2) + x];
			v2 = (s32)buf[CU_SIZE * (y + 2) + x] - (s32)buf[CU_SIZE * (y + 1) + x];
			v3 = (s32)buf[CU_SIZE * (y + 3) + x] - (s32)buf[CU_SIZE * (y + 0) + x];
			buf[CU_SIZE * (y + 0) + x] = (pel)((64 * (v0 + v1)) >> 8);
			buf[CU_SIZE * (y + 1) + x] = (pel)((-36 * v2 - 83 * v3) >> 8);
			buf[CU_SIZE * (y + 2) + x] = (pel)((64 * (v0 - v1)) >> 8);
			buf[CU_SIZE * (y + 3) + x] = (pel)((83 * v2 - 36 * v3) >> 8);
		}
	}
}

void idct(pel *tb, const u8 qp)
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
		tb[CU_SIZE * 0 + x] = (pel)((v0 - v3) >> 12);
		tb[CU_SIZE * 1 + x] = (pel)((v1 - v2) >> 12);
		tb[CU_SIZE * 2 + x] = (pel)((v1 + v2) >> 12);
		tb[CU_SIZE * 3 + x] = (pel)((v0 + v3) >> 12);
	}
}

void idct(const cuStruct *cu, const u8 sWidth, const u8 qp, const u8 shift)
{
	// Assign the height of all CU scales to be equal to the width
	const u8 sHeight = sWidth;

	// Transform the luma component
	idct(&cu->pLuma[sWidth / 2], sWidth, sHeight, qp, shift);
	idct(&cu->pLuma[CU_SIZE * sHeight / 2], sWidth, sHeight, qp, shift);
	idct(&cu->pLuma[CU_SIZE * sHeight / 2 + sWidth / 2], sWidth, sHeight, qp, shift);

	// Transform the chroma components
	if (cu->chromaSub != CHROMA_400)
	{
		// Transform the first chroma component
		idct(&cu->pChroma1[sWidth / 2], sWidth, sHeight, qp, shift);
		idct(&cu->pChroma1[CU_SIZE * sHeight / 2], sWidth, sHeight, qp, shift);
		idct(&cu->pChroma1[CU_SIZE * sHeight / 2 + sWidth / 2], sWidth, sHeight, qp, shift);

		// Transform the second chroma component
		idct(&cu->pChroma2[sWidth / 2], sWidth, sHeight, qp, shift);
		idct(&cu->pChroma2[CU_SIZE * sHeight / 2], sWidth, sHeight, qp, shift);
		idct(&cu->pChroma2[CU_SIZE * sHeight / 2 + sWidth / 2], sWidth, sHeight, qp, shift);
	}
}

void idct(pel *tb, const u8 sWidth, const u8 sHeight, const u8 qp, const u8 shift)
{
	s32 t0, t1, t2, t3;
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < sHeight / 2; y++)
	{
		for (u8 x = 0; x < sWidth / 2; x += 4)
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
	for (u8 x = 0; x < sWidth / 2; x++)
	{
		for (u8 y = 0; y < sHeight / 2; y += 4)
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
