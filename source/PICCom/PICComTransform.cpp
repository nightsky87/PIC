#include "PICComTransform.h"

void dct(cpStruct *cp)
{
	// Transform the luma component
	dct(cp->pLuma, cp->width, cp->height, cp->bitshift);

	// Transform the chroma components
	if (cp->chromaSub != CHROMA_400)
	{
		dct(cp->pChroma1, cp->width, cp->height, cp->bitshift);
		dct(cp->pChroma2, cp->width, cp->height, cp->bitshift);
	}
}

void dct(s32 *tb)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < 4; y++)
	{
		v0 = (s32)tb[y * 4 + 0] + (s32)tb[y * 4 + 3];
		v1 = (s32)tb[y * 4 + 1] + (s32)tb[y * 4 + 2];
		v2 = (s32)tb[y * 4 + 2] - (s32)tb[y * 4 + 1];
		v3 = (s32)tb[y * 4 + 3] - (s32)tb[y * 4 + 0];
		tb[y * 4 + 0] = (s32)((64 * (v0 + v1)) >> 9);
		tb[y * 4 + 1] = (s32)((-36 * v2 - 83 * v3) >> 9);
		tb[y * 4 + 2] = (s32)((64 * (v0 - v1)) >> 9);
		tb[y * 4 + 3] = (s32)((83 * v2 - 36 * v3) >> 9);
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 x = 0; x < 4; x ++)
	{
		v0 = (s32)tb[0 + x] + (s32)tb[12 + x];
		v1 = (s32)tb[4 + x] + (s32)tb[8 + x];
		v2 = (s32)tb[8 + x] - (s32)tb[4 + x];
		v3 = (s32)tb[12 + x] - (s32)tb[0 + x];
		tb[0 + x] = (s32)((64 * (v0 + v1)) >> 8);
		tb[4 + x] = (s32)((-36 * v2 - 83 * v3) >> 8);
		tb[8 + x] = (s32)((64 * (v0 - v1)) >> 8);
		tb[12 + x] = (s32)((83 * v2 - 36 * v3) >> 8);
	}
}

void dct(s32 *tb, u8 width, u8 height, u8 bitshift)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < height; y += 2)
	{
		for (u8 x = 0; x < width; x += 8)
		{
			v0 = (s32)tb[y * width + x + 1] + (s32)tb[y * width + x + 7];
			v1 = (s32)tb[y * width + x + 3] + (s32)tb[y * width + x + 5];
			v2 = (s32)tb[y * width + x + 5] - (s32)tb[y * width + x + 3];
			v3 = (s32)tb[y * width + x + 7] - (s32)tb[y * width + x + 1];
			tb[y * width + x + 1] = (s32)((64 * (v0 + v1)) >> (bitshift + 1));
			tb[y * width + x + 3] = (s32)((-36 * v2 - 83 * v3) >> (bitshift + 1));
			tb[y * width + x + 5] = (s32)((64 * (v0 - v1)) >> (bitshift + 1));
			tb[y * width + x + 7] = (s32)((83 * v2 - 36 * v3) >> (bitshift + 1));

			v0 = (s32)tb[(y + 1) * width + x + 0] + (s32)tb[(y + 1) * width + x + 6];
			v1 = (s32)tb[(y + 1) * width + x + 2] + (s32)tb[(y + 1) * width + x + 4];
			v2 = (s32)tb[(y + 1) * width + x + 4] - (s32)tb[(y + 1) * width + x + 2];
			v3 = (s32)tb[(y + 1) * width + x + 6] - (s32)tb[(y + 1) * width + x + 0];
			tb[(y + 1) * width + x + 0] = (s32)((64 * (v0 + v1)) >> (bitshift + 1));
			tb[(y + 1) * width + x + 2] = (s32)((-36 * v2 - 83 * v3) >> (bitshift + 1));
			tb[(y + 1) * width + x + 4] = (s32)((64 * (v0 - v1)) >> (bitshift + 1));
			tb[(y + 1) * width + x + 6] = (s32)((83 * v2 - 36 * v3) >> (bitshift + 1));

			v0 = (s32)tb[(y + 1) * width + x + 1] + (s32)tb[(y + 1) * width + x + 7];
			v1 = (s32)tb[(y + 1) * width + x + 3] + (s32)tb[(y + 1) * width + x + 5];
			v2 = (s32)tb[(y + 1) * width + x + 5] - (s32)tb[(y + 1) * width + x + 3];
			v3 = (s32)tb[(y + 1) * width + x + 7] - (s32)tb[(y + 1) * width + x + 1];
			tb[(y + 1) * width + x + 1] = (s32)((64 * (v0 + v1)) >> (bitshift + 1));
			tb[(y + 1) * width + x + 3] = (s32)((-36 * v2 - 83 * v3) >> (bitshift + 1));
			tb[(y + 1) * width + x + 5] = (s32)((64 * (v0 - v1)) >> (bitshift + 1));
			tb[(y + 1) * width + x + 7] = (s32)((83 * v2 - 36 * v3) >> (bitshift + 1));
		}
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 y = 0; y < height; y += 8)
	{
		for (u8 x = 0; x < width; x += 2)
		{
			v0 = (s32)tb[(y + 0) * width + x + 1] + (s32)tb[(y + 6) * width + x + 1];
			v1 = (s32)tb[(y + 2) * width + x + 1] + (s32)tb[(y + 4) * width + x + 1];
			v2 = (s32)tb[(y + 4) * width + x + 1] - (s32)tb[(y + 2) * width + x + 1];
			v3 = (s32)tb[(y + 6) * width + x + 1] - (s32)tb[(y + 0) * width + x + 1];
			tb[(y + 0) * width + x + 1] = (s32)((64 * (v0 + v1)) >> 8);
			tb[(y + 2) * width + x + 1] = (s32)((-36 * v2 - 83 * v3) >> 8);
			tb[(y + 4) * width + x + 1] = (s32)((64 * (v0 - v1)) >> 8);
			tb[(y + 6) * width + x + 1] = (s32)((83 * v2 - 36 * v3) >> 8);

			v0 = (s32)tb[(y + 1) * width + x] + (s32)tb[(y + 7) * width + x];
			v1 = (s32)tb[(y + 3) * width + x] + (s32)tb[(y + 5) * width + x];
			v2 = (s32)tb[(y + 5) * width + x] - (s32)tb[(y + 3) * width + x];
			v3 = (s32)tb[(y + 7) * width + x] - (s32)tb[(y + 1) * width + x];
			tb[(y + 1) * width + x] = (s32)((64 * (v0 + v1)) >> 8);
			tb[(y + 3) * width + x] = (s32)((-36 * v2 - 83 * v3) >> 8);
			tb[(y + 5) * width + x] = (s32)((64 * (v0 - v1)) >> 8);
			tb[(y + 7) * width + x] = (s32)((83 * v2 - 36 * v3) >> 8);

			v0 = (s32)tb[(y + 1) * width + x + 1] + (s32)tb[(y + 7) * width + x + 1];
			v1 = (s32)tb[(y + 3) * width + x + 1] + (s32)tb[(y + 5) * width + x + 1];
			v2 = (s32)tb[(y + 5) * width + x + 1] - (s32)tb[(y + 3) * width + x + 1];
			v3 = (s32)tb[(y + 7) * width + x + 1] - (s32)tb[(y + 1) * width + x + 1];
			tb[(y + 1) * width + x + 1] = (s32)((64 * (v0 + v1)) >> 8);
			tb[(y + 3) * width + x + 1] = (s32)((-36 * v2 - 83 * v3) >> 8);
			tb[(y + 5) * width + x + 1] = (s32)((64 * (v0 - v1)) >> 8);
			tb[(y + 7) * width + x + 1] = (s32)((83 * v2 - 36 * v3) >> 8);
		}
	}
}

void idct(cpStruct *cp)
{
	// Transform the luma component
	idct(cp->pLuma, cp->width, cp->height, cp->bitshift);

	// Transform the chroma components
	if (cp->chromaSub != CHROMA_400)
	{
		idct(cp->pChroma1, cp->width, cp->height, cp->bitshift);
		idct(cp->pChroma2, cp->width, cp->height, cp->bitshift);
	}
}

void idct(s32 *tb)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < 4; y++)
	{
		v0 = 64 * ((s32)tb[y * 4 + 0] + (s32)tb[y * 4 + 2]);
		v1 = 64 * ((s32)tb[y * 4 + 0] - (s32)tb[y * 4 + 2]);
		v2 = -36 * (s32)tb[y * 4 + 1] + 83 * (s32)tb[y * 4 + 3];
		v3 = -83 * (s32)tb[y * 4 + 1] - 36 * (s32)tb[y * 4 + 3];
		tb[y * 4 + 0] = (s32)((v0 - v3) >> 7);
		tb[y * 4 + 1] = (s32)((v1 - v2) >> 7);
		tb[y * 4 + 2] = (s32)((v1 + v2) >> 7);
		tb[y * 4 + 3] = (s32)((v0 + v3) >> 7);
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 x = 0; x < 4; x++)
	{
		v0 = 64 * ((s32)tb[0 + x] + (s32)tb[8 + x]);
		v1 = 64 * ((s32)tb[0 + x] - (s32)tb[8 + x]);
		v2 = -36 * (s32)tb[4 + x] + 83 * (s32)tb[12 + x];
		v3 = -83 * (s32)tb[4 + x] - 36 * (s32)tb[12 + x];
		tb[0 + x] = (s32)((v0 - v3) >> 4);
		tb[4 + x] = (s32)((v1 - v2) >> 4);
		tb[8 + x] = (s32)((v1 + v2) >> 4);
		tb[12 + x] = (s32)((v0 + v3) >> 4);
	}
}

void idct(s32 *tb, u8 width, u8 height, u8 bitshift)
{
	s32 v0, v1, v2, v3;

	// Perform the horizontal pass of the 4-point even-odd DCT
	for (u8 y = 0; y < height; y += 2)
	{
		for (u8 x = 0; x < width; x += 8)
		{
			v0 = 64 * ((s32)tb[y * width + x + 1] + (s32)tb[y * width + x + 5]);
			v1 = 64 * ((s32)tb[y * width + x + 1] - (s32)tb[y * width + x + 5]);
			v2 = -36 * (s32)tb[y * width + x + 3] + 83 * (s32)tb[y * width + x + 7];
			v3 = -83 * (s32)tb[y * width + x + 3] - 36 * (s32)tb[y * width + x + 7];
			tb[y * width + x + 1] = (s32)((v0 - v3) >> 7);
			tb[y * width + x + 3] = (s32)((v1 - v2) >> 7);
			tb[y * width + x + 5] = (s32)((v1 + v2) >> 7);
			tb[y * width + x + 7] = (s32)((v0 + v3) >> 7);

			v0 = 64 * ((s32)tb[(y + 1) * width + x + 0] + (s32)tb[(y + 1) * width + x + 4]);
			v1 = 64 * ((s32)tb[(y + 1) * width + x + 0] - (s32)tb[(y + 1) * width + x + 4]);
			v2 = -36 * (s32)tb[(y + 1) * width + x + 2] + 83 * (s32)tb[(y + 1) * width + x + 6];
			v3 = -83 * (s32)tb[(y + 1) * width + x + 2] - 36 * (s32)tb[(y + 1) * width + x + 6];
			tb[(y + 1) * width + x + 0] = (s32)((v0 - v3) >> 7);
			tb[(y + 1) * width + x + 2] = (s32)((v1 - v2) >> 7);
			tb[(y + 1) * width + x + 4] = (s32)((v1 + v2) >> 7);
			tb[(y + 1) * width + x + 6] = (s32)((v0 + v3) >> 7);

			v0 = 64 * ((s32)tb[(y + 1) * width + x + 1] + (s32)tb[(y + 1) * width + x + 5]);
			v1 = 64 * ((s32)tb[(y + 1) * width + x + 1] - (s32)tb[(y + 1) * width + x + 5]);
			v2 = -36 * (s32)tb[(y + 1) * width + x + 3] + 83 * (s32)tb[(y + 1) * width + x + 7];
			v3 = -83 * (s32)tb[(y + 1) * width + x + 3] - 36 * (s32)tb[(y + 1) * width + x + 7];
			tb[(y + 1) * width + x + 1] = (s32)((v0 - v3) >> 7);
			tb[(y + 1) * width + x + 3] = (s32)((v1 - v2) >> 7);
			tb[(y + 1) * width + x + 5] = (s32)((v1 + v2) >> 7);
			tb[(y + 1) * width + x + 7] = (s32)((v0 + v3) >> 7);
		}
	}

	// Perform the vertical pass of the 4-point even-odd DCT
	for (u8 y = 0; y < height; y += 8)
	{
		for (u8 x = 0; x < width; x += 2)
		{
			v0 = 64 * ((s32)tb[(y + 0) * width + x + 1] + (s32)tb[(y + 4) * width + x + 1]);
			v1 = 64 * ((s32)tb[(y + 0) * width + x + 1] - (s32)tb[(y + 4) * width + x + 1]);
			v2 = -36 * (s32)tb[(y + 2) * width + x + 1] + 83 * (s32)tb[(y + 6) * width + x + 1];
			v3 = -83 * (s32)tb[(y + 2) * width + x + 1] - 36 * (s32)tb[(y + 6) * width + x + 1];
			tb[(y + 0) * width + x + 1] = (s32)((v0 - v3) >> (12 - bitshift));
			tb[(y + 2) * width + x + 1] = (s32)((v1 - v2) >> (12 - bitshift));
			tb[(y + 4) * width + x + 1] = (s32)((v1 + v2) >> (12 - bitshift));
			tb[(y + 6) * width + x + 1] = (s32)((v0 + v3) >> (12 - bitshift));

			v0 = 64 * ((s32)tb[(y + 1) * width + x] + (s32)tb[(y + 5) * width + x]);
			v1 = 64 * ((s32)tb[(y + 1) * width + x] - (s32)tb[(y + 5) * width + x]);
			v2 = -36 * (s32)tb[(y + 3) * width + x] + 83 * (s32)tb[(y + 7) * width + x];
			v3 = -83 * (s32)tb[(y + 3) * width + x] - 36 * (s32)tb[(y + 7) * width + x];
			tb[(y + 1) * width + x] = (s32)((v0 - v3) >> (12 - bitshift));
			tb[(y + 3) * width + x] = (s32)((v1 - v2) >> (12 - bitshift));
			tb[(y + 5) * width + x] = (s32)((v1 + v2) >> (12 - bitshift));
			tb[(y + 7) * width + x] = (s32)((v0 + v3) >> (12 - bitshift));

			v0 = 64 * ((s32)tb[(y + 1) * width + x + 1] + (s32)tb[(y + 5) * width + x + 1]);
			v1 = 64 * ((s32)tb[(y + 1) * width + x + 1] - (s32)tb[(y + 5) * width + x + 1]);
			v2 = -36 * (s32)tb[(y + 3) * width + x + 1] + 83 * (s32)tb[(y + 7) * width + x + 1];
			v3 = -83 * (s32)tb[(y + 3) * width + x + 1] - 36 * (s32)tb[(y + 7) * width + x + 1];
			tb[(y + 1) * width + x + 1] = (s32)((v0 - v3) >> (12 - bitshift));
			tb[(y + 3) * width + x + 1] = (s32)((v1 - v2) >> (12 - bitshift));
			tb[(y + 5) * width + x + 1] = (s32)((v1 + v2) >> (12 - bitshift));
			tb[(y + 7) * width + x + 1] = (s32)((v0 + v3) >> (12 - bitshift));
		}
	}
}
