#include "PICEncCU.h"
#include "PICComBasis.h"
#include "PICComTransform.h"
#include "PICComQuant.h"
#include "PICComBAC.h"
#include <cstring>

void PICEncCU(s16 *img, u16 width, u16 height, paramStruct param)
{
	// Allocate static memory for the coefficient data
	static pel *cuLuma = new pel[CU_SIZE * CU_SIZE];
	static pel *cuChroma1 = new pel[CU_SIZE * CU_SIZE];
	static pel *cuChroma2 = new pel[CU_SIZE * CU_SIZE];

	// Allocate space for the basis selectors
	static u8 *basisSelect = new u8[NUM_BASIS_SELECTORS];

	// Create a static Coding Unit structure
	static cuStruct cu = { cuLuma, cuChroma1, cuChroma2, basisSelect, 0, CU_SIZE, CU_SIZE, param.chromaSub };

	// Copy the pixels to the CU
	if (param.chromaSub == CHROMA_400)
	{
		for (u8 y = 0; y < CU_SIZE; y++)
		{
			for (u8 x = 0; x < CU_SIZE; x++)
			{
				cu.pLuma[CU_SIZE * y + x] = img[width * y + x];
			}
		}
	}
	else
	{
		const u32 chOffset = width * height;
		for (u8 y = 0; y < CU_SIZE; y++)
		{
			for (u8 x = 0; x < CU_SIZE; x++)
			{
				cu.pLuma[CU_SIZE * y + x] = img[width * y + x];
				cu.pChroma1[CU_SIZE * y + x] = img[width * y + x + chOffset];
				cu.pChroma2[CU_SIZE * y + x] = img[width * y + x + 2 * chOffset];
			}
		}
	}

	// Allocate static memory for a buffer of each scale
	static pel *bufLuma = new pel[CU_SIZE * CU_SIZE];
	static pel *bufChroma1 = new pel[CU_SIZE * CU_SIZE];
	static pel *bufChroma2 = new pel[CU_SIZE * CU_SIZE];

	// Copy all the data
	if (param.chromaSub == CHROMA_400)
	{
		memcpy(bufLuma, cu.pLuma, CU_SIZE * CU_SIZE * sizeof(pel));
	}
	else
	{
		memcpy(bufLuma, cu.pLuma, CU_SIZE * CU_SIZE * sizeof(pel));
		memcpy(bufChroma1, cu.pChroma1, CU_SIZE * CU_SIZE * sizeof(pel));
		memcpy(bufChroma2, cu.pChroma2, CU_SIZE * CU_SIZE * sizeof(pel));
	}

	// Begin with finest scale and progressively transform each scale
	cu.basisIndex = 0;
	u8 qp = param.qp + 24;
	u8 shift = 0;
	for (u8 sWidth = CU_SIZE; sWidth > 4; sWidth >>= 1)
	{
		// Transform the current scale
		dct(&cu, sWidth, sWidth, shift, qp);

		// Find the locally-optimal basis functions
		basisSearch(&cu, sWidth, sWidth, qp);
		reduce(&cu, bufLuma, bufChroma1, bufChroma2, sWidth, sWidth);

		shift += 2;
		qp -= 6;
	}

	// Transform and encode the core block
	dct(cu.pLuma, shift, qp);
	coreEnc(cu.pLuma);
	idct(cu.pLuma, shift, qp);

	// Encode and reconstruct from the smallest scale
	for (u8 sWidth = 8; sWidth <= 64; sWidth <<= 1)
	{
		// Encode the coefficients at the current scale
		// Encode the basis selectors at the current scale
		// Reverse the transformation
		qp += 6;
		shift -= 2;
		idct(&cu, sWidth, sWidth, shift, qp);

		// Reconstruct the next level
		basisInverse(&cu, sWidth, sWidth);

		// Apply an optional filter to the reconstructed scale
		//filter(&cu, sWidth, sWidth, 8);
	}

	// Copy to the image
	for (u8 y = 0; y < CU_SIZE; y++)
	{
		for (u8 x = 0; x < CU_SIZE; x++)
		{
			img[width * y + x] = cu.pLuma[CU_SIZE * y + x];
		}
	}
}

void reduce(const cuStruct *cu, pel *pLuma, pel *pChroma1, pel *pChroma2, const u8 sWidth, const u8 sHeight)
{
	// Process each 2 x 2 region
	if (cu->chromaSub == CHROMA_400)
	{
		for (u8 y = 0; y < sHeight; y += 2)
		{
			for (u8 x = 0; x < sWidth; x += 2)
			{
				cu->pLuma[CU_SIZE * y / 2 + x / 2] = pLuma[CU_SIZE * y + x] + pLuma[CU_SIZE * y + x + 1] + pLuma[CU_SIZE * (y + 1) + x] + pLuma[CU_SIZE * (y + 1) + x + 1];
			}
		}

		for (u8 y = 0; y < sHeight / 2; y++)
		{
			memcpy(&pLuma[CU_SIZE * y], &cu->pLuma[CU_SIZE * y], sWidth / 2 * sizeof(pel));
		}
	}
	else
	{
		for (u8 y = 0; y < sHeight; y += 2)
		{
			for (u8 x = 0; x < sWidth; x += 2)
			{
				cu->pLuma[CU_SIZE * y / 2 + x / 2] = pLuma[CU_SIZE * y + x] + pLuma[CU_SIZE * y + x + 1] + pLuma[CU_SIZE * (y + 1) + x] + pLuma[CU_SIZE * (y + 1) + x + 1];
				cu->pChroma1[CU_SIZE * y / 2 + x / 2] = pChroma1[CU_SIZE * y + x] + pChroma1[CU_SIZE * y + x + 1] + pChroma1[CU_SIZE * (y + 1) + x] + pChroma1[CU_SIZE * (y + 1) + x + 1];
				cu->pChroma2[CU_SIZE * y / 2 + x / 2] = pChroma2[CU_SIZE * y + x] + pChroma2[CU_SIZE * y + x + 1] + pChroma2[CU_SIZE * (y + 1) + x] + pChroma2[CU_SIZE * (y + 1) + x + 1];
			}
		}

		for (u8 y = 0; y < sHeight / 2; y++)
		{
			memcpy(&pLuma[CU_SIZE * y], &cu->pLuma[CU_SIZE * y], sWidth / 2 * sizeof(pel));
			memcpy(&pChroma1[CU_SIZE * y], &cu->pChroma1[CU_SIZE * y], sWidth / 2 * sizeof(pel));
			memcpy(&pChroma2[CU_SIZE * y], &cu->pChroma2[CU_SIZE * y], sWidth / 2 * sizeof(pel));
		}
	}
}

void filter(cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 preserve)
{
	for (u8 y = 1; y < sHeight - 1; y++)
	{
		for (u8 x = 1; x < sWidth - 1; x++)
		{
			cu->pLuma[CU_SIZE * y + x] *= preserve;
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * (y - 1) + (x - 1)];
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * (y - 1) + x];
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * (y - 1) + (x + 1)];
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * y + (x - 1)];
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * y + (x + 1)];
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * (y + 1) + (x - 1)];
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * (y + 1) + x];
			cu->pLuma[CU_SIZE * y + x] += cu->pLuma[CU_SIZE * (y + 1) + (x + 1)];
			cu->pLuma[CU_SIZE * y + x] /= (preserve + 8);
		}
	}
}