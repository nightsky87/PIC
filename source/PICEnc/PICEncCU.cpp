#include "PICEncCU.h"
#include "PICComBasis.h"
#include "PICComTransform.h"
#include "PICComQuant.h"
#include "PICComBAC.h"
#include <cstring>

pel PICEncCU(s16 *img, u16 width, u16 height, paramStruct param, pel prevDC)
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

	u8 quadInc = 6;

	// Begin with finest scale and progressively transform each scale
	cu.basisIndex = 0;
	u8 qp = param.qp + 4 * quadInc;
	u8 shift = 0;
	for (u8 sWidth = CU_SIZE; sWidth > 8; sWidth >>= 1)
	{
		// Find the locally-optimal interpolants and basis functions
		basisSearch(&cu, sWidth, sWidth);

		// Transform the current scale
		dct(&cu, sWidth, sWidth, shift, qp);

		shift += 2;
		qp -= quadInc;
	}

	// Transform the core block
	dct(cu.pLuma, 8, 8, shift, qp);

	// Express the DC value relative to the previous CU
	cu.pLuma[0] -= prevDC;

	// Encode the core block
	coreEnc(cu.pLuma);

	// Reconstruct the core block
	cu.pLuma[0] += prevDC;
	prevDC = cu.pLuma[0];
	idct(cu.pLuma, 8, 8, shift, qp);

	// Encode and reconstruct from the smallest scale
	u8 quadID = 0;
	for (u8 sWidth = 16; sWidth <= 64; sWidth <<= 1)
	{
		qp += quadInc;
		shift -= 2;

		// Encode the coefficients at the current scale
		quadEnc(&cu, sWidth, sWidth, quadID);
		quadID += 3;

		// Encode the basis selectors at the current scale
		
		// Reverse the transformation
		idct(&cu, sWidth, sWidth, shift, qp);

		// Reconstruct the next level
		basisInverse(&cu, sWidth, sWidth);
	}

	// Copy to the image
	for (u8 y = 0; y < CU_SIZE; y++)
	{
		for (u8 x = 0; x < CU_SIZE; x++)
		{
			img[width * y + x] = cu.pLuma[CU_SIZE * y + x];
		}
	}

	return prevDC;
}
