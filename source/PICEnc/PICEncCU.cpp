#include "PICEncCU.h"
#include "PICComBasis.h"
#include "PICComTransform.h"
#include "PICComQuant.h"
//#include "PICComBAC.h"
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
	u8 qp = param.qp;
	u8 shift = 0;
	for (u8 sWidth = CU_SIZE; sWidth > 4; sWidth >>= 1)
	{
		// Transform the current scale
		dct(&cu, sWidth, shift);
		shift += 2;

		// Find the locally-optimal basis functions
		basisSearch(&cu, sWidth, qp);

		// Fill the averages
		fillAverage(&cu, bufLuma, bufChroma1, bufChroma2, sWidth);
	}

	// Transform the smallest scale
	//dct(cu.pLuma, qp);
	//idct(cu.pLuma, qp);

	// Encode and reconstruct from the smallest scale
	for (u8 sWidth = 8; sWidth <= 64; sWidth <<= 1)
	{
		// Encode the coefficients at the current scale
		// Encode the basis selectors at the current scale
		// Reverse the transformation
		idct(&cu, sWidth, qp, shift);
		shift -= 2;

		// Reconstruct the next level
		basisInverse(&cu, sWidth);
	}
	//// Generate the coefficients in the CP
	//generateCP(img, width, height, cpLargest);

	//// Encode the coarsest scale
	//dct(cpSmallest->pLuma);
	//quantConst(cpSmallest->pLuma, param.qp);
	//coreEnc(cpSmallest->pLuma);

	//// Reconstruct the coarsest scale
	//dequantConst(cpSmallest->pLuma, param.qp);
	//idct(cpSmallest->pLuma);

	//// Encode each Coding Pyramid Slice (CPS) starting from the coarsest scale
	//cp = cpSmallest;
	//u8 qp = param.qp;
	//bool isHorizontalPass = true;
	//do
	//{
	//	if (isHorizontalPass)
	//	{
	//		// Search for the best predictors
	//		predSearchHorz(cp);

	//		// Transform and quantize the residuals
	//		dct(cp->larger);
	//		quantConst(cp->larger, qp);

	//		// Encode the residuals and prediction modes
	//		//resEnc(cp->larger);
	//		//modeEnc(cp->larger);

	//		// Reconstruct the larger slice
	//		dequantConst(cp->larger, qp);
	//		idct(cp->larger);
	//		predHorz(cp);
	//	}
	//	else
	//	{
	//		// Search for the best predictors
	//		predSearchVert(cp);

	//		// Transform and quantize the residuals
	//		dct(cp->larger);
	//		quantConst(cp->larger, qp);

	//		// Encode the residuals and prediction modes
	//		//resEnc(cp->larger);
	//		//modeEnc(cp->larger);

	//		// Reconstruct the larger slice
	//		dequantConst(cp->larger, qp);
	//		idct(cp->larger);
	//		predVert(cp);
	//	}
	//	isHorizontalPass = !isHorizontalPass;

	//	// Apply an optional filter to the reconstructed slice
	//	//sliceFilt(cp->larger, 4);

	//	// Proceed to the next slice
	//	cp = cp->larger;
	//	qp += 3;
	//} while (cp->larger != NULL);

	// Copy to the image
	for (u8 y = 0; y < CU_SIZE; y++)
	{
		for (u8 x = 0; x < CU_SIZE; x++)
		{
			img[width * y + x] = cu.pLuma[CU_SIZE * y + x];
		}
	}
}

void fillAverage(const cuStruct *cu, pel *pLuma, pel *pChroma1, pel *pChroma2, const u8 sWidth)
{
	// Assign the height of all CU scales to be equal to the width
	const u8 sHeight = sWidth;

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