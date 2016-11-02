#include "PICEncCU.h"
#include "PICComPred.h"
#include "PICComTransform.h"
#include "PICComQuant.h"
#include "PICComBAC.h"
#include <cstring>

void PICEncCU(s16 *img, u16 width, u16 height, paramStruct param)
{
	// Create the Coding Pyramid (CP) structure
	cpStruct *cpSmallest, *cpLargest, *cp;

	// Begin with a 4x4 Coding Pyramid Slice (CPS) and initialize successive levels
	cp = new cpStruct;
	cpLargest = cp;
	cpSmallest = cp;
	cpSmallest->smaller = NULL;
	for (u8 i = 4, bs = 8; i <= CU_SIZE; i *= 2)
	{
		// Allocate space for the CPS coefficients
		u16 cpsSize = i * i;
		cp->pLuma = new s32[cpsSize];
		cp->pChroma1 = new s32[cpsSize];
		cp->pChroma2 = new s32[cpsSize];
		cp->predMode = new u8[cpsSize / 16];

		// Write the CPS parameters
		cp->width = i;
		cp->height = i;
		cp->chromaSub = param.chromaSub;

		// Assign the bitshift for the DCT
		cp->bitshift = bs;
		bs -= 2;

		// Allocate the next CPS
		cpLargest = cp;
		cp = new cpStruct;
		cp->smaller = cpLargest;
		cpLargest->larger = cp;
	}
	cpLargest->larger = NULL;
	delete(cp);

	// Generate the coefficients in the CP
	generateCP(img, width, height, cpLargest);

	// Encode the coarsest scale
	dct(cpSmallest->pLuma);
	quantConst(cpSmallest->pLuma, param.qp);
	coreEnc(cpSmallest->pLuma);

	// Reconstruct the coarsest scale
	dequantConst(cpSmallest->pLuma, param.qp);
	idct(cpSmallest->pLuma);

	// Encode each Coding Pyramid Slice (CPS) starting from the coarsest scale
	cp = cpSmallest;
	u8 qp = param.qp;
	while (cp->width < CU_SIZE)
	{
		// Search for the best predictors for the luma component
		predSearch(cp);

		// Transform and quantize the residuals
		dct(cp->larger);
		quantConst(cp->larger, qp);

		// Encode the residuals and prediction modes
		resEnc(cp->larger);
		//modeEnc(cp->larger);

		// Reconstruct the larger slice
		dequantConst(cp->larger, qp);
		idct(cp->larger);
		pred(cp);

		// Apply an optional filter to the reconstructed slice
		sliceFilt(cp->larger, 8);

		// Proceed to the next slice
		cp = cp->larger;
		qp += 6;
	}

	// Copy to the image
	for (u8 y = 0; y < CU_SIZE; y++)
	{
		for (u8 x = 0; x < CU_SIZE; x++)
		{
			img[width * y + x] = cpLargest->pLuma[cpLargest->width * y + x];
		}
	}
}

void generateCP(s16 *img, u16 width, u16 height, cpStruct *cpLargest)
{
	// Copy all luma pixels to the finest scale
	for (u8 y = 0; y < CU_SIZE; y++)
	{
		for (u8 x = 0; x < CU_SIZE; x++)
		{
			cpLargest->pLuma[cpLargest->width * y + x] = img[width * y + x];
		}
	}

	// Generate each slice of the CP
	cpStruct *cp = cpLargest;
	while (cp->width > 4)
	{
		for (u8 y = 0; y < cp->smaller->height; y++)
		{
			for (u8 x = 0; x < cp->smaller->width; x++)
			{
				cp->smaller->pLuma[cp->smaller->width * y + x] = cp->pLuma[cp->width * (2 * y) + (2 * x)] + cp->pLuma[cp->width * (2 * y) + (2 * x + 1)] + cp->pLuma[cp->width * (2 * y + 1) + (2 * x)] + cp->pLuma[cp->width * (2 * y + 1) + (2 * x + 1)];
			}
		}
		cp = cp->smaller;
	}
}

void sliceFilt(cpStruct *cp, u8 preserve)
{
	// Filter the CPS
	for (u8 y = 1; y < cp->height - 1; y++)
	{
		for (u8 x = 1; x < cp->width - 1; x++)
		{
			cp->pLuma[cp->width * y + x] *= preserve;
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * (y - 1) + (x - 1)];
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * (y - 1) + x];
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * (y - 1) + (x + 1)];
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * y + (x - 1)];
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * y + (x + 1)];
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * (y + 1) + (x - 1)];
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * (y + 1) + x];
			cp->pLuma[cp->width * y + x] += cp->pLuma[cp->width * (y + 1) + (x + 1)];
			cp->pLuma[cp->width * y + x] /= (preserve + 8);
		}
	}
}

