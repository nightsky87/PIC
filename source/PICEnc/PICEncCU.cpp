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
	cp->larger = new cpStruct;
	cp->larger->smaller = cp;
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

		cp->larger->pLuma = new s32[2 * cpsSize];
		cp->larger->pChroma1 = new s32[2 * cpsSize];
		cp->larger->pChroma2 = new s32[2 * cpsSize];
		cp->larger->predMode = new u8[cpsSize / 8];

		// Write the CPS parameters
		cp->width = i;
		cp->height = i;
		cp->chromaSub = param.chromaSub;

		cp->larger->width = 2 * i;
		cp->larger->height = i;
		cp->larger->chromaSub = param.chromaSub;

		// Assign the bitshift for the DCT
		cp->bitshift = bs;
		cp->larger->bitshift = bs - 1;
		bs -= 2;

		// Allocate the next two CPS
		cp->larger->larger = new cpStruct;
		cp->larger->larger->larger = new cpStruct;

		// Link in the reverse direction
		cp->larger->larger->larger->smaller = cp->larger->larger;
		cp->larger->larger->smaller = cp->larger;
		cp = cp->larger->larger;
	}
	cpLargest = cp->smaller->smaller;
	delete(cp->smaller);
	delete(cp);
	cpLargest->larger = NULL;

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
		sliceFilt(cp->larger, 32);

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
		u8 width = cp->width;
		u8 height = cp->height / 2;

		for (u8 y = 0; y < height; y++)
		{
			for (u8 x = 0; x < width; x++)
			{
				cp->smaller->pLuma[width * y + x] = cp->pLuma[width * (2 * y) + x] + cp->pLuma[width * (2 * y + 1) + x];
			}
		}
		cp = cp->smaller;

		for (u8 y = 0; y < height; y++)
		{
			for (u8 x = 0; x < width / 2; x++)
			{
				cp->smaller->pLuma[width / 2 * y + x] = cp->pLuma[width * y + (2 * x)] + cp->pLuma[width * y + (2 * x + 1)];
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

