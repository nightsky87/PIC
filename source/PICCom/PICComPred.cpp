#include "PICComPred.h"
#include "PICComPredLookup.h"
#include <cmath>

void predHorz(cpStruct *cp)
{
	// Create a padded copy of the current CP pixels and transpose
	u8 pWidth = cp->height + 2;
	u8 pHeight = 2 * cp->width + 1;
	s32 *cpPix = new s32[pWidth * pHeight];
	memset(cpPix, 0, pWidth * pHeight * sizeof(s32));
	for (u8 y = 0; y < cp->height; y++)
	{
		for (u8 x = 0; x < cp->width; x++)
		{
			cpPix[pWidth * 2 * y + x + 1] = cp->pLuma[cp->width * x + y];
		}
		cpPix[pWidth * 2 * y] = cpPix[pWidth * 2 * y + 1];
		cpPix[pWidth * 2 * y + pWidth - 1] = cpPix[pWidth * 2 * y + pWidth - 2];
	}
	memcpy(&cpPix[pWidth * (pHeight - 1)], &cpPix[pWidth * (pHeight - 3)], pWidth * sizeof(s32));

	// Process each prediction block
	u8 ind = 0;
	for (u8 y = 0; y < pHeight - 1; y += 8)
	{
		for (u8 x = 1; x < pWidth - 1; x += 4)
		{
			predBlock(&cpPix[pWidth * y + x], pWidth, cp->predMode[ind]);
			ind++;
		}
	}

	// Calculate the final prediction residual
	for (u8 y = 0; y < cp->larger->height; y++)
	{
		for (u8 x = 0; x < cp->larger->width; x += 2)
		{
			cp->larger->pLuma[cp->larger->width * y + x + 1] += cpPix[pWidth * (x + 1) + y + 1];
			cp->larger->pLuma[cp->larger->width * y + x] -= cp->larger->pLuma[cp->larger->width * y + x + 1];
		}
	}

	// Garbage collection
	delete(cpPix);
}

void predVert(cpStruct *cp)
{
	// Create a padded copy of the current CP pixels
	u8 pWidth = cp->width + 2;
	u8 pHeight = 2 * cp->height + 1;
	s32 *cpPix = new s32[pWidth * pHeight];
	memset(cpPix, 0, pWidth * pHeight * sizeof(s32));
	for (u8 y = 0; y < cp->height; y++)
	{
		for (u8 x = 0; x < cp->width; x++)
		{
			cpPix[pWidth * 2 * y + x + 1] = cp->pLuma[cp->width * y + x];
		}
		cpPix[pWidth * 2 * y] = cpPix[pWidth * 2 * y + 1];
		cpPix[pWidth * 2 * y + pWidth - 1] = cpPix[pWidth * 2 * y + pWidth - 2];
	}
	memcpy(&cpPix[pWidth * (pHeight - 1)], &cpPix[pWidth * (pHeight - 3)], pWidth * sizeof(s32));

	// Process each prediction block
	u8 ind = 0;
	for (u8 y = 0; y < pHeight - 1; y += 8)
	{
		for (u8 x = 1; x < pWidth - 1; x += 4)
		{
			predBlock(&cpPix[pWidth * y + x], pWidth, cp->predMode[ind]);
			ind++;
		}
	}

	// Calculate the final prediction residual
	for (u8 y = 0; y < cp->larger->height; y += 2)
	{
		for (u8 x = 0; x < cp->larger->width; x++)
		{
			cp->larger->pLuma[cp->larger->width * (y + 1) + x] += cpPix[pWidth * (y + 1) + x + 1];
			cp->larger->pLuma[cp->larger->width * y + x] -= cp->larger->pLuma[cp->larger->width * (y + 1) + x];
		}
	}
}

void predBlock(s32 *pb, u8 stride, u8 mode)
{
	if (mode == 0)
		predNN(pb, stride);
	else
		predInt(pb, stride, mode - 1);
}

void predInt(s32 *pb, u8 stride, u8 intNum)
{
	const s16 *predTable = &pred_weights[6 * intNum];

	for (u8 y = 0; y < 8; y += 2)
	{
		for (u8 x = 0; x < 4; x++)
		{
			pb[stride * (y + 1) + x] = predTable[0] * pb[stride * y + x - 1];
			pb[stride * (y + 1) + x] += predTable[1] * pb[stride * y + x];
			pb[stride * (y + 1) + x] += predTable[2] * pb[stride * y + x + 1];
			pb[stride * (y + 1) + x] += predTable[3] * pb[stride * (y + 2) + x - 1];
			pb[stride * (y + 1) + x] += predTable[4] * pb[stride * (y + 2) + x];
			pb[stride * (y + 1) + x] += predTable[5] * pb[stride * (y + 2) + x + 1];
			pb[stride * (y + 1) + x] = (pb[stride * (y + 1) + x] + 32) / 64;
		}
	}
}

void predNN(s32 *pb, u8 stride)
{
	// Copy the averages
	for (u8 y = 0; y < 8; y += 2)
	{
		for (u8 x = 0; x < 4; x++)
		{
			pb[stride * (y + 1) + x] = (pb[stride * y + x] + 1) / 2;
		}
	}
}

void predSearchHorz(cpStruct *cp)
{
	// Create a padded copy of the current CP pixels and transpose
	u8 pWidth = cp->height + 2;
	u8 pHeight = 2 * cp->width + 1;
	s32 *cpPix = new s32[pWidth * pHeight];
	memset(cpPix, 0, pWidth * pHeight * sizeof(s32));
	for (u8 y = 0; y < cp->height; y++)
	{
		for (u8 x = 0; x < cp->width; x++)
		{
			cpPix[pWidth * 2 * y + x + 1] = cp->pLuma[cp->width * x + y];
		}
		cpPix[pWidth * 2 * y] = cpPix[pWidth * 2 * y + 1];
		cpPix[pWidth * 2 * y + pWidth - 1] = cpPix[pWidth * 2 * y + pWidth - 2];
	}
	memcpy(&cpPix[pWidth * (pHeight - 1)], &cpPix[pWidth * (pHeight - 3)], pWidth * sizeof(s32));

	// Create a transposed copy of the target pixels
	s32 *cpTarget = new s32[cp->larger->width * cp->larger->height];
	for (u8 y = 0; y < cp->larger->height; y++)
	{
		for (u8 x = 0; x < cp->larger->width; x++)
		{
			cpTarget[cp->height * x + y] = cp->larger->pLuma[cp->larger->width * y + x];
		}
	}

	// Process each prediction block
	u8 ind = 0;
	for (u8 y = 0; y < pHeight - 1; y += 8)
	{
		for (u8 x = 1; x < pWidth - 1; x += 4)
		{
			cp->predMode[ind] = predSearchBlock(&cpPix[pWidth * y + x], pWidth, &cpTarget[cp->height * y + x - 1], cp->height);
			predBlock(&cpPix[pWidth * y + x], pWidth, cp->predMode[ind]);
			ind++;
		}
	}

	// Calculate the final prediction residual
	for (u8 y = 0; y < cp->larger->height; y++)
	{
		for (u8 x = 0; x < cp->larger->width; x += 2)
		{
			cp->larger->pLuma[cp->larger->width * y + x] = cpPix[pWidth * x + y + 1];
			cp->larger->pLuma[cp->larger->width * y + x + 1] -= cpPix[pWidth * (x + 1) + y + 1];
		}
	}
}

void predSearchVert(cpStruct *cp)
{
	// Create a padded copy of the current CP pixels
	u8 pWidth = cp->width + 2;
	u8 pHeight = 2 * cp->height + 1;
	s32 *cpPix = new s32[pWidth * pHeight];
	memset(cpPix, 0, pWidth * pHeight * sizeof(s32));
	for (u8 y = 0; y < cp->height; y++)
	{
		for (u8 x = 0; x < cp->width; x++)
		{
			cpPix[pWidth * 2 * y + x + 1] = cp->pLuma[cp->width * y + x];
		}
		cpPix[pWidth * 2 * y] = cpPix[pWidth * 2 * y + 1];
		cpPix[pWidth * 2 * y + pWidth - 1] = cpPix[pWidth * 2 * y + pWidth - 2];
	}
	memcpy(&cpPix[pWidth * (pHeight - 1)], &cpPix[pWidth * (pHeight - 3)], pWidth * sizeof(s32));

	// Process each prediction block
	u8 ind = 0;
	for (u8 y = 0; y < pHeight - 1; y += 8)
	{
		for (u8 x = 1; x < pWidth - 1; x += 4)
		{
			cp->predMode[ind] = predSearchBlock(&cpPix[pWidth * y + x], pWidth, &cp->larger->pLuma[cp->larger->width * y + x - 1], cp->larger->width);
			predBlock(&cpPix[pWidth * y + x], pWidth, cp->predMode[ind]);
			ind++;
		}
	}

	// Calculate the final prediction residual
	for (u8 y = 0; y < cp->larger->height; y += 2)
	{
		for (u8 x = 0; x < cp->larger->width; x++)
		{
			cp->larger->pLuma[cp->larger->width * y + x] = cpPix[pWidth * y + x + 1];
			cp->larger->pLuma[cp->larger->width * (y + 1) + x] -= cpPix[pWidth * (y + 1) + x + 1];
		}
	}
}

u8 predSearchBlock(s32 *pb, u8 stride, s32 *pbTrue, u8 strideTrue)
{
	// Initialize the best mode to nearest-neighbor
	u8 bestMode = 0;

	// Test with nearest-neighbor interpolation
	predNN(pb, stride);
	u32 bestMetric = hadamardMetric(pb, stride, pbTrue, strideTrue);

	// Test with the other interpolants
	for (u8 mode = 1; mode <= 16; mode++)
	{
		predInt(pb, stride, mode - 1);
		u32 metric = hadamardMetric(pb, stride, pbTrue, strideTrue);

		if (metric < bestMetric)
		{
			bestMetric = metric;
			bestMode = mode;
		}
	}

	return bestMode;
}

u32 hadamardMetric(s32 *pb, u8 stride, s32 *pbTrue, u8 strideTrue)
{
	// Calculate the residual
	for (u8 y = 1; y < 8; y += 2)
	{
		for (u8 x = 0; x < 4; x++)
		{
			pb[stride * y + x] -= pbTrue[strideTrue * y + x];
		}
	}
	// Metric calculation
	u32 metric2 = 0;
	for (u8 y = 1; y < 8; y += 2)
	{
		for (u8 x = 0; x < 4; x++)
		{
			s32 tmp = pb[stride * y + x] < 0 ? -pb[stride * y + x] : pb[stride * y + x];
			metric2 += (u32)tmp;
		}
	}
	return metric2;

	// Perform the row-wise Hadamard transform
	for (u8 y = 1; y < 8; y += 2)
	{
		// First pass
		for (u8 x = 0; x < 2; x++)
		{
			s32 tmp = pb[stride * y + x];
			pb[stride * y + x] += pb[stride * y + x + 2];
			pb[stride * y + x + 2] = tmp - pb[stride * y + x + 2];
		}

		// Second pass
		for (u8 x = 0; x < 4; x += 2)
		{
			s32 tmp = pb[stride * y + x];
			pb[stride * y + x] += pb[stride * y + x + 1];
			pb[stride * y + x + 1] = tmp - pb[stride * y + x + 1];
		}
	}

	// Perform the column-wise Hadamard transform
	u32 metric = 0;
	for (u8 x = 0; x < 4; x++)
	{
		// First pass
		for (u8 y = 1; y < 4; y += 2)
		{
			s32 tmp = pb[stride * y + x];
			pb[stride * y + x] += pb[stride * (y + 4) + x];
			pb[stride * (y + 4) + x] = tmp - pb[stride * (y + 4) + x];
		}

		// Second pass
		for (u8 y = 1; y < 8; y += 4)
		{
			s32 tmp = pb[stride * y + x];
			pb[stride * y + x] += pb[stride * (y + 2) + x];
			pb[stride * (y + 2) + x] = tmp - pb[stride * (y + 2) + x];
		}

		// Metric calculation
		for (u8 y = 1; y < 8; y += 2)
		{
			s32 tmp = pb[stride * y + x] < 0 ? -pb[stride * y + x] : pb[stride * y + x];
			metric += (u32)tmp;
		}
	}

	return metric;
}
