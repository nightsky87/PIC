#include "PICComPred.h"
#include "PICComPredLookup.h"
#include <cmath>

void pred(cpStruct *cp)
{
	// Perform predictions on the luma component
	predLuma(cp);

	// Perform predictions on the chroma components
	if (cp->chromaSub != CHROMA_400)
		predChroma(cp);
}

void predLuma(cpStruct *cp)
{
	// Create a padded copy of the current CP pixels
	u8 pWidth = 2 * cp->width + 3;
	s32 *cpPix = new s32[(2 * cp->height + 3) * pWidth];
	memset(cpPix, 0, (2 * cp->height + 3) * pWidth * sizeof(s32));
	for (u8 y = 0; y < cp->height; y++)
	{
		cpPix[pWidth * 2 * (y + 1)] = cp->pLuma[cp->width * y];
		for (u8 x = 0; x < cp->width; x++)
		{
			cpPix[pWidth * 2 * (y + 1) + 2 * (x + 1)] = cp->pLuma[cp->width * y + x];
		}
		cpPix[pWidth * 2 * (y + 1) + pWidth - 1] = cp->pLuma[cp->width * y + cp->width - 1];
	}
	memcpy(cpPix, &cpPix[2 * pWidth], pWidth * sizeof(s32));
	memcpy(&cpPix[pWidth * (2 * cp->height + 2)], &cpPix[pWidth * (2 * cp->height)], pWidth * sizeof(s32));

	// Process each prediction block
	u8 ind = 0;
	for (u8 y = 2; y < cp->larger->height + 2; y += 8)
	{
		for (u8 x = 2; x < cp->larger->width + 2; x += 8)
		{
			predBlock(&cpPix[pWidth * y + x], pWidth, cp->predMode[ind]);
			ind++;
		}
	}

	// Generate the final prediction with residual correction
	for (u8 y = 0; y < cp->larger->height; y += 2)
	{
		for (u8 x = 0; x < cp->larger->width; x += 2)
		{
			cp->larger->pLuma[cp->larger->width * y + x + 1] += cpPix[pWidth * (y + 2) + x + 3];
			cp->larger->pLuma[cp->larger->width * (y + 1) + x] += cpPix[pWidth * (y + 3) + x + 2];
			cp->larger->pLuma[cp->larger->width * (y + 1) + x + 1] += cpPix[pWidth * (y + 3) + x + 3];

			cp->larger->pLuma[cp->larger->width * y + x] -= (cp->larger->pLuma[cp->larger->width * y + x + 1] + cp->larger->pLuma[cp->larger->width * (y + 1) + x] + cp->larger->pLuma[cp->larger->width * (y + 1) + x + 1]);
		}
	}
}

void predChroma(cpStruct *cp)
{
	// Copy the width and height of the current level
	u8 width = cp->width;
	u8 height = cp->height;

	// Always use nearest-neighbor interpolation for the chroma channel
	u8 lwidth = cp->larger->width;
	cpStruct *cpLarger = cp->larger;
	for (u8 y = 0; y < height; y++)
	{
		for (u8 x = 0; x < width; x++)
		{
			cpLarger->pChroma1[lwidth * (2 * y) + (2 * x)] = cp->pChroma1[width * y + x];
			cpLarger->pChroma2[lwidth * (2 * y) + (2 * x)] = cp->pChroma2[width * y + x];
			cpLarger->pChroma1[lwidth * (2 * y) + (2 * x + 1)] = cp->pChroma1[width * y + x];
			cpLarger->pChroma2[lwidth * (2 * y) + (2 * x + 1)] = cp->pChroma2[width * y + x];
			cpLarger->pChroma1[lwidth * (2 * y + 1) + (2 * x)] = cp->pChroma1[width * y + x];
			cpLarger->pChroma2[lwidth * (2 * y + 1) + (2 * x)] = cp->pChroma2[width * y + x];
			cpLarger->pChroma1[lwidth * (2 * y + 1) + (2 * x + 1)] = cp->pChroma1[width * y + x];
			cpLarger->pChroma2[lwidth * (2 * y + 1) + (2 * x + 1)] = cp->pChroma2[width * y + x];
		}
	}
}

void predChroma(cpStruct *cp, cpStruct *cpDiff)
{
	// Copy the width and height of the current level
	u8 width = cp->width;
	u8 height = cp->height;

	// Copy the known pixels
	memcpy(cp->larger->pChroma1, cpDiff->pChroma1, 4 * width * height * sizeof(s16));
	memcpy(cp->larger->pChroma2, cpDiff->pChroma2, 4 * width * height * sizeof(s16));

	// Always use nearest-neighbor interpolation for the chroma channel
	u8 lwidth = cp->larger->width;
	for (u8 y = 0; y < height; y++)
	{
		for (u8 x = 0; x < width; x++)
		{
			cpDiff->pChroma1[lwidth * (2 * y) + (2 * x)] -= cp->pChroma1[width * y + x];
			cpDiff->pChroma2[lwidth * (2 * y) + (2 * x)] -= cp->pChroma2[width * y + x];
			cpDiff->pChroma1[lwidth * (2 * y) + (2 * x + 1)] -= cp->pChroma1[width * y + x];
			cpDiff->pChroma2[lwidth * (2 * y) + (2 * x + 1)] -= cp->pChroma2[width * y + x];
			cpDiff->pChroma1[lwidth * (2 * y + 1) + (2 * x)] -= cp->pChroma1[width * y + x];
			cpDiff->pChroma2[lwidth * (2 * y + 1) + (2 * x)] -= cp->pChroma2[width * y + x];
			cpDiff->pChroma1[lwidth * (2 * y + 1) + (2 * x + 1)] -= cp->pChroma1[width * y + x];
			cpDiff->pChroma2[lwidth * (2 * y + 1) + (2 * x + 1)] -= cp->pChroma2[width * y + x];
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
	const s16 *predTable1 = &pred_w1[6 * intNum];
	const s16 *predTable2 = &pred_w2[6 * intNum];
	const s16 *predTable3 = &pred_w3[4 * intNum];

	for (u8 y = 0; y < 8; y += 2)
	{
		for (u8 x = 0; x < 8; x += 2)
		{
			pb[stride * y + (x + 1)] = predTable1[0] * pb[stride * (y - 2) + x];
			pb[stride * y + (x + 1)] += predTable1[1] * pb[stride * (y - 2) + x + 2];
			pb[stride * y + (x + 1)] += predTable1[2] * pb[stride * y + x];
			pb[stride * y + (x + 1)] += predTable1[3] * pb[stride * y + x + 2];
			pb[stride * y + (x + 1)] += predTable1[4] * pb[stride * (y + 2) + x];
			pb[stride * y + (x + 1)] += predTable1[5] * pb[stride * (y + 2) + x + 2];
			pb[stride * y + (x + 1)] = (pb[stride * y + (x + 1)] + 32) / 64;

			pb[stride * (y + 1) + x] = predTable2[0] * pb[stride * y + x - 2];
			pb[stride * (y + 1) + x] += predTable2[1] * pb[stride * y + x];
			pb[stride * (y + 1) + x] += predTable2[2] * pb[stride * y + x + 2];
			pb[stride * (y + 1) + x] += predTable2[3] * pb[stride * (y + 2) + x - 2];
			pb[stride * (y + 1) + x] += predTable2[4] * pb[stride * (y + 2) + x];
			pb[stride * (y + 1) + x] += predTable2[5] * pb[stride * (y + 2) + x + 2];
			pb[stride * (y + 1) + x] = (pb[stride * (y + 1) + x] + 32) / 64;

			pb[stride * (y + 1) + (x + 1)] = predTable3[0] * pb[stride * y + x];
			pb[stride * (y + 1) + (x + 1)] += predTable3[1] * pb[stride * y + x + 2];
			pb[stride * (y + 1) + (x + 1)] += predTable3[2] * pb[stride * (y + 2) + x];
			pb[stride * (y + 1) + (x + 1)] += predTable3[3] * pb[stride * (y + 2) + x + 2];
			pb[stride * (y + 1) + (x + 1)] = (pb[stride * (y + 1) + (x + 1)] + 32) / 64;
		}
	}
}

void predNN(s32 *pb, u8 stride)
{
	// Copy the averages
	for (u8 y = 0; y < 8; y += 2)
	{
		for (u8 x = 0; x < 8; x += 2)
		{
			pb[stride * y + (x + 1)] = (pb[stride * y + x] + 2) / 4;
			pb[stride * (y + 1) + x] = (pb[stride * y + x] + 2) / 4;
			pb[stride * (y + 1) + (x + 1)] = (pb[stride * y + x] + 2) / 4;
		}
	}
}

void predSearch(cpStruct *cp)
{
	// Create a padded copy of the current CP pixels
	u8 pWidth = 2 * cp->width + 3;
	s32 *cpPix = new s32[(2 * cp->height + 3) * pWidth];
	memset(cpPix, 0, (2 * cp->height + 3) * pWidth * sizeof(s32));
	for (u8 y = 0; y < cp->height; y++)
	{
		cpPix[pWidth * 2 * (y + 1)] = cp->pLuma[cp->width * y];
		for (u8 x = 0; x < cp->width; x++)
		{
			cpPix[pWidth * 2 * (y + 1) + 2 * (x + 1)] = cp->pLuma[cp->width * y + x];
		}
		cpPix[pWidth * 2 * (y + 1) + pWidth - 1] = cp->pLuma[cp->width * y + cp->width - 1];
	}
	memcpy(cpPix, &cpPix[2 * pWidth], pWidth * sizeof(s32));
	memcpy(&cpPix[pWidth * (2 * cp->height + 2)], &cpPix[pWidth * (2 * cp->height)], pWidth * sizeof(s32));

	// Process each prediction block
	u8 ind = 0;
	for (u8 y = 2; y < cp->larger->height + 2; y += 8)
	{
		for (u8 x = 2; x < cp->larger->width + 2; x += 8)
		{
			cp->predMode[ind] = predSearchBlock(&cpPix[pWidth * y + x], pWidth, &cp->larger->pLuma[cp->larger->width * (y - 2) + x - 2], cp->larger->width);
			predBlock(&cpPix[pWidth * y + x], pWidth, cp->predMode[ind]);
			ind++;
		}
	}

	// Calculate the final prediction residual
	for (u8 y = 0; y < cp->larger->height; y += 2)
	{
		for (u8 x = 0; x < cp->larger->width; x += 2)
		{
			cp->larger->pLuma[cp->larger->width * y + x] = cpPix[pWidth * (y + 2) + x + 2];

			cp->larger->pLuma[cp->larger->width * y + x + 1] -= cpPix[pWidth * (y + 2) + x + 3];
			cp->larger->pLuma[cp->larger->width * (y + 1) + x] -= cpPix[pWidth * (y + 3) + x + 2];
			cp->larger->pLuma[cp->larger->width * (y + 1) + x + 1] -= cpPix[pWidth * (y + 3) + x + 3];
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
	for (u8 y = 0; y < 8; y += 2)
	{
		for (u8 x = 0; x < 8; x += 2)
		{
			pb[stride * y + (x + 1)] -= pbTrue[strideTrue * y + (x + 1)];
			pb[stride * (y + 1) + x] -= pbTrue[strideTrue * (y + 1) + x];
			pb[stride * (y + 1) + (x + 1)] -= pbTrue[strideTrue * (y + 1) + (x + 1)];
		}
	}

	// Perform the row-wise Hadamard transform
	for (u8 y = 0; y < 8; y += 2)
	{
		// First pass
		for (u8 x = 0; x < 4; x += 2)
		{
			s32 tmp1 = pb[stride * y + (x + 1)];
			s32 tmp2 = pb[stride * (y + 1) + x];
			s32 tmp3 = pb[stride * (y + 1) + (x + 1)];

			pb[stride * y + (x + 1)] += pb[stride * y + (x + 1) + 4];
			pb[stride * (y + 1) + x] += pb[stride * (y + 1) + x + 4];
			pb[stride * (y + 1) + (x + 1)] += pb[stride * (y + 1) + (x + 1) + 4];

			pb[stride * y + (x + 1) + 4] = tmp1 - pb[stride * y + (x + 1) + 4];
			pb[stride * (y + 1) + x + 4] = tmp2 - pb[stride * (y + 1) + x + 4];
			pb[stride * (y + 1) + (x + 1) + 4] = tmp3 - pb[stride * (y + 1) + (x + 1) + 4];
		}

		// Second pass
		for (u8 x = 0; x < 8; x += 4)
		{
			s32 tmp1 = pb[stride * y + (x + 1)];
			s32 tmp2 = pb[stride * (y + 1) + x];
			s32 tmp3 = pb[stride * (y + 1) + (x + 1)];

			pb[stride * y + (x + 1)] += pb[stride * y + (x + 1) + 2];
			pb[stride * (y + 1) + x] += pb[stride * (y + 1) + x + 2];
			pb[stride * (y + 1) + (x + 1)] += pb[stride * (y + 1) + (x + 1) + 2];

			pb[stride * y + (x + 1) + 2] = tmp1 - pb[stride * y + (x + 1) + 2];
			pb[stride * (y + 1) + x + 2] = tmp2 - pb[stride * (y + 1) + x + 2];
			pb[stride * (y + 1) + (x + 1) + 2] = tmp3 - pb[stride * (y + 1) + (x + 1) + 2];
		}
	}

	// Perform the row-wise Hadamard transform
	u32 metric = 0;
	for (u8 x = 0; x < 8; x += 2)
	{
		// First pass
		for (u8 y = 0; y < 4; y += 2)
		{
			s32 tmp1 = pb[stride * y + (x + 1)];
			s32 tmp2 = pb[stride * (y + 1) + x];
			s32 tmp3 = pb[stride * (y + 1) + (x + 1)];

			pb[stride * y + (x + 1)] += pb[stride * (y + 4) + (x + 1)];
			pb[stride * (y + 1) + x] += pb[stride * (y + 1 + 4) + x];
			pb[stride * (y + 1) + (x + 1)] += pb[stride * (y + 1 + 4) + (x + 1)];

			pb[stride * (y + 4) + (x + 1)] = tmp1 - pb[stride * (y + 4) + (x + 1)];
			pb[stride * (y + 1 + 4) + x] = tmp2 - pb[stride * (y + 1 + 4) + x];
			pb[stride * (y + 1 + 4) + (x + 1)] = tmp3 - pb[stride * (y + 1 + 4) + (x + 1)];
		}

		// Second pass
		for (u8 y = 0; y < 8; y += 4)
		{
			s32 tmp1 = pb[stride * y + (x + 1)];
			s32 tmp2 = pb[stride * (y + 1) + x];
			s32 tmp3 = pb[stride * (y + 1) + (x + 1)];

			pb[stride * y + (x + 1)] += pb[stride * (y + 2) + (x + 1)];
			pb[stride * (y + 1) + x] += pb[stride * (y + 1 + 2) + x];
			pb[stride * (y + 1) + (x + 1)] += pb[stride * (y + 1 + 2) + (x + 1)];

			pb[stride * (y + 2) + (x + 1)] = tmp1 - pb[stride * (y + 2) + (x + 1)];
			pb[stride * (y + 1 + 2) + x] = tmp2 - pb[stride * (y + 1 + 2) + x];
			pb[stride * (y + 1 + 2) + (x + 1)] = tmp3 - pb[stride * (y + 1 + 2) + (x + 1)];
		}

		// Metric calculation
		for (u8 y = 0; y < 8; y += 2)
		{
			s32 tmp1 = pb[stride * y + (x + 1)] < 0 ? -pb[stride * y + (x + 1)] : pb[stride * y + (x + 1)];
			s32 tmp2 = pb[stride * (y + 1) + x] < 0 ? -pb[stride * (y + 1) + x] : pb[stride * (y + 1) + x];
			s32 tmp3 = pb[stride * (y + 1) + (x + 1)] < 0 ? -pb[stride * (y + 1) + (x + 1)] : pb[stride * (y + 1) + (x + 1)];

			metric += (u32)(tmp1 + tmp2 + tmp3);
		}
	}

	return metric;
}

//void predSearch(cpStruct *cp)
//{
//	// Perform predictions on the luma component
//	predSearchLuma(pu);
//
//	// Perform predictions on the chroma components
//	if (pu.chromaSub != CHROMA_400)
//	{
//		predSearchChroma(pu);
//	}
//}
//
//void predSearchLuma(puStruct &pu)
//{
//	// Allocate space for the prediction block
//	static s16 pb[81];
//
//	// Process the top-left prediction block
//	for (u8 i = 1; i <= 8; i++)
//	{
//		memcpy(&pb[9 * i + 1], &pu.cu->cbLuma[pu.cu->stride * (i - 1)], 8 * sizeof(s16));
//	}
//	pb[8] = pb[17];
//	pb[72] = pb[73];
//	pb[0] = (pb[17] + pb[73] + 1) / 2;
//	for (u8 i = 1; i < 8; i++)
//	{
//		pb[i] = ((8 - i) * pb[0] + i * pb[8] + 4) / 8;
//		pb[9 * i] = ((8 - i) * pb[0] + i * pb[72] + 4) / 8;
//	}
//	predSearchBlock(pb, &pu.modeLuma[0], &pu.scanLuma[0]);
//	for (u8 i = 1; i < 8; i++)
//	{
//		memcpy(&pu.cu->cbLuma[pu.cu->stride * (i - 1)], &pb[9 * i + 1], 7 * sizeof(s16));
//	}
//
//	// Process the top row of prediction blocks
//	for (u8 x = 8; x < pu.cu->width; x += 8)
//	{
//		for (u8 i = 1; i <= 8; i++)
//		{
//			memcpy(&pb[9 * i], &pu.cu->cbLuma[pu.cu->stride * (i - 1) + x - 1], 9 * sizeof(s16));
//		}
//		pb[0] = pb[9];
//		pb[8] = pb[17];
//		for (u8 i = 1; i < 8; i++)
//		{
//			pb[i] = ((8 - i) * pb[0] + i * pb[8] + 4) / 8;
//		}
//		predSearchBlock(pb, &pu.modeLuma[x / 8], &pu.scanLuma[x / 8]);
//		for (u8 i = 1; i < 8; i++)
//		{
//			memcpy(&pu.cu->cbLuma[pu.cu->stride * (i - 1) + x], &pb[9 * i + 1], 7 * sizeof(s16));
//		}
//	}
//
//	// Process the remaining rows of prediction blocks
//	for (u8 y = 8; y < pu.cu->height; y += 8)
//	{
//		// Process the left-most prediction block of the current row
//		for (u8 i = 0; i <= 8; i++)
//		{
//			memcpy(&pb[9 * i + 1], &pu.cu->cbLuma[pu.cu->stride * (y + i - 1)], 8 * sizeof(s16));
//		}
//		pb[0] = pb[1];
//		pb[72] = pb[73];
//		for (u8 i = 1; i < 8; i++)
//		{
//			pb[9 * i] = ((8 - i) * pb[0] + i * pb[72] + 4) / 8;
//		}
//		predSearchBlock(pb, &pu.modeLuma[pu.stride * (y / 8)], &pu.scanLuma[pu.stride * (y / 8)]);
//		for (u8 i = 1; i < 8; i++)
//		{
//			memcpy(&pu.cu->cbLuma[pu.cu->stride * (y + i - 1)], &pb[9 * i + 1], 7 * sizeof(s16));
//		}
//
//		// Process the remaining blocks of the current row
//		for (u8 x = 8; x < pu.cu->width; x += 8)
//		{
//			for (u8 i = 0; i <= 8; i++)
//			{
//				memcpy(&pb[9 * i], &pu.cu->cbLuma[pu.cu->stride * (y + i - 1) + x - 1], 9 * sizeof(s16));
//			}
//			predSearchBlock(pb, &pu.modeLuma[pu.stride * (y / 8) + x / 8], &pu.scanLuma[pu.stride * (y / 8) + x / 8]);
//			for (u8 i = 1; i < 8; i++)
//			{
//				memcpy(&pu.cu->cbLuma[pu.cu->stride * (y + i - 1) + x], &pb[9 * i + 1], 7 * sizeof(s16));
//			}
//		}
//	}
//}
//
//void predSearchChroma(puStruct &pu)
//{
//	// Allocate space for the prediction block
//	static s16 pb1[81], pb2[81];
//
//	// Determine the size factors based on chroma subsampling
//	const u8 widthDiv = (pu.chromaSub == CHROMA_444) ? 1 : 2;
//	const u8 heightDiv = (pu.chromaSub == CHROMA_444) ? 1 : 2;
//
//	// Process the top-left prediction block
//	for (u8 i = 1; i <= 8; i++)
//	{
//		memcpy(&pb1[9 * i + 1], &pu.cu->cbChroma1[pu.cu->stride * (i - 1)], 8 * sizeof(s16));
//		memcpy(&pb2[9 * i + 1], &pu.cu->cbChroma2[pu.cu->stride * (i - 1)], 8 * sizeof(s16));
//	}
//	pb1[8] = pb1[17];
//	pb2[8] = pb2[17];
//	pb1[72] = pb1[73];
//	pb2[72] = pb2[73];
//	pb1[0] = (pb1[17] + pb1[73] + 1) / 2;
//	pb2[0] = (pb2[17] + pb2[73] + 1) / 2;
//	for (u8 i = 1; i < 8; i++)
//	{
//		pb1[i] = ((8 - i) * pb1[0] + i * pb1[8] + 4) / 8;
//		pb2[i] = ((8 - i) * pb2[0] + i * pb2[8] + 4) / 8;
//		pb1[9 * i] = ((8 - i) * pb1[0] + i * pb1[72] + 4) / 8;
//		pb2[9 * i] = ((8 - i) * pb2[0] + i * pb2[72] + 4) / 8;
//	}
//	predSearchBlock(pb1, pb2, &pu.modeLuma[0], &pu.modeChroma[0], &pu.scanChroma[0]);
//	for (u8 i = 1; i < 8; i++)
//	{
//		memcpy(&pu.cu->cbChroma1[pu.cu->stride * (i - 1)], &pb1[9 * i + 1], 7 * sizeof(s16));
//		memcpy(&pu.cu->cbChroma2[pu.cu->stride * (i - 1)], &pb2[9 * i + 1], 7 * sizeof(s16));
//	}
//
//	// Process the top row of prediction blocks
//	for (u8 x = 8; x < pu.cu->width / widthDiv; x += 8)
//	{
//		for (u8 i = 1; i <= 8; i++)
//		{
//			memcpy(&pb1[9 * i], &pu.cu->cbChroma1[pu.cu->stride * (i - 1) + x - 1], 9 * sizeof(s16));
//			memcpy(&pb2[9 * i], &pu.cu->cbChroma2[pu.cu->stride * (i - 1) + x - 1], 9 * sizeof(s16));
//		}
//		pb1[0] = pb1[9];
//		pb2[0] = pb2[9];
//		pb1[8] = pb1[17];
//		pb2[8] = pb2[17];
//		for (u8 i = 1; i < 8; i++)
//		{
//			pb1[i] = ((8 - i) * pb1[0] + i * pb1[8] + 4) / 8;
//			pb2[i] = ((8 - i) * pb2[0] + i * pb2[8] + 4) / 8;
//		}
//		predSearchBlock(pb1, pb2, &pu.modeLuma[x / 4], &pu.modeChroma[x / 8], &pu.scanChroma[x / 8]);
//		for (u8 i = 1; i < 8; i++)
//		{
//			memcpy(&pu.cu->cbChroma1[pu.cu->stride * (i - 1) + x], &pb1[9 * i + 1], 7 * sizeof(s16));
//			memcpy(&pu.cu->cbChroma2[pu.cu->stride * (i - 1) + x], &pb2[9 * i + 1], 7 * sizeof(s16));
//		}
//	}
//
//	// Process the remaining rows of prediction blocks
//	for (u8 y = 8; y < pu.cu->height / heightDiv; y += 8)
//	{
//		// Process the left-most prediction block of the current row
//		for (u8 i = 0; i <= 8; i++)
//		{
//			memcpy(&pb1[9 * i + 1], &pu.cu->cbChroma1[pu.cu->stride * (y + i - 1)], 8 * sizeof(s16));
//			memcpy(&pb2[9 * i + 1], &pu.cu->cbChroma2[pu.cu->stride * (y + i - 1)], 8 * sizeof(s16));
//		}
//		pb1[0] = pb1[1];
//		pb2[0] = pb2[1];
//		pb1[72] = pb1[73];
//		pb2[72] = pb2[73];
//		for (u8 i = 1; i < 8; i++)
//		{
//			pb1[9 * i] = ((8 - i) * pb1[0] + i * pb1[72] + 4) / 8;
//			pb2[9 * i] = ((8 - i) * pb2[0] + i * pb2[72] + 4) / 8;
//		}
//		predSearchBlock(pb1, pb2, &pu.modeLuma[pu.stride * (y / 4)], &pu.modeChroma[pu.stride * (y / 8)], &pu.scanChroma[pu.stride * (y / 8)]);
//		for (u8 i = 1; i < 8; i++)
//		{
//			memcpy(&pu.cu->cbChroma1[pu.cu->stride * (y + i - 1)], &pb1[9 * i + 1], 7 * sizeof(s16));
//			memcpy(&pu.cu->cbChroma2[pu.cu->stride * (y + i - 1)], &pb2[9 * i + 1], 7 * sizeof(s16));
//		}
//
//		// Process the remaining blocks of the current row
//		for (u8 x = 8; x < pu.cu->width / widthDiv; x += 8)
//		{
//			for (u8 i = 0; i <= 8; i++)
//			{
//				memcpy(&pb1[9 * i], &pu.cu->cbChroma1[pu.cu->stride * (y + i - 1) + x - 1], 9 * sizeof(s16));
//				memcpy(&pb2[9 * i], &pu.cu->cbChroma2[pu.cu->stride * (y + i - 1) + x - 1], 9 * sizeof(s16));
//			}
//			predSearchBlock(pb1, pb2, &pu.modeLuma[pu.stride * (y / 4) + x / 4], &pu.modeChroma[pu.stride * (y / 8) + x / 8], &pu.scanChroma[pu.stride * (y / 8) + x / 8]);
//			for (u8 i = 1; i < 8; i++)
//			{
//				memcpy(&pu.cu->cbChroma1[pu.cu->stride * (y + i - 1) + x], &pb1[9 * i + 1], 7 * sizeof(s16));
//				memcpy(&pu.cu->cbChroma2[pu.cu->stride * (y + i - 1) + x], &pb2[9 * i + 1], 7 * sizeof(s16));
//			}
//		}
//	}
//}
//
//void predSearchBlock(s16 *pb, u8 *mode, ScanDir *scan)
//{
//	// Allocate space for the prediction residuals
//	static s16 pbEst[81];
//
//	// Find the best prediction mode
//	u16 minMetric = 16384;
//	for (u8 i = 0; i < 29; i++)
//	{
//		// Perform the prediction
//		memcpy(pbEst, pb, 81 * sizeof(s16));
//		predBlock(pbEst, i);
//
//		// Determine if current mode is optimal
//		u16 metric = hadamardMetric(pb, pbEst);
//		if (metric < minMetric)
//		{
//			minMetric = metric;
//			*mode = i;
//		}
//	}
//
//	// Assign the scan direction
//	if (*mode == 0)
//	{
//		*scan = SCAN_DIAG;
//	}
//	else
//	{
//		u8 dir = (*mode <= 16) ? (*mode - 1) : (*mode <= 24 ? (*mode - 17) : (*mode - 25));
//		dir = (dir + 1) % 8;
//
//		if (dir < 3)
//		{
//			*scan = SCAN_HORZ;
//		}
//		else if (dir == 3 || dir == 7)
//		{
//			*scan = SCAN_DIAG;
//		}
//		else
//		{
//			*scan = SCAN_VERT;
//		}
//	}
//
//	// Perform the final prediction
//	predBlock(pb, *mode);
//}
//
//void predSearchBlock(s16 *pbChroma1, s16 *pbChroma2, u8 *modeLuma, u8 *modeChroma, ScanDir *scan)
//{
//	// Allocate space for the prediction residuals
//	static s16 pbEstChroma1[81], pbEstChroma2[81];
//
//	// Generate a unique list of modes for searching
//	static u8 modeSearch[4];
//
//	// Copy the corresponding luma modes
//	modeSearch[0] = modeLuma[0];
//	modeSearch[1] = modeLuma[1];
//	modeSearch[2] = modeLuma[CU_SIZE / 8];
//	modeSearch[3] = modeLuma[CU_SIZE / 8 + 1];
//
//	// Remove duplicates from the search list
//	bool hasDC = (modeLuma[0] == 0);
//	bool hasVert = (modeLuma[0] == 21);
//	bool hasHorz = (modeLuma[0] == 17);
//	for (u8 i = 1; i < 4; i++)
//	{
//		bool isDuplicate = false;
//		for (u8 j = 0; j < i; j++)
//		{
//			isDuplicate |= (modeSearch[i] == modeSearch[j]);
//		}
//
//		if (isDuplicate)
//		{
//			if (!hasDC)
//			{
//				modeSearch[i] = 0;
//				hasDC = true;
//			}
//			else if (!hasVert)
//			{
//				modeSearch[i] = 21;
//				hasVert = true;
//			}
//			else if (!hasHorz)
//			{
//				modeSearch[i] = 17;
//				hasHorz = true;
//			}
//			else
//			{
//				modeSearch[i] = 25;
//			}
//		}
//	}
//
//	// Find the best prediction mode in the search list
//	u16 minMetric = 16384;
//	for (u8 i = 0; i < 4; i++)
//	{
//		// Perform the prediction
//		memcpy(pbEstChroma1, pbChroma1, 81 * sizeof(s16));
//		memcpy(pbEstChroma2, pbChroma2, 81 * sizeof(s16));
//		predBlock(pbEstChroma1, modeSearch[i]);
//		predBlock(pbEstChroma2, modeSearch[i]);
//
//		// Determine if current mode is optimal
//		u16 metric = hadamardMetric(pbChroma1, pbEstChroma1) + hadamardMetric(pbChroma2, pbEstChroma2);
//		if (metric < minMetric)
//		{
//			minMetric = metric;
//			*modeChroma = i;
//		}
//	}
//
//	// Assign the scan direction
//	if (modeSearch[*modeChroma] == 0)
//	{
//		*scan = SCAN_DIAG;
//	}
//	else
//	{
//		u8 dir = (modeSearch[*modeChroma] <= 16) ? (modeSearch[*modeChroma] - 1) : (modeSearch[*modeChroma] <= 24 ? (modeSearch[*modeChroma] - 17) : (modeSearch[*modeChroma] - 25));
//		dir = (dir + 1) % 8;
//
//		if (dir < 3)
//		{
//			*scan = SCAN_HORZ;
//		}
//		else if (dir == 3 || dir == 7)
//		{
//			*scan = SCAN_DIAG;
//		}
//		else
//		{
//			*scan = SCAN_VERT;
//		}
//	}
//
//	// Perform the final prediction
//	predBlock(pbChroma1, modeSearch[*modeChroma]);
//	predBlock(pbChroma2, modeSearch[*modeChroma]);
//}
//
//void predBlock(s16 *pb, u8 mode)
//{
//	// Pass to the appropriate handler
//	if (mode == 0)
//	{
//		predDC(pb);
//	}
//	else if (mode <= 16)
//	{
//		predUPM(pb, mode - 1);
//	}
//	else if (mode <= 24)
//	{
//		predLPM(pb, mode - 17);
//	}
//	else
//	{
//		predBPM(pb, mode - 24);
//	}
//}
//
//void predBlock(s16 *pb1, s16 *pb2, u8 *modeLuma, u8 modeChroma, ScanDir *scanChroma)
//{
//	// Generate a unique list of modes for searching
//	static u8 modeSearch[4];
//
//	// Copy the corresponding luma modes
//	modeSearch[0] = modeLuma[0];
//	modeSearch[1] = modeLuma[1];
//	modeSearch[2] = modeLuma[CU_SIZE / 8];
//	modeSearch[3] = modeLuma[CU_SIZE / 8 + 1];
//
//	// Remove duplicates from the search list
//	bool hasDC = (modeLuma[0] == 0);
//	bool hasVert = (modeLuma[0] == 21);
//	bool hasHorz = (modeLuma[0] == 17);
//	for (u8 i = 1; i < 4; i++)
//	{
//		bool isDuplicate = false;
//		for (u8 j = 0; j < i; j++)
//		{
//			isDuplicate |= (modeSearch[i] == modeSearch[j]);
//		}
//
//		if (isDuplicate)
//		{
//			if (!hasDC)
//			{
//				modeSearch[i] = 0;
//				hasDC = true;
//			}
//			else if (!hasVert)
//			{
//				modeSearch[i] = 21;
//				hasVert = true;
//			}
//			else if (!hasHorz)
//			{
//				modeSearch[i] = 17;
//				hasHorz = true;
//			}
//			else
//			{
//				modeSearch[i] = 25;
//			}
//		}
//	}
//
//	predBlock(pb1, modeSearch[modeChroma]);
//	predBlock(pb2, modeSearch[modeChroma]);
//
//	// Assign the scan direction
//	if (modeSearch[modeChroma] == 0)
//	{
//		*scanChroma = SCAN_DIAG;
//	}
//	else
//	{
//		u8 dir = (modeSearch[modeChroma] <= 16) ? (modeSearch[modeChroma] - 1) : (modeSearch[modeChroma] <= 24 ? (modeSearch[modeChroma] - 17) : (modeSearch[modeChroma] - 25));
//		dir = (dir + 1) % 8;
//
//		if (dir < 3)
//		{
//			*scanChroma = SCAN_HORZ;
//		}
//		else if (dir == 3 || dir == 7)
//		{
//			*scanChroma = SCAN_DIAG;
//		}
//		else
//		{
//			*scanChroma = SCAN_VERT;
//		}
//	}
//}
//
//
//void predUPM(s16 *pb, u8 dir)
//{
//	// Define a lookup for angular factors (128 * tan(angle) and 128 / tan(angle))
//	static const s16 mulFact[8] = { 0, 53, 128, 309, 0, -309, -128, -53 };
//	static const s16 divFact[8] = { 0, 309, 128, 53, 0, -53, -128, -309 };
//
//	// Determine the current factors
//	const s16 mf = mulFact[dir % 8];
//	const s16 df = divFact[dir % 8];
//
//	// Handle each quadrant and the horizontal and vertical cases separately
//	dir = dir % 16;
//	if (dir == 0)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				pb[9 * y + x] = pb[y * 9 + 8];
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 9] = (3 * pb[x + 9] + pb[x] + 2) / 4; // Top border
//			pb[x + 63] = (3 * pb[x + 63] + pb[x + 72] + 2) / 4; // Bottom border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 1] = (3 * pb[9 * y + 1] + pb[9 * y] + 2) / 4; // Left border
//		}
//	}
//	else if (dir < 4)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				// Calculate the reference positions
//				const s16 posRight = 128 * y + ((8 - x) * mf);
//				const s16 posBottom = 128 * x + ((8 - y) * df);
//
//				// Determine the source for the current pixel
//				if (posRight >= 0 && posRight <= 1024)
//				{
//					const s16 w = posRight % 128;
//					const s16 ref1 = pb[9 * (posRight / 128) + 8];
//					const s16 ref2 = pb[9 * (posRight / 128 + 1) + 8];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//				else
//				{
//					const s16 w = posBottom % 128;
//					const s16 ref1 = pb[posBottom / 128 + 72];
//					const s16 ref2 = pb[posBottom / 128 + 73];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 9] = (3 * pb[x + 9] + pb[x] + 2) / 4; // Top border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 1] = (3 * pb[9 * y + 1] + pb[9 * y] + 2) / 4; // Left border
//		}
//	}
//	else if (dir == 4)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				pb[9 * y + x] = pb[72 + x];
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 9] = (3 * pb[x + 9] + pb[x] + 2) / 4; // Top border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 1] = (3 * pb[9 * y + 1] + pb[9 * y] + 2) / 4; // Left border
//			pb[9 * y + 7] = (3 * pb[9 * y + 7] + pb[9 * y + 8] + 2) / 4; // Right border
//		}
//	}
//	else if (dir < 8)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				// Calculate the reference positions
//				const s16 posLeft = 128 * y - (x * mf);
//				const s16 posBottom = 128 * x + ((8 - y) * df);
//
//				// Determine the source for the current pixel
//				if (posLeft >= 0 && posLeft <= 1024)
//				{
//					const s16 w = posLeft % 128;
//					const s16 ref1 = pb[9 * (posLeft / 128)];
//					const s16 ref2 = pb[9 * (posLeft / 128 + 1)];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//				else
//				{
//					const s16 w = posBottom % 128;
//					const s16 ref1 = pb[posBottom / 128 + 72];
//					const s16 ref2 = pb[posBottom / 128 + 73];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 9] = (3 * pb[x + 9] + pb[x] + 2) / 4; // Top border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 7] = (3 * pb[9 * y + 7] + pb[9 * y + 8] + 2) / 4; // Right border
//		}
//	}
//	else if (dir == 8)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				pb[9 * y + x] = pb[y * 9];
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 9] = (3 * pb[x + 9] + pb[x] + 2) / 4; // Top border
//			pb[x + 63] = (3 * pb[x + 63] + pb[x + 72] + 2) / 4; // Bottom border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 7] = (3 * pb[9 * y + 7] + pb[9 * y + 8] + 2) / 4; // Right border
//		}
//	}
//	else if (dir < 12)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				// Calculate the reference positions
//				const s16 posLeft = 128 * y - (x * mf);
//				const s16 posTop = 128 * x - (y * df);
//
//				// Determine the source for the current pixel
//				if (posLeft >= 0 && posLeft <= 1024)
//				{
//					const s16 w = posLeft % 128;
//					const s16 ref1 = pb[9 * (posLeft / 128)];
//					const s16 ref2 = pb[9 * (posLeft / 128 + 1)];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//				else
//				{
//					const s16 w = posTop % 128;
//					const s16 ref1 = pb[posTop / 128];
//					const s16 ref2 = pb[posTop / 128 + 1];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 63] = (3 * pb[x + 63] + pb[x + 72] + 2) / 4; // Bottom border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 7] = (3 * pb[9 * y + 7] + pb[9 * y + 8] + 2) / 4; // Right border
//		}
//	}
//	else if (dir == 12)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				pb[9 * y + x] = pb[x];
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 63] = (3 * pb[x + 63] + pb[x + 72] + 2) / 4; // Bottom border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 1] = (3 * pb[9 * y + 1] + pb[9 * y] + 2) / 4; // Left border
//			pb[9 * y + 7] = (3 * pb[9 * y + 7] + pb[9 * y + 8] + 2) / 4; // Right border
//		}
//	}
//	else
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				// Calculate the reference positions
//				const s16 posRight = 128 * y + ((8 - x) * mf);
//				const s16 posTop = 128 * x - (y * df);
//
//				// Determine the source for the current pixel
//				if (posRight >= 0 && posRight <= 1024)
//				{
//					const s16 w = posRight % 128;
//					const s16 ref1 = pb[9 * (posRight / 128) + 8];
//					const s16 ref2 = pb[9 * (posRight / 128 + 1) + 8];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//				else
//				{
//					const s16 w = posTop % 128;
//					const s16 ref1 = pb[posTop / 128];
//					const s16 ref2 = pb[posTop / 128 + 1];
//
//					pb[9 * y + x] = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//				}
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 63] = (3 * pb[x + 63] + pb[x + 72] + 2) / 4; // Bottom border
//		}
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 1] = (3 * pb[9 * y + 1] + pb[9 * y] + 2) / 4; // Left border
//		}
//	}
//}
//
//void predLPM(s16 *pb, u8 dir)
//{
//	// Define a lookup for angular factors (128 * tan(angle) and 128 / tan(angle))
//	static const s16 mulFact[8] = { 0, 53, 128, 309, 0, -309, -128, -53 };
//	static const s16 divFact[8] = { 0, 309, 128, 53, 0, -53, -128, -309 };
//
//	// Determine the current factors
//	dir = dir % 8;
//	const s16 mf = mulFact[dir];
//	const s16 df = divFact[dir];
//
//	// Handle each quadrant and the horizontal and vertical cases separately
//	if (dir == 0)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				pb[9 * y + x] = ((8 - x) * pb[y * 9] + x * pb[y * 9 + 8] + 4) / 8;
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[x + 9] = (3 * pb[x + 9] + pb[x] + 2) / 4; // Top border
//			pb[x + 63] = (3 * pb[x + 63] + pb[x + 72] + 2) / 4; // Bottom border
//		}
//	}
//	else if (dir < 4)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				// Calculate the reference positions
//				const s16 posLeft = 128 * y - (x * mf);
//				const s16 posRight = 128 * y + ((8 - x) * mf);
//				const s16 posTop = 128 * x - (y * df);
//				const s16 posBottom = 128 * x + ((8 - y) * df);
//
//				// Determine the source for the current pixel
//				s16 refTopLeft, refBottomRight;
//				s16 distTopLeft, distBottomRight;
//
//				if (posLeft >= 0 && posLeft <= 1024)
//				{
//					const s16 w = posLeft % 128;
//					const s16 ref1 = pb[9 * (posLeft / 128)];
//					const s16 ref2 = pb[9 * (posLeft / 128 + 1)];
//
//					refTopLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopLeft = 8 * x + 8 * y - posLeft / 16;
//				}
//				else
//				{
//					const s16 w = posTop % 128;
//					const s16 ref1 = pb[posTop / 128];
//					const s16 ref2 = pb[posTop / 128 + 1];
//
//					refTopLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopLeft = 8 * x - posTop / 16 + 8 * y;
//				}
//
//				if (posRight >= 0 && posRight <= 1024)
//				{
//					const s16 w = posRight % 128;
//					const s16 ref1 = pb[9 * (posRight / 128) + 8];
//					const s16 ref2 = pb[9 * (posRight / 128 + 1) + 8];
//
//					refBottomRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomRight = 8 * (8 - x) + posRight / 16 - 8 * y;
//				}
//				else
//				{
//					const s16 w = posBottom % 128;
//					const s16 ref1 = pb[posBottom / 128 + 72];
//					const s16 ref2 = pb[posBottom / 128 + 73];
//
//					refBottomRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomRight = posBottom / 16 - 8 * x + 8 * (8 - y);
//				}
//
//				pb[9 * y + x] = (distBottomRight * refTopLeft + distTopLeft * refBottomRight) / (distBottomRight + distTopLeft);
//			}
//		}
//	}
//	else if (dir == 4)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				pb[9 * y + x] = ((8 - y) * pb[x] + y * pb[72 + x] + 4) / 8;
//			}
//		}
//
//		// Apply weak filtering to the borders
//		for (u8 y = 1; y < 8; y++)
//		{
//			pb[9 * y + 1] = (3 * pb[9 * y + 1] + pb[9 * y] + 2) / 4; // Left border
//			pb[9 * y + 7] = (3 * pb[9 * y + 7] + pb[9 * y + 8] + 2) / 4; // Right border
//		}
//	}
//	else
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				// Calculate the reference positions
//				const s16 posLeft = 128 * y - (x * mf);
//				const s16 posRight = 128 * y + ((8 - x) * mf);
//				const s16 posTop = 128 * x - (y * df);
//				const s16 posBottom = 128 * x + ((8 - y) * df);
//
//				// Determine the source for the current pixel
//				s16 refTopRight, refBottomLeft;
//				s16 distTopRight, distBottomLeft;
//
//				if (posRight >= 0 && posRight <= 1024)
//				{
//					const s16 w = posRight % 128;
//					const s16 ref1 = pb[9 * (posRight / 128) + 8];
//					const s16 ref2 = pb[9 * (posRight / 128 + 1) + 8];
//
//					refTopRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopRight = 8 * (8 - x) + 8 * y - posRight / 16;
//				}
//				else
//				{
//					const s16 w = posTop % 128;
//					const s16 ref1 = pb[posTop / 128];
//					const s16 ref2 = pb[posTop / 128 + 1];
//
//					refTopRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopRight = posTop / 16 - 8 * x + 8 * y;
//				}
//
//				if (posLeft >= 0 && posLeft <= 1024)
//				{
//					const s16 w = posLeft % 128;
//					const s16 ref1 = pb[9 * (posLeft / 128)];
//					const s16 ref2 = pb[9 * (posLeft / 128 + 1)];
//
//					refBottomLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomLeft = 8 * x + posLeft / 16 - 8 * y;
//				}
//				else
//				{
//					const s16 w = posBottom % 128;
//					const s16 ref1 = pb[posBottom / 128 + 72];
//					const s16 ref2 = pb[posBottom / 128 + 73];
//
//					refBottomLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomLeft = 8 * x - posBottom / 16 + 8 * (8 - y);
//				}
//
//				pb[9 * y + x] = (distTopRight * refBottomLeft + distBottomLeft * refTopRight) / (distTopRight + distBottomLeft);
//			}
//		}
//	}
//}
//
//void predBPM(s16 *pb, u8 dir)
//{
//	// Define a lookup for angular factors (128 * tan(angle) and 128 / tan(angle))
//	static const s16 mulFact[8] = { 0, 53, 128, 309 };
//	static const s16 divFact[8] = { 0, 309, 128, 53 };
//
//	// Determine the current factors
//	dir = dir % 4;
//	const s16 mf = mulFact[dir];
//	const s16 df = divFact[dir];
//
//	// Handle each quadrant and the horizontal and vertical cases separately
//	if (dir == 0)
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				pb[9 * y + x] = ((8 - x) * pb[y * 9] + x * pb[y * 9 + 8] + (8 - y) * pb[x] + y * pb[72 + x] + 8) / 16;
//			}
//		}
//	}
//	else
//	{
//		// Predict the pixels
//		for (u8 y = 1; y < 8; y++)
//		{
//			for (u8 x = 1; x < 8; x++)
//			{
//				// Calculate the reference positions
//				const s16 posLeft = 128 * y - (x * mf);
//				const s16 posRight = 128 * y + ((8 - x) * mf);
//				const s16 posTop = 128 * x - (y * df);
//				const s16 posBottom = 128 * x + ((8 - y) * df);
//
//				const s16 posLeft2 = 128 * y + (x * df);
//				const s16 posRight2 = 128 * y - ((8 - x) * df);
//				const s16 posTop2 = 128 * x + (y * mf);
//				const s16 posBottom2 = 128 * x - ((8 - y) * mf);
//
//				// Determine the source for the current pixel
//				s16 refTopLeft, refBottomRight;
//				s16 distTopLeft, distBottomRight;
//
//				s16 refTopRight, refBottomLeft;
//				s16 distTopRight, distBottomLeft;
//
//				// Main references
//				if (posLeft >= 0 && posLeft <= 1024)
//				{
//					const s16 w = posLeft % 128;
//					const s16 ref1 = pb[9 * (posLeft / 128)];
//					const s16 ref2 = pb[9 * (posLeft / 128 + 1)];
//
//					refTopLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopLeft = 8 * x + 8 * y - posLeft / 16;
//				}
//				else
//				{
//					const s16 w = posTop % 128;
//					const s16 ref1 = pb[posTop / 128];
//					const s16 ref2 = pb[posTop / 128 + 1];
//
//					refTopLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopLeft = 8 * x - posTop / 16 + 8 * y;
//				}
//
//				if (posRight >= 0 && posRight <= 1024)
//				{
//					const s16 w = posRight % 128;
//					const s16 ref1 = pb[9 * (posRight / 128) + 8];
//					const s16 ref2 = pb[9 * (posRight / 128 + 1) + 8];
//
//					refBottomRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomRight = 8 * (8 - x) + posRight / 16 - 8 * y;
//				}
//				else
//				{
//					const s16 w = posBottom % 128;
//					const s16 ref1 = pb[posBottom / 128 + 72];
//					const s16 ref2 = pb[posBottom / 128 + 73];
//
//					refBottomRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomRight = posBottom / 16 - 8 * x + 8 * (8 - y);
//				}
//
//				// Perpendicular references
//				if (posRight2 >= 0 && posRight2 <= 1024)
//				{
//					const s16 w = posRight2 % 128;
//					const s16 ref1 = pb[9 * (posRight2 / 128) + 8];
//					const s16 ref2 = pb[9 * (posRight2 / 128 + 1) + 8];
//
//					refTopRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopRight = 8 * (8 - x) + 8 * y - posRight2 / 16;
//				}
//				else
//				{
//					const s16 w = posTop2 % 128;
//					const s16 ref1 = pb[posTop2 / 128];
//					const s16 ref2 = pb[posTop2 / 128 + 1];
//
//					refTopRight = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distTopRight = posTop2 / 16 - 8 * x + 8 * y;
//				}
//
//				if (posLeft2 >= 0 && posLeft2 <= 1024)
//				{
//					const s16 w = posLeft2 % 128;
//					const s16 ref1 = pb[9 * (posLeft2 / 128)];
//					const s16 ref2 = pb[9 * (posLeft2 / 128 + 1)];
//
//					refBottomLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomLeft = 8 * x + posLeft2 / 16 - 8 * y;
//				}
//				else
//				{
//					const s16 w = posBottom2 % 128;
//					const s16 ref1 = pb[posBottom2 / 128 + 72];
//					const s16 ref2 = pb[posBottom2 / 128 + 73];
//
//					refBottomLeft = ((128 - w) * ref1 + w * ref2 + 64) / 128;
//					distBottomLeft = 8 * x - posBottom2 / 16 + 8 * (8 - y);
//				}
//
//				pb[9 * y + x] = (distBottomRight * refTopLeft + distTopLeft * refBottomRight) / (distBottomRight + distTopLeft);
//				pb[9 * y + x] += (distTopRight * refBottomLeft + distBottomLeft * refTopRight) / (distTopRight + distBottomLeft);
//				pb[9 * y + x] = (pb[9 * y + x] + 1) / 2;
//			}
//		}
//	}
//}
//
//void predNN(s16 *pb)
//{
//	// Calculate the DC value of the known pixels
//	s16 dcVal = 0;
//	for (u8 y = 1; y < 8; y++)
//	{
//		dcVal += (pb[9 * y] + pb[9 * y + 8]);
//	}
//	for (u8 x = 0; x <= 8; x++)
//	{
//		dcVal += (pb[x] + pb[x + 72]);
//	}
//	dcVal /= 32;
//
//	// Fill all pixels
//	for (u8 y = 1; y < 8; y++)
//	{
//		for (u8 x = 1; x < 8; x++)
//		{
//			pb[9 * y + x] = dcVal;
//		}
//	}
//
//	// Apply weak filtering to the borders
//	for (u8 x = 1; x < 8; x++)
//	{
//		pb[x + 9] = (3 * pb[x + 9] + pb[x] + 2) / 4; // Top border
//		pb[x + 63] = (3 * pb[x + 63] + pb[x + 72] + 2) / 4; // Bottom border
//	}
//	for (u8 y = 1; y < 8; y++)
//	{
//		pb[9 * y + 1] = (3 * pb[9 * y + 1] + pb[9 * y] + 2) / 4; // Left border
//		pb[9 * y + 7] = (3 * pb[9 * y + 7] + pb[9 * y + 8] + 2) / 4; // Right border
//	}
//}
//
