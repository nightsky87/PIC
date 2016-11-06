#include "PICComBasis.h"
#include "PICComBasisLookup.h"
#include "PICComQuant.h"
#include <cstring>

void basisSearch(cuStruct *cu, const u8 sWidth, const u8 qp)
{
	// Assign the height of all CU scales to be equal to the width
	const u8 sHeight = sWidth;

	// Calculate the horizontal and vertical offsets
	const u8 hOffset = sWidth / 2;
	const u8 vOffset = sHeight / 2;

	// Find a suitable basis for every 8x8 block
	for (u8 y = 0; y < sHeight / 2; y += 4)
	{
		for (u8 x = 0; x < sWidth / 2; x += 4)
		{
			// Find the l0 metric for residuals within the block
			static u32 metric[NUM_BASIS];
			memset(metric, 0, NUM_BASIS * sizeof(u32));
			for (u8 i = 0; i < 4; i++)
			{
				for (u8 j = 0; j < 4; j++)
				{
					for (u8 basis = 0; basis < NUM_BASIS; basis++)
					{
						const s32 *basisTable = &basis_weights[12 * basis];

						const s32 p1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)];
						const s32 p2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + hOffset)];
						const s32 p3 = cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j)];
						const s32 p4 = cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j + hOffset)];

						const s32 t1 = (basisTable[0] * p1 + basisTable[1] * p2 + basisTable[2] * p3 + basisTable[3] * p4 + 32) / 64;
						const s32 t2 = (basisTable[4] * p1 + basisTable[5] * p2 + basisTable[6] * p3 + basisTable[7] * p4 + 32) / 64;
						const s32 t3 = (basisTable[8] * p1 + basisTable[9] * p2 + basisTable[10] * p3 + basisTable[11] * p4 + 32) / 64;

						metric[basis] += t1 < 0 ? -t1 : t1;
						metric[basis] += t2 < 0 ? -t2 : t2;
						metric[basis] += t3 < 0 ? -t3 : t3;
					}
				}
			}

			// Use a greedy minimization criteria to select the best basis
			u8 bestBasis = 0;
			u32 bestMetric = metric[0];
			for (u8 basis = 1; basis < NUM_BASIS; basis++)
			{
				if (metric[basis] < bestMetric)
				{
					bestBasis = basis;
					bestMetric = metric[basis];
				}
			}
			cu->basisSelect[cu->basisIndex] = bestBasis;
			cu->basisIndex++;
			const s32 *basisTable = &basis_weights[12 * bestBasis];

			// Perform the final transform for the given block
			for (u8 i = 0; i < 4; i++)
			{
				for (u8 j = 0; j < 4; j++)
				{
					const s32 p1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)];
					const s32 p2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + hOffset)];
					const s32 p3 = cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j)];
					const s32 p4 = cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j + hOffset)];

					cu->pLuma[CU_SIZE * (y + i) + (x + j)] = 0;
					cu->pLuma[CU_SIZE * (y + i) + (x + j + hOffset)] = (pel)quantVal((basisTable[0] * p1 + basisTable[1] * p2 + basisTable[2] * p3 + basisTable[3] * p4 + 32) / 64, qp);
					cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j)] = (pel)quantVal((basisTable[4] * p1 + basisTable[5] * p2 + basisTable[6] * p3 + basisTable[7] * p4 + 32) / 64, qp);
					cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j + hOffset)] = (pel)quantVal((basisTable[8] * p1 + basisTable[9] * p2 + basisTable[10] * p3 + basisTable[11] * p4 + 32) / 64, qp);
				}
			}
		}
	}
}

void basisInverse(cuStruct *cu, const u8 sWidth)
{
	// Allocate a static buffer for the pixels
	static pel buf[CU_SIZE * CU_SIZE];

	// Assign the height of all CU scales to be equal to the width
	const u8 sHeight = sWidth;

	// Calculate the horizontal and vertical offsets
	const u8 hOffset = sWidth / 2;
	const u8 vOffset = sHeight / 2;

	// Reconstruct each 8x8 block
	for (s8 y = sHeight / 2 - 4; y >= 0; y -= 4)
	{
		for (s8 x = sWidth / 2 - 4; x >= 0; x -= 4)
		{
			// Determine the basis for reconstruction
			cu->basisIndex--;
			u8 basis = cu->basisSelect[cu->basisIndex];
			const s32 *basisTable = &basis_weights[12 * basis];

			// Perform the final transform for the given block
			for (u8 i = 0; i < 4; i++)
			{
				for (u8 j = 0; j < 4; j++)
				{
					const s32 p1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)];
					const s32 p2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + hOffset)];
					const s32 p3 = cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j)];
					const s32 p4 = cu->pLuma[CU_SIZE * (y + i + vOffset) + (x + j + hOffset)];

					buf[CU_SIZE * 2 * (y + i) + 2 * (x + j)] = (pel)((64 * p1 + basisTable[0] * p2 + basisTable[4] * p3 + basisTable[8] * p4 + 128) / 256);
					buf[CU_SIZE * 2 * (y + i) + 2 * (x + j) + 1] = (pel)((64 * p1 + basisTable[1] * p2 + basisTable[5] * p3 + basisTable[9] * p4 + 128) / 256);
					buf[CU_SIZE * 2 * (y + i) + 2 * (x + j) + CU_SIZE] = (pel)((64 * p1 + basisTable[2] * p2 + basisTable[6] * p3 + basisTable[10] * p4 + 128) / 256);
					buf[CU_SIZE * 2 * (y + i) + 2 * (x + j) + CU_SIZE + 1] = (pel)((64 * p1 + basisTable[3] * p2 + basisTable[7] * p3 + basisTable[11] * p4 + 128) / 256);
				}
			}
		}
	}

	// Copy the new values to the CU
	for (u8 y = 0; y < sHeight; y++)
		memcpy(&cu->pLuma[CU_SIZE * y], &buf[CU_SIZE * y], sWidth * sizeof(pel));
}
