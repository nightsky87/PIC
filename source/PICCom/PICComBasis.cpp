#include "PICComBasis.h"
#include "PICComBasisLookup.h"
#include "PICComQuant.h"
#include <cstring>

void basisSearch(cuStruct *cu, const u8 sWidth, const u8 sHeight)
{
	// Allocate a buffer for the pixel sums
	const u16 bufWidth = CU_SIZE + 1;
	static pel sumLuma[bufWidth * bufWidth];

	// Generate the sums
	for (u8 y = 0; y < sHeight; y += 2)
	{
		for (u8 x = 0; x < sWidth; x += 2)
		{
			sumLuma[bufWidth * y + x] = cu->pLuma[CU_SIZE * y + x];
			sumLuma[bufWidth * y + x] += cu->pLuma[CU_SIZE * y + x + 1];
			sumLuma[bufWidth * y + x] += cu->pLuma[CU_SIZE * (y + 1) + x];
			sumLuma[bufWidth * y + x] += cu->pLuma[CU_SIZE * (y + 1) + x + 1];
		}
		sumLuma[bufWidth * y + sWidth] = sumLuma[bufWidth * y + sWidth - 4];
	}
	memcpy(&sumLuma[bufWidth * sHeight], &sumLuma[bufWidth * (sHeight - 4)], (sWidth + 1) * sizeof(pel));

	// Find a suitable basis for every 8x8 block
	u8 ind = cu->basisIndex;
	for (u8 y = 0; y < sHeight; y += 8)
	{
		for (u8 x = 0; x < sWidth; x += 8)
		{
			static u32 metric[NUM_BASIS];
			for (u8 basis = 0; basis < NUM_BASIS; basis++)
			{
				const s32 *basisTable = &basis_weights[12 * basis];
				const s32 *interpTable = &interp_weights[16 * basis];

				// Find the interpolation residuals mapped to the microbasis
				for (u8 i = 0; i < 8; i += 2)
				{
					for (u8 j = 0; j < 8; j += 2)
					{
						// Calculate the prediction errors
						const s32 e1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)] - (interpTable[0] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[1] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[2] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[3] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
						const s32 e2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + 1)] - (interpTable[4] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[5] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[6] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[7] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
						const s32 e3 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j)] - (interpTable[8] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[9] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[10] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[11] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
						const s32 e4 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j + 1)] - (interpTable[12] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[13] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[14] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[15] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;

						sumLuma[bufWidth * (y + i) + (x + j + 1)] = (basisTable[0] * e1 + basisTable[1] * e2 + basisTable[2] * e3 + basisTable[3] * e4 + 2048) / 4096;
						sumLuma[bufWidth * (y + i + 1) + (x + j)] = (basisTable[4] * e1 + basisTable[5] * e2 + basisTable[6] * e3 + basisTable[7] * e4 + 2048) / 4096;
						sumLuma[bufWidth * (y + i + 1) + (x + j + 1)] = (basisTable[8] * e1 + basisTable[9] * e2 + basisTable[10] * e3 + basisTable[11] * e4 + 2048) / 4096;
					}
				}

				// Perform the Hadamard transform on the microbasis coefficients
				metric[basis] = hadamardMetric(&sumLuma[bufWidth * y + x + 1], bufWidth);
				metric[basis] += hadamardMetric(&sumLuma[bufWidth * (y + 1) + x], bufWidth);
				metric[basis] += hadamardMetric(&sumLuma[bufWidth * (y + 1) + x + 1], bufWidth);
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
			cu->basisSelect[ind] = bestBasis;
			ind++;

			const s32 *basisTable = &basis_weights[12 * bestBasis];
			const s32 *interpTable = &interp_weights[16 * bestBasis];

			// Find the interpolation residuals mapped to the microbasis
			for (u8 i = 0; i < 8; i += 2)
			{
				for (u8 j = 0; j < 8; j += 2)
				{
					// Calculate the prediction errors
					const s32 e1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)] - (interpTable[0] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[1] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[2] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[3] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
					const s32 e2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + 1)] - (interpTable[4] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[5] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[6] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[7] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
					const s32 e3 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j)] - (interpTable[8] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[9] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[10] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[11] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
					const s32 e4 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j + 1)] - (interpTable[12] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[13] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[14] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[15] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;

					sumLuma[bufWidth * (y + i) + (x + j + 1)] = (basisTable[0] * e1 + basisTable[1] * e2 + basisTable[2] * e3 + basisTable[3] * e4 + 2048) / 4096;
					sumLuma[bufWidth * (y + i + 1) + (x + j)] = (basisTable[4] * e1 + basisTable[5] * e2 + basisTable[6] * e3 + basisTable[7] * e4 + 2048) / 4096;
					sumLuma[bufWidth * (y + i + 1) + (x + j + 1)] = (basisTable[8] * e1 + basisTable[9] * e2 + basisTable[10] * e3 + basisTable[11] * e4 + 2048) / 4096;
				}
			}
		}
	}

	// Copy the final coefficients to the CU structure
	for (u8 y = 0; y < sHeight; y += 2)
	{
		for (u8 x = 0; x < sWidth; x += 2)
		{
			cu->pLuma[CU_SIZE * y / 2 + x / 2] = sumLuma[bufWidth * y + x];
			cu->pLuma[CU_SIZE * y / 2 + x / 2 + sWidth / 2] = sumLuma[bufWidth * y + x + 1];
			cu->pLuma[CU_SIZE * (y / 2 + sHeight / 2) + x / 2] = sumLuma[bufWidth * (y + 1) + x];
			cu->pLuma[CU_SIZE * (y / 2 + sHeight / 2) + x / 2 + sWidth / 2] = sumLuma[bufWidth * (y + 1) + x + 1];
		}
	}

	// Adjust the basis index
	cu->basisIndex = ind;
}

void basisInverse(cuStruct *cu, const u8 sWidth, const u8 sHeight)
{
	// Allocate a buffer for the pixel sums
	const u16 bufWidth = CU_SIZE + 1;
	static pel sumLuma[bufWidth * bufWidth];

	// Copy the sums and the residual
	for (u8 y = 0; y < sHeight; y += 2)
	{
		for (u8 x = 0; x < sWidth; x += 2)
		{
			sumLuma[bufWidth * y + x] = cu->pLuma[CU_SIZE * y / 2 + x / 2];
			sumLuma[bufWidth * y + x + 1] = cu->pLuma[CU_SIZE * y / 2 + x / 2 + sWidth / 2];
			sumLuma[bufWidth * (y + 1) + x] = cu->pLuma[CU_SIZE * (y / 2 + sHeight / 2) + x / 2];
			sumLuma[bufWidth * (y + 1) + x + 1] = cu->pLuma[CU_SIZE * (y / 2 + sHeight / 2) + x / 2 + sWidth / 2];
		}
		sumLuma[bufWidth * y + sWidth] = sumLuma[bufWidth * y + sWidth - 4];
	}
	memcpy(&sumLuma[bufWidth * sHeight], &sumLuma[bufWidth * (sHeight - 4)], (sWidth + 1) * sizeof(pel));

	// Apply the chosen bases
	u8 ind = cu->basisIndex;
	for (s16 y = sHeight - 8; y >= 0; y -= 8)
	{
		for (s16 x = sWidth - 8; x >= 0; x -= 8)
		{
			ind--;
			const s32 *basisTable = &basis_weights[12 * cu->basisSelect[ind]];
			const s32 *interpTable = &interp_weights[16 * cu->basisSelect[ind]];

			for (u8 i = 0; i < 8; i += 2)
			{
				for (u8 j = 0; j < 8; j += 2)
				{
					// Calculate the prediction errors
					const s32 p1 = (interpTable[0] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[1] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[2] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[3] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
					const s32 p2 = (interpTable[4] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[5] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[6] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[7] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
					const s32 p3 = (interpTable[8] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[9] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[10] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[11] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;
					const s32 p4 = (interpTable[12] * sumLuma[bufWidth * (y + i) + (x + j)] + interpTable[13] * sumLuma[bufWidth * (y + i) + (x + j + 2)] + interpTable[14] * sumLuma[bufWidth * (y + i + 2) + (x + j)] + interpTable[15] * sumLuma[bufWidth * (y + i + 2) + (x + j + 2)] + 2048) / 4096;

					const s32 t1 = sumLuma[bufWidth * (y + i) + (x + j + 1)];
					const s32 t2 = sumLuma[bufWidth * (y + i + 1) + (x + j)];
					const s32 t3 = sumLuma[bufWidth * (y + i + 1) + (x + j + 1)];

					cu->pLuma[CU_SIZE * (y + i) + (x + j)] = p1 + (basisTable[0] * t1 + basisTable[4] * t2 + basisTable[8] * t3 + 2048) / 4096;
					cu->pLuma[CU_SIZE * (y + i) + (x + j + 1)] = p2 + (basisTable[1] * t1 + basisTable[5] * t2 + basisTable[9] * t3 + 2048) / 4096;
					cu->pLuma[CU_SIZE * (y + i + 1) + (x + j)] = p3 + (basisTable[2] * t1 + basisTable[6] * t2 + basisTable[10] * t3 + 2048) / 4096;
					cu->pLuma[CU_SIZE * (y + i + 1) + (x + j + 1)] = p4 + (basisTable[3] * t1 + basisTable[7] * t2 + basisTable[11] * t3 + 2048) / 4096;
				}
			}
		}
	}

	// Adjust the basis index
	cu->basisIndex = ind;
}

u32 hadamardMetric(pel *tb, const u8 stride)
{
	// Horizontal pass
	for (u8 y = 0; y < 8; y += 2)
	{
		pel t0 = tb[stride * y + 0];
		pel t1 = tb[stride * y + 2];
		tb[stride * y + 0] += tb[stride * y + 4];
		tb[stride * y + 2] += tb[stride * y + 6];
		tb[stride * y + 4] = t0 - tb[stride * y + 4];
		tb[stride * y + 6] = t1 - tb[stride * y + 6];

		t0 = tb[stride * y + 0];
		t1 = tb[stride * y + 4];
		tb[stride * y + 0] += tb[stride * y + 2];
		tb[stride * y + 4] += tb[stride * y + 6];
		tb[stride * y + 2] = t0 - tb[stride * y + 2];
		tb[stride * y + 6] = t1 - tb[stride * y + 6];
	}

	// Vertical pass
	for (u8 x = 0; x < 8; x += 2)
	{
		pel t0 = tb[stride * 0 + x];
		pel t1 = tb[stride * 2 + x];
		tb[stride * 0 + x] += tb[stride * 4 + x];
		tb[stride * 2 + x] += tb[stride * 6 + x];
		tb[stride * 4 + x] = t0 - tb[stride * 4 + x];
		tb[stride * 6 + x] = t1 - tb[stride * 6 + x];

		t0 = tb[stride * 0 + x];
		t1 = tb[stride * 4 + x];
		tb[stride * 0 + x] += tb[stride * 2 + x];
		tb[stride * 4 + x] += tb[stride * 6 + x];
		tb[stride * 2 + x] = t0 - tb[stride * 2 + x];
		tb[stride * 6 + x] = t1 - tb[stride * 6 + x];
	}

	// Metric calculation
	u32 metric = 0;
	for (u8 y = 0; y < 8; y += 2)
	{
		for (u8 x = 0; x < 8; x += 2)
		{
			metric += (tb[stride * y + x] < 0 ? -tb[stride * y + x] : tb[stride * y + x]);
		}
	}

	return metric;
}