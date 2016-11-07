#include "PICComBasis.h"
#include "PICComBasisLookup.h"
#include "PICComQuant.h"
#include <cstring>

void basisSearch(cuStruct *cu, const u8 sWidth, const u8 sHeight, const u8 qp)
{
	// Allocate a buffer for the transform coefficients
	static pel luma[CU_SIZE * CU_SIZE];
	static pel chroma1[CU_SIZE * CU_SIZE];
	static pel chroma2[CU_SIZE * CU_SIZE];

	// Calculate the offsets
	const u8 xoff = sWidth / 2;
	const u8 yoff = sHeight / 2;

	// Find a suitable basis for every 8x8 block
	u8 ind = cu->basisIndex;
	if (cu->chromaSub == CHROMA_400)
	{
		for (u8 y = 0; y < sHeight; y += 8)
		{
			for (u8 x = 0; x < sWidth; x += 8)
			{
				// Find the l0 metric for residuals within the block
				static u32 metric[NUM_BASIS];
				memset(metric, 0, NUM_BASIS * sizeof(s32));
				for (u8 basis = 0; basis < NUM_BASIS; basis++)
				{
					const s32 *basisTable = &basis_weights[12 * basis];

					for (u8 i = 0; i < 8; i += 2)
					{
						for (u8 j = 0; j < 8; j += 2)
						{
							const s32 p1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)];
							const s32 p2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + 1)];
							const s32 p3 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j)];
							const s32 p4 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j + 1)];

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
				cu->basisSelect[ind] = bestBasis;
				ind++;

				// Apply the chosen basis
				const s32 *basisTable = &basis_weights[12 * bestBasis];
				for (u8 i = 0; i < 8; i += 2)
				{
					for (u8 j = 0; j < 8; j += 2)
					{
						const s32 p1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)];
						const s32 p2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + 1)];
						const s32 p3 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j)];
						const s32 p4 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j + 1)];

						luma[CU_SIZE * (y + i) / 2 + (x + j) / 2] = 0;
						luma[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff] = quantVal((basisTable[0] * p1 + basisTable[1] * p2 + basisTable[2] * p3 + basisTable[3] * p4 + 128) / 256, qp);
						luma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2] = quantVal((basisTable[4] * p1 + basisTable[5] * p2 + basisTable[6] * p3 + basisTable[7] * p4 + 128) / 256, qp);
						luma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff] = quantVal((basisTable[8] * p1 + basisTable[9] * p2 + basisTable[10] * p3 + basisTable[11] * p4 + 128) / 256, qp);
					}
				}
			}
		}


		for (u8 y = 0; y < sHeight; y++)
			memcpy(&cu->pLuma[CU_SIZE * y], &luma[CU_SIZE * y], sWidth * sizeof(pel));
	}
	else
	{
		for (u8 y = 0; y < sHeight; y += 8)
		{
			for (u8 x = 0; x < sWidth; x += 8)
			{
				// Find the l0 metric for residuals within the block
				static u32 metric[NUM_BASIS];
				memset(metric, 0, NUM_BASIS * sizeof(s32));
				for (u8 basis = 0; basis < NUM_BASIS; basis++)
				{
					const s32 *basisTable = &basis_weights[12 * basis];

					for (u8 i = 0; i < 8; i += 2)
					{
						for (u8 j = 0; j < 8; j += 2)
						{
							const s32 p1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)];
							const s32 p2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + 1)];
							const s32 p3 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j)];
							const s32 p4 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j + 1)];

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
				cu->basisSelect[ind] = bestBasis;
				ind++;

				// Apply the chosen basis
				const s32 *basisTable = &basis_weights[12 * bestBasis];
				for (u8 i = 0; i < 8; i += 2)
				{
					for (u8 j = 0; j < 8; j += 2)
					{
						s32 p1 = cu->pLuma[CU_SIZE * (y + i) + (x + j)];
						s32 p2 = cu->pLuma[CU_SIZE * (y + i) + (x + j + 1)];
						s32 p3 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j)];
						s32 p4 = cu->pLuma[CU_SIZE * (y + i + 1) + (x + j + 1)];

						luma[CU_SIZE * (y + i) / 2 + (x + j) / 2] = 0;
						luma[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff] = quantVal((basisTable[0] * p1 + basisTable[1] * p2 + basisTable[2] * p3 + basisTable[3] * p4 + 128) / 256, qp);
						luma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2] = quantVal((basisTable[4] * p1 + basisTable[5] * p2 + basisTable[6] * p3 + basisTable[7] * p4 + 128) / 256, qp);
						luma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff] = quantVal((basisTable[8] * p1 + basisTable[9] * p2 + basisTable[10] * p3 + basisTable[11] * p4 + 128) / 256, qp);

						p1 = cu->pChroma1[CU_SIZE * (y + i) + (x + j)];
						p2 = cu->pChroma1[CU_SIZE * (y + i) + (x + j + 1)];
						p3 = cu->pChroma1[CU_SIZE * (y + i + 1) + (x + j)];
						p4 = cu->pChroma1[CU_SIZE * (y + i + 1) + (x + j + 1)];

						chroma1[CU_SIZE * (y + i) / 2 + (x + j) / 2] = 0;
						chroma1[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff] = quantVal((basisTable[0] * p1 + basisTable[1] * p2 + basisTable[2] * p3 + basisTable[3] * p4 + 128) / 256, qp);
						chroma1[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2] = quantVal((basisTable[4] * p1 + basisTable[5] * p2 + basisTable[6] * p3 + basisTable[7] * p4 + 128) / 256, qp);
						chroma1[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff] = quantVal((basisTable[8] * p1 + basisTable[9] * p2 + basisTable[10] * p3 + basisTable[11] * p4 + 128) / 256, qp);

						p1 = cu->pChroma2[CU_SIZE * (y + i) + (x + j)];
						p2 = cu->pChroma2[CU_SIZE * (y + i) + (x + j + 1)];
						p3 = cu->pChroma2[CU_SIZE * (y + i + 1) + (x + j)];
						p4 = cu->pChroma2[CU_SIZE * (y + i + 1) + (x + j + 1)];

						chroma2[CU_SIZE * (y + i) / 2 + (x + j) / 2] = 0;
						chroma2[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff] = quantVal((basisTable[0] * p1 + basisTable[1] * p2 + basisTable[2] * p3 + basisTable[3] * p4 + 128) / 256, qp);
						chroma2[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2] = quantVal((basisTable[4] * p1 + basisTable[5] * p2 + basisTable[6] * p3 + basisTable[7] * p4 + 128) / 256, qp);
						chroma2[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff] = quantVal((basisTable[8] * p1 + basisTable[9] * p2 + basisTable[10] * p3 + basisTable[11] * p4 + 128) / 256, qp);
					}
				}
			}
		}


		for (u8 y = 0; y < sHeight; y++)
		{
			memcpy(&cu->pLuma[CU_SIZE * y], &luma[CU_SIZE * y], sWidth * sizeof(pel));
			memcpy(&cu->pChroma1[CU_SIZE * y], &chroma1[CU_SIZE * y], sWidth * sizeof(pel));
			memcpy(&cu->pChroma2[CU_SIZE * y], &chroma2[CU_SIZE * y], sWidth * sizeof(pel));
		}
	}

	// Adjust the basis index
	cu->basisIndex = ind;
}

void basisInverse(cuStruct *cu, const u8 sWidth, const u8 sHeight)
{
	// Allocate a buffer for the transform coefficients
	static pel luma[CU_SIZE * CU_SIZE];
	static pel chroma1[CU_SIZE * CU_SIZE];
	static pel chroma2[CU_SIZE * CU_SIZE];

	// Calculate the horizontal and vertical offsets
	const u8 hOffset = sWidth / 2;
	const u8 vOffset = sHeight / 2;

	// Reconstruct each 8x8 block
	// Calculate the offsets
	const u8 xoff = sWidth / 2;
	const u8 yoff = sHeight / 2;

	// Apply the chosen bases
	u8 ind = cu->basisIndex;
	if (cu->chromaSub == CHROMA_400)
	{
		for (s16 y = sHeight - 8; y >= 0; y -= 8)
		{
			for (s16 x = sWidth - 8; x >= 0; x -= 8)
			{
				ind--;
				const s32 *basisTable = &basis_weights[12 * cu->basisSelect[ind]];

				for (u8 i = 0; i < 8; i += 2)
				{
					for (u8 j = 0; j < 8; j += 2)
					{
						const s32 p1 = cu->pLuma[CU_SIZE * (y + i) / 2 + (x + j) / 2];
						const s32 p2 = cu->pLuma[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff];
						const s32 p3 = cu->pLuma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2];
						const s32 p4 = cu->pLuma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff];

						luma[CU_SIZE * (y + i) + (x + j)] = (64 * p1 + basisTable[0] * p2 + basisTable[4] * p3 + basisTable[8] * p4 + 128) / 256;
						luma[CU_SIZE * (y + i) + (x + j + 1)] = (64 * p1 + basisTable[1] * p2 + basisTable[5] * p3 + basisTable[9] * p4 + 128) / 256;
						luma[CU_SIZE * (y + i + 1) + (x + j)] = (64 * p1 + basisTable[2] * p2 + basisTable[6] * p3 + basisTable[10] * p4 + 128) / 256;
						luma[CU_SIZE * (y + i + 1) + (x + j + 1)] = (64 * p1 + basisTable[3] * p2 + basisTable[7] * p3 + basisTable[11] * p4 + 128) / 256;
					}
				}
			}
		}

		for (u8 y = 0; y < sHeight; y++)
			memcpy(&cu->pLuma[CU_SIZE * y], &luma[CU_SIZE * y], sWidth * sizeof(pel));
	}
	else
	{
		for (s16 y = sHeight - 8; y >= 0; y -= 8)
		{
			for (s16 x = sWidth - 8; x >= 0; x -= 8)
			{
				ind--;
				const s32 *basisTable = &basis_weights[12 * cu->basisSelect[ind]];

				for (u8 i = 0; i < 8; i += 2)
				{
					for (u8 j = 0; j < 8; j += 2)
					{
						s32 p1 = cu->pLuma[CU_SIZE * (y + i) / 2 + (x + j) / 2];
						s32 p2 = cu->pLuma[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff];
						s32 p3 = cu->pLuma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2];
						s32 p4 = cu->pLuma[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff];

						luma[CU_SIZE * (y + i) + (x + j)] = (64 * p1 + basisTable[0] * p2 + basisTable[4] * p3 + basisTable[8] * p4 + 128) / 256;
						luma[CU_SIZE * (y + i) + (x + j + 1)] = (64 * p1 + basisTable[1] * p2 + basisTable[5] * p3 + basisTable[9] * p4 + 128) / 256;
						luma[CU_SIZE * (y + i + 1) + (x + j)] = (64 * p1 + basisTable[2] * p2 + basisTable[6] * p3 + basisTable[10] * p4 + 128) / 256;
						luma[CU_SIZE * (y + i + 1) + (x + j + 1)] = (64 * p1 + basisTable[3] * p2 + basisTable[7] * p3 + basisTable[11] * p4 + 128) / 256;

						p1 = cu->pChroma1[CU_SIZE * (y + i) / 2 + (x + j) / 2];
						p2 = cu->pChroma1[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff];
						p3 = cu->pChroma1[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2];
						p4 = cu->pChroma1[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff];

						chroma1[CU_SIZE * (y + i) + (x + j)] = (64 * p1 + basisTable[0] * p2 + basisTable[4] * p3 + basisTable[8] * p4 + 128) / 256;
						chroma1[CU_SIZE * (y + i) + (x + j + 1)] = (64 * p1 + basisTable[1] * p2 + basisTable[5] * p3 + basisTable[9] * p4 + 128) / 256;
						chroma1[CU_SIZE * (y + i + 1) + (x + j)] = (64 * p1 + basisTable[2] * p2 + basisTable[6] * p3 + basisTable[10] * p4 + 128) / 256;
						chroma1[CU_SIZE * (y + i + 1) + (x + j + 1)] = (64 * p1 + basisTable[3] * p2 + basisTable[7] * p3 + basisTable[11] * p4 + 128) / 256;

						p1 = cu->pChroma2[CU_SIZE * (y + i) / 2 + (x + j) / 2];
						p2 = cu->pChroma2[CU_SIZE * (y + i) / 2 + (x + j) / 2 + xoff];
						p3 = cu->pChroma2[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2];
						p4 = cu->pChroma2[CU_SIZE * ((y + i) / 2 + yoff) + (x + j) / 2 + xoff];

						chroma2[CU_SIZE * (y + i) + (x + j)] = (64 * p1 + basisTable[0] * p2 + basisTable[4] * p3 + basisTable[8] * p4 + 128) / 256;
						chroma2[CU_SIZE * (y + i) + (x + j + 1)] = (64 * p1 + basisTable[1] * p2 + basisTable[5] * p3 + basisTable[9] * p4 + 128) / 256;
						chroma2[CU_SIZE * (y + i + 1) + (x + j)] = (64 * p1 + basisTable[2] * p2 + basisTable[6] * p3 + basisTable[10] * p4 + 128) / 256;
						chroma2[CU_SIZE * (y + i + 1) + (x + j + 1)] = (64 * p1 + basisTable[3] * p2 + basisTable[7] * p3 + basisTable[11] * p4 + 128) / 256;
					}
				}
			}
		}

		for (u8 y = 0; y < sHeight; y++)
		{
			memcpy(&cu->pLuma[CU_SIZE * y], &luma[CU_SIZE * y], sWidth * sizeof(pel));
			memcpy(&cu->pChroma1[CU_SIZE * y], &chroma1[CU_SIZE * y], sWidth * sizeof(pel));
			memcpy(&cu->pChroma2[CU_SIZE * y], &chroma2[CU_SIZE * y], sWidth * sizeof(pel));
		}
	}

	// Adjust the basis index
	cu->basisIndex = ind;
}
