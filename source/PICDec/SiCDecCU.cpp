#include "SiCDecCU.h"

void SiCDecCU(s16 *img, u16 width, u16 height, paramStruct param)
{
	// Create a buffer for the coding blocks (CB) and the sub-image coding blocks (SCB)
	const u16 cbNumPix = CU_SIZE * CU_SIZE;

	static s16 cbLuma[cbNumPix], cbResLuma[cbNumPix];
	static s16 cbChroma1[cbNumPix], cbResChroma1[cbNumPix];
	static s16 cbChroma2[cbNumPix], cbResChroma2[cbNumPix];

	static s16 vscbLuma[cbNumPix / 8];
	static s16 vscbChroma1[cbNumPix / 8];
	static s16 vscbChroma2[cbNumPix / 8];

	static s16 hscbLuma[cbNumPix / 8];
	static s16 hscbChroma1[cbNumPix / 8];
	static s16 hscbChroma2[cbNumPix / 8];

	// Define the horizontal and vertical sub-coding units
	static cuStruct cu = { cbLuma, cbChroma1, cbChroma2, CU_SIZE, CU_SIZE, CU_SIZE, param.chromaSub };
	static scuStruct hscu = { hscbLuma, hscbChroma1, hscbChroma2, CU_SIZE / 8, CU_SIZE, CU_SIZE / 8, param.chromaSub };
	static scuStruct vscu = { vscbLuma, vscbChroma1, vscbChroma2, CU_SIZE, CU_SIZE / 8, CU_SIZE, param.chromaSub };

	// Fill all components with zeros
	memset(cbLuma, 0, cbNumPix * sizeof(s16));
	memset(cbChroma1, 0, cbNumPix * sizeof(s16));
	memset(cbChroma2, 0, cbNumPix * sizeof(s16));

	memset(vscbLuma, 0, cbNumPix / 8 * sizeof(s16));
	memset(vscbChroma1, 0, cbNumPix / 8 * sizeof(s16));
	memset(vscbChroma2, 0, cbNumPix / 8 * sizeof(s16));

	memset(hscbLuma, 0, cbNumPix / 8 * sizeof(s16));
	memset(hscbChroma1, 0, cbNumPix / 8 * sizeof(s16));
	memset(hscbChroma2, 0, cbNumPix / 8 * sizeof(s16));

	// Decode the sub-image transform units (STUs)
	static stuStruct stu = { &vscu, &hscu, CU_SIZE, CU_SIZE, CU_SIZE, param.chromaSub };
	stuDec(stu);
	dequantTable(stu, param.q1);

	// Apply the inverse DCT to the sub-images
	tuInverse(stu);

	// Copy back the reconstructed vertical SCU pixels
	copyFromSCU(cu, hscu, vscu);

	// Define the prediction unit (PU)
	static u8 pbModeLuma[cbNumPix / 64];
	static u8 pbModeChroma[cbNumPix / 64];
	static ScanDir scanLuma[cbNumPix / 64];
	static ScanDir scanChroma[cbNumPix / 64];
	static puStruct pu = { &cu, pbModeLuma, pbModeChroma, scanLuma, scanChroma, CU_SIZE / 8, CU_SIZE / 8, CU_SIZE / 8, param.chromaSub };

	// Decode the prediction modes and apply prediction
	puDec(pu);
	pred(pu);

	// Decode and transform the residual
	static rtuStruct rtu = { cbResLuma, cbResChroma1, cbResChroma2, scanLuma, scanChroma, CU_SIZE, CU_SIZE, CU_SIZE, param.chromaSub };
	rtuDec(rtu);
	dequantConst(rtu, param.q2);
	tuInverse(rtu);

	// Reconstruct the CU
	static cuStruct cuRes = { cbResLuma, cbResChroma1, cbResChroma2, CU_SIZE, CU_SIZE, CU_SIZE, param.chromaSub };
	cu = cu + cuRes;

	// Copy pixels to the image
	copyFromCU(img, width, height, cu);
}

void copyFromCU(s16 *img, u16 width, u16 height, cuStruct &cu)
{
	// Calculate the channel offset
	const u32 chOffset = width * height;

	// Copy luma pixels to the CBs and the SCBs
	for (u8 y = 0; y < cu.height; y++)
	{
		for (u8 x = 0; x < cu.width; x++)
		{
			img[width * y + x] = cu.cbLuma[cu.stride * y + x];
		}
	}

	// Copy chroma pixels to the CBs and the SCBs
	if (cu.chromaSub == CHROMA_444)
	{
		for (u8 y = 0; y < cu.height; y++)
		{
			for (u8 x = 0; x < cu.width; x++)
			{
				img[width * y + x + chOffset] = cu.cbChroma1[cu.stride * y + x];
				img[width * y + x + 2 * chOffset] = cu.cbChroma2[cu.stride * y + x];
			}
		}
	}
	else if (cu.chromaSub == CHROMA_420)
	{
		for (u8 y = 0; y < cu.height; y += 2)
		{
			for (u8 x = 0; x < cu.width; x += 2)
			{
				img[width * y + x + chOffset] = cu.cbChroma1[cu.stride * (y / 2) + x / 2];
				img[width * y + x + 1 + chOffset] = cu.cbChroma1[cu.stride * (y / 2) + x / 2];
				img[width * (y + 1) + x + chOffset] = cu.cbChroma1[cu.stride * (y / 2) + x / 2];
				img[width * (y + 1) + x + 1 + chOffset] = cu.cbChroma1[cu.stride * (y / 2) + x / 2];


				img[width * y + x + 2 * chOffset] = cu.cbChroma2[cu.stride * (y / 2) + x / 2];
				img[width * y + x + 1 + 2 * chOffset] = cu.cbChroma2[cu.stride * (y / 2) + x / 2];
				img[width * (y + 1) + x + 2 * chOffset] = cu.cbChroma2[cu.stride * (y / 2) + x / 2];
				img[width * (y + 1) + x + 1 + 2 * chOffset] = cu.cbChroma2[cu.stride * (y / 2) + x / 2];
			}
		}
	}
}

void copyFromSCU(cuStruct &cu, scuStruct &hscu, scuStruct &vscu)
{
	// Copy pixels from the CU to the SCU
	for (u8 y = 7; y < cu.height; y += 8)
	{
		for (u8 x = 0; x < cu.width; x++)
		{
			cu.cbLuma[cu.stride * x + y] = hscu.scbLuma[hscu.stride * x + y / 8];
			cu.cbChroma1[cu.stride * x + y] = hscu.scbChroma1[hscu.stride * x + y / 8];
			cu.cbChroma2[cu.stride * x + y] = hscu.scbChroma2[hscu.stride * x + y / 8];

			if (x % 8 < 7)
			{
				cu.cbLuma[cu.stride * y + x] = vscu.scbLuma[vscu.stride * (y / 8) + x];
				cu.cbChroma1[cu.stride * y + x] = vscu.scbChroma1[vscu.stride * (y / 8) + x];
				cu.cbChroma2[cu.stride * y + x] = vscu.scbChroma2[vscu.stride * (y / 8) + x];
			}
			else
			{
				cu.cbLuma[cu.stride * y + x] = (cu.cbLuma[cu.stride * y + x] + vscu.scbLuma[vscu.stride * (y / 8) + x] + 1) / 2;
				cu.cbChroma1[cu.stride * y + x] = (cu.cbChroma1[cu.stride * y + x] + vscu.scbChroma1[vscu.stride * (y / 8) + x] + 1) / 2;
				cu.cbChroma2[cu.stride * y + x] = (cu.cbChroma2[cu.stride * y + x] + vscu.scbChroma2[vscu.stride * (y / 8) + x] + 1) / 2;
			}
		}
	}
}
