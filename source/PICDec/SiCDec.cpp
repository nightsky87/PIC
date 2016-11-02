#include "CImg.h"
#include "ComDef.h"
#include "SiCDecCU.h"
#include "SiCComBACCore.h"
#include <cmath>

using namespace cimg_library;

int main(int argc, char *argv[])
{
	// For now, always require two arguments - input file and output file
	if (argc != 3)
	{
		printf("Usage:\n");
		printf("  SiCDec <input file> <output file>\n");
		printf("  Note: Only bitmaps are currently supported.\n");
		return -1;
	}

	u8 channels;
	u16 trueWidth, trueHeight;

	paramStruct param;

	// Read the bitstream to memory and initialize the decoding process
	ReadBitstream(argv[1], &trueWidth, &trueHeight, &channels, &param);
	InitializeDecoder();

	// Calculate the modulo-64 dimensions
	u16 width = trueWidth;
	u16 height = trueHeight;
	if (trueWidth % CU_SIZE)
		width = trueWidth + (CU_SIZE - (trueWidth % CU_SIZE));
	if (trueHeight % CU_SIZE)
		height = trueHeight + (CU_SIZE - (trueHeight % CU_SIZE));

	// Allocate space for the YCbCr values
	CImg<s16> img(width, height, 1, channels, 0);

	// Process each coding unit (CU)
	for (u32 y = 0; y < height; y += CU_SIZE)
	{
		for (u32 x = 0; x < width; x += CU_SIZE)
		{
			SiCDecCU(&img(x, y), width, height, param);
		}
	}

	// Terminate the encoding process and write to file
	DecodeTerminate();

	if (param.chromaSub != CHROMA_400)
	{
		// Apply the Rec. 601 inverse color transform
		const u32 chOffset = width * height;
		for (u32 i = 0; i < chOffset; i++)
		{
			double r = round((double)img(i) + 1.402 * (double)(img(i + 2 * chOffset) - 128));
			double g = round((double)img(i) - 0.344136 * (double)(img(i + chOffset) - 128) - 0.714136 * (double)(img(i + 2 * chOffset) - 128));
			double b = round((double)img(i) + 1.772 * (double)(img(i + chOffset) - 128));

			img(i) = (s16)((r < 0) ? 0 : ((r > 255) ? 255 : r));
			img(i + chOffset) = (s16)((g < 0) ? 0 : ((g > 255) ? 255 : g));
			img(i + 2 * chOffset) = (s16)((b < 0) ? 0 : ((b > 255) ? 255 : b));
		}
	}

	// Save the image to the output file
	img.save(argv[2]);

	return 0;
}
