#include "CImg.h"
#include "ComDef.h"
#include "PICEncCU.h"
#include "PICComBACCore.h"

#include <cassert>

using namespace cimg_library;

int main(int argc, char *argv[])
{
	// For now, always require four arguments - input file, output file, chroma subsampling mode, and the quantizer
	if (argc != 5)
	{
		printf("Usage:\n");
		printf("  SiCEnc <input file> <output file> <440 | 420 | 400> <quant>\n");
		printf("  Note: Only bitmaps are currently supported.\n");
		return -1;
	}

	// Load a test image
	CImg<s16> img(argv[1]);

	// Determine the original dimensions
	u16 trueWidth = img.width();
	u16 trueHeight = img.height();
	u8 channels = img.spectrum();

	// Calculate the modulo-64 dimensions
	u16 width = trueWidth;
	u16 height = trueHeight;
	if (trueWidth % CU_SIZE)
		width = trueWidth + (CU_SIZE - (trueWidth % CU_SIZE));
	if (trueHeight % CU_SIZE)
		height = trueHeight + (CU_SIZE - (trueHeight % CU_SIZE));

	// Pad the image
	img.resize(width, height, 1, channels , 0, 1);

	// Define the parameters
	paramStruct param;
	param.qp = (u16)strtol(argv[4], NULL, 10);
	param.chromaSub = strcmp(argv[3], "400") == 0 ? CHROMA_400 : (strcmp(argv[3], "420") ? CHROMA_420 : CHROMA_444);

	// Automatically detect grayscale images
	if (channels == 1)
		param.chromaSub = CHROMA_400;

	// Apply the Rec. 601 YCbCr transform
	const u32 chOffset = width * height;
	if (channels == 3)
	{
		for (u32 i = 0; i < chOffset; i++)
		{
			double y = round(0.299 * (double)img(i) + 0.587 * (double)img(i + chOffset) + 0.114 * (double)img(i + 2 * chOffset));
			double cb = round(128 - 0.168736 * (double)img(i) - 0.331264 * (double)img(i + chOffset) + 0.5 * (double)img(i + 2 * chOffset));
			double cr = round(128 + 0.5 * (double)img(i) - 0.418688 * (double)img(i + chOffset) - 0.081312 * (double)img(i + 2 * chOffset));

			img(i) = (s16)((y < 0) ? 0 : ((y > 255) ? 255 : y));
			img(i + chOffset) = (s16)((cb < 0) ? 0 : ((cb > 255) ? 255 : cb));
			img(i + 2 * chOffset) = (s16)((cr < 0) ? 0 : ((cr > 255) ? 255 : cr));
		}
	}
	CImg<s16> imgRef = img.get_channel(0);

	// Initialize the arithmetic encoder with 8 bpp for grayscale images or 24 bpp for color images
	if (param.chromaSub == CHROMA_400)
		InitializeEncoder(8 * width * height);
	else
		InitializeEncoder(24 * width * height);

	// Process each coding unit (CU)
	for (u16 y = 0; y < height; y += CU_SIZE)
	{
		for (u16 x = 0; x < width; x += CU_SIZE)
		{
			PICEncCU(&img(x, y), width, height, param);
		}
	}

	if (param.chromaSub == CHROMA_400)
	{
		channels = 1;
		img.channel(0);
	}

	printf("%.4f\n", img.PSNR(imgRef, 255));
	//CImgDisplay disp(img, "", 0);
	img.save("tmp.bmp");
	img.display();

	// Terminate the encoding process and write to file
	EncodeTerminate();
	//WriteBitstream(argv[2], trueWidth, trueHeight, channels, param);

	return 0;
}