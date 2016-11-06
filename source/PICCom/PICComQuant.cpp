#include "PICComQuant.h"

s32 quantVal(const s32 val, const u8 qp)
{
	static s32 mul[6] = { 26214, 23302, 20560, 18396, 16384, 14564 };

	const u8 shift = qp / 6 + 19;
	const u8 ind = qp % 6;

	// Calculate the offset for rounding
	const s32 offset = (1 << shift) / 2;

	// Quantize and return the value
	return ((val * mul[ind] + offset) >> shift);
}

s32 dequantVal(const s32 val, const u8 qp)
{
	static s32 mul[6] = { 40, 45, 51, 57, 64, 72 };

	const u8 shift = qp / 6;
	const u8 ind = qp % 6;

	return (((mul[ind] << shift) * val) >> 1);
}