#include "PICComQuant.h"
#include <cmath>

void quantConst(s32 *tb, u16 qp)
{
	// Calculate the scaling factor
	double qStep = pow(2, ((double)qp - 4) / 6);

	// Quantize all coefficients with the same factor
	for (u8 y = 0; y < 4; y++)
	{
		for (u8 x = 0; x < 4; x++)
		{
			tb[4 * y + x] = (s32)round((double)tb[4 * y + x] / qStep);
		}
	}
}

void quantConst(cpStruct *cp, u16 qp)
{
	// Calculate the scaling factor
	double qStep = pow(2, ((double)qp - 4) / 6);

	// Quantize all coefficients with the same factor
	if (cp->width != cp->height)
	{
		for (u8 y = 0; y < cp->height; y++)
		{
			for (u8 x = 0; x < cp->width; x += 2)
			{
				cp->pLuma[cp->width * y + x + 1] = (s32)round((double)cp->pLuma[cp->width * y + x + 1] / qStep);
			}
		}
	}
	else
	{
		for (u8 y = 0; y < cp->height; y += 2)
		{
			for (u8 x = 0; x < cp->width; x++)
			{
				cp->pLuma[cp->width * (y + 1) + x] = (s32)round((double)cp->pLuma[cp->width * (y + 1) + x] / qStep);
			}
		}
	}
}

void dequantConst(s32 *tb, u16 qp)
{
	// Calculate the scaling factor
	double qStep = pow(2, ((double)qp - 4) / 6);

	// Quantize all coefficients with the same factor
	for (u8 y = 0; y < 4; y++)
	{
		for (u8 x = 0; x < 4; x++)
		{
			tb[4 * y + x] = coeffCast((double)tb[4 * y + x] * qStep);
		}
	}
}

void dequantConst(cpStruct *cp, u16 qp)
{
	// Calculate the scaling factor
	double qStep = pow(2, ((double)qp - 4) / 6);

	// Dequantize all coefficients with the same factor
	if (cp->width != cp->height)
	{
		for (u8 y = 0; y < cp->height; y++)
		{
			for (u8 x = 0; x < cp->width; x += 2)
			{
				cp->pLuma[cp->width * y + x + 1] = coeffCast((double)cp->pLuma[cp->width * y + x + 1] * qStep);
			}
		}
	}
	else
	{
		for (u8 y = 0; y < cp->height; y += 2)
		{
			for (u8 x = 0; x < cp->width; x++)
			{
				cp->pLuma[cp->width * (y + 1) + x] = coeffCast((double)cp->pLuma[cp->width * (y + 1) + x] * qStep);
			}
		}
	}
}

s16 coeffCast(double val)
{
	val = round(val);
	return (s32)((val < -32768) ? -32768 : ((val > 32767) ? 32767 : val));
}
