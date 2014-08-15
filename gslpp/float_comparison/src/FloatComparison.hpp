/*
 * FloatComparison.hpp
 *
 *  Created on: Jul 5, 2014
 *      Author: alinsch
 */

#include <stdint.h>
#include <cstring>
#include <cmath>

namespace gslpp{
namespace float_comparison{

///bitwise mask to get the sign bit, shift 1 to the first bit
const uint64_t signBitMaskDouble = (uint64_t)1 << (64-1);

///bitwise mask to get the mantissa bits, set all bits to one and shift exponent+sign bits back
const uint64_t mantissaBitMaskDouble = (~(uint64_t)0) >> 12;

///bitwise mask to get the exponent bits, flip all bits to 1 that are not 1 by sign or mantissa bit
const uint64_t exponentBitMaskDouble = ~(signBitMaskDouble|mantissaBitMaskDouble);

inline bool is_NaN_in_binary_double( uint64_t a_AsIntergerFromByteOrder );

template<>
bool equal_up_to_significant_digits_binary<double>(double a, double b, int significantDigits){
	//
	// #float steps required by the significant digits
	const uint64_t numberOfAllowedFloatStepsAAndBMayDiffer = (uint64_t)1 << significantDigits;
	//
	//convert a to its bitwise representation
	uint64_t a_AsIntergerFromByteOrder;
	std::memcpy(&a_AsIntergerFromByteOrder,&a,sizeof(double));
	uint64_t b_AsIntergerFromByteOrder;
	std::memcpy(&b_AsIntergerFromByteOrder,&b,sizeof(double));
	//
	//Check if a or b is NAN, NAN means the exponent is all 1 and some mantissa bits are not
	if ( ((a_AsIntergerFromByteOrder & exponentBitMaskDouble) == exponentBitMaskDouble) and
			((a_AsIntergerFromByteOrder & mantissaBitMaskDouble) != 0)){
		return false;
	}
	if ( ((b_AsIntergerFromByteOrder & exponentBitMaskDouble) == exponentBitMaskDouble) and
			((b_AsIntergerFromByteOrder & mantissaBitMaskDouble) != 0)){
		return false;
	}
	//
	//Check if a is infinity, return true if b is also infinity with the same sign
	if ( ((a_AsIntergerFromByteOrder & exponentBitMaskDouble) == exponentBitMaskDouble) and
			((a_AsIntergerFromByteOrder & mantissaBitMaskDouble) == 0)){
		//a is infinity
		if ( (a_AsIntergerFromByteOrder & signBitMaskDouble) == (b_AsIntergerFromByteOrder & signBitMaskDouble)){
			if ( ((b_AsIntergerFromByteOrder & exponentBitMaskDouble) == exponentBitMaskDouble) &&
					((b_AsIntergerFromByteOrder & mantissaBitMaskDouble) == 0)){
				//b is infinity with the same sign
				return true;
			}
		} else {
			return false;
		}
	}
	//
	//for the same sign representable floats AND their integer representation are adjacent to each other
	//	go to a representation where this is also the case across zero
	//	http://en.wikipedia.org/wiki/Two%27s_complement
	if ( (a_AsIntergerFromByteOrder & signBitMaskDouble) != 0 )//a negative
		a_AsIntergerFromByteOrder = (~(uint64_t)0)/2 - a_AsIntergerFromByteOrder;
	if ( (b_AsIntergerFromByteOrder & signBitMaskDouble) != 0 )//b negative
		b_AsIntergerFromByteOrder = (~(uint64_t)0)/2 - b_AsIntergerFromByteOrder;
	return (a_AsIntergerFromByteOrder <= b_AsIntergerFromByteOrder) ?
			( numberOfAllowedFloatStepsAAndBMayDiffer  > (b_AsIntergerFromByteOrder-a_AsIntergerFromByteOrder) )
				:
			( numberOfAllowedFloatStepsAAndBMayDiffer  > (a_AsIntergerFromByteOrder-b_AsIntergerFromByteOrder) );
}

template<>
bool equal_up_to_significant_digits_decimal<double>(double a, double b, int significantDigits){
	const double log2=0.3010299957;
	const int significantDigitsBase2 = std::ceil(significantDigits/log2);
	//
	//Determine the number of digits in binary of the mantissa that is allowed to vary
	return equal_up_to_significant_digits_binary(a,b,52+1-significantDigitsBase2);//1 stored implicitly
}


template<typename T>
bool equal_up_to_significant_digits_decimal_above_threshold(
		T a, T b, int significantDigits, T threshold){
	//
	//impose threshold condition
	if ( (std::fabs(b) < threshold) and (std::fabs(a) < threshold) )
		return true;
	//
	//threshold passed, check if a agrees with b on the first digits
	return equal_up_to_significant_digits_decimal(a,b,significantDigits);
}


template<>
bool is_NaN<double>( double a) {
	uint64_t a_AsIntergerFromByteOrder;
	std::memcpy(&a_AsIntergerFromByteOrder,&a,sizeof(double));
	return is_NaN_in_binary_double(a_AsIntergerFromByteOrder);
}

inline bool is_NaN_in_binary_double( uint64_t a_AsIntergerFromByteOrder ) {
	//
	//Check if a or b is NAN, NAN means the exponent is all 1 and some mantissa bits are not
	if ( ((a_AsIntergerFromByteOrder & exponentBitMaskDouble) == exponentBitMaskDouble) and
			((a_AsIntergerFromByteOrder & mantissaBitMaskDouble) != 0)){
		return true;
	}
	return false;
}

} /*namespace float_comparison*/
} /*namespace gslpp*/
