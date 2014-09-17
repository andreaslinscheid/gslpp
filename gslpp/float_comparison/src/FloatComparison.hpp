/*
 * FloatComparison.hpp
 *
 *  Created on: Jul 5, 2014
 *      Author: alinsch
 */

#include <stdint.h>
#include <cstring>
#include <cmath>
#include <complex>

namespace gslpp{
namespace float_comparison{

//	here follow implementation details
/// @cond Doxygen_suppress

//trait to cast a float type into its corresponding bit pattern container
template <typename T> struct bitPatternType {};
template<> struct bitPatternType<double> {
	typedef uint64_t type;
};
template<> struct bitPatternType<float> {
	typedef uint32_t type;
};

//bitwise mask trait to get the sign bit, shift 1 to the first bit
template <typename T> struct signBitMask {}; //no value defined in the default case which must not compile!
template<> struct signBitMask<double> {
	static const uint64_t value = (uint64_t)1 << (64-1);
};
template<> struct signBitMask<float> {
	static const uint32_t value = (uint32_t)1 << (32-1);
};

//bitwise mask to get the mantissa bits, set all bits to one and shift exponent+sign bits back
template <typename T> struct mantissaBitMask {}; //no value defined in the default case which must not compile!
template<> struct mantissaBitMask<double> {
	static const uint64_t value = (~(uint64_t)0) >> 12;
};
template<> struct mantissaBitMask<float> {
	static const uint32_t value = (~(uint32_t)0) >> 9;
};

//bitwise mask to get the exponent bits, flip all bits to 1 that are not 1 by sign or mantissa bit
template <typename T> struct  exponentBitMask {
	static const typename bitPatternType<T>::type value = ~(signBitMask<T>::value | mantissaBitMask<T>::value );
};

template <typename T> struct numMantissaBits {}; //no value defined in the default case which must not compile!
template <> struct numMantissaBits<double> {
	static const uint64_t value = 52;
};
template <> struct numMantissaBits<float> {
	static const uint64_t value = 23;
};

//Check the bit pattern if it signals a NaN
//NaN means the exponent is all 1 and some mantissa bits are not
template<typename T>
bool bit_pattern_check_is_NaN( typename bitPatternType<T>::type a_AsIntergerFromByteOrder ){
	if ( ((a_AsIntergerFromByteOrder & exponentBitMask<T>::value) == exponentBitMask<T>::value) and
			((a_AsIntergerFromByteOrder & mantissaBitMask<T>::value) != 0)){
		return true;
	}
	return false;
};

//We emulate function template partial specialization that is not allowed by the C++ standard
//by creating an auxiliary class with a member functions that does the job.
template<typename T>
struct is_NaN_impl {
	static bool is_NaN_call( T a ) {
		typename bitPatternType<T>::type a_AsIntergerFromByteOrder;
		std::memcpy(&a_AsIntergerFromByteOrder,&a,sizeof( typename bitPatternType<T>::type));
		return bit_pattern_check_is_NaN<T>(a_AsIntergerFromByteOrder);
	};
};
//Partial template specialization for complex in terms of the type of complex
template<typename T>
struct is_NaN_impl< std::complex<T> > {
	static bool is_NaN_call( std::complex<T> a ) {
		return is_NaN<T>( a.real() ) or is_NaN<T>( a.imag() );
	};
};

//Implementation of the main template function that deduces to the member of the class
//which can be partially specialized.
template<typename T>
bool is_NaN(T a) {
	return is_NaN_impl<T>::is_NaN_call(a);
};

template<typename T>
struct equal_up_to_significant_digits_binary_impl {
	static bool equal_up_to_significant_digits_binary_call(T a, T b, int significantDigits) {

		//#float steps required by the significant digits
		const typename bitPatternType<T>::type
			numberOfAllowedFloatStepsAAndBMayDiffer =
					static_cast<typename bitPatternType<T>::type>(1) << significantDigits;

		//convert a and b to their bitwise representation
		typename bitPatternType<T>::type a_AsIntergerFromByteOrder;
		std::memcpy(&a_AsIntergerFromByteOrder,&a,sizeof(typename bitPatternType<T>::type));
		typename bitPatternType<T>::type b_AsIntergerFromByteOrder;
		std::memcpy(&b_AsIntergerFromByteOrder,&b,sizeof(typename bitPatternType<T>::type));

		//Check if the bit pattern of a or b signals NAN
		if ( bit_pattern_check_is_NaN<T>(a_AsIntergerFromByteOrder) )
			return false;
		if ( bit_pattern_check_is_NaN<T>(b_AsIntergerFromByteOrder) )
			return false;

		//for the same sign representable floats AND their integer representation are adjacent to each other
		//	go to a representation where this is also the case across zero
		//	http://en.wikipedia.org/wiki/Two%27s_complement
		if ( (a_AsIntergerFromByteOrder & signBitMask<T>::value) != 0 )//a negative
			a_AsIntergerFromByteOrder = (~static_cast<typename bitPatternType<T>::type>(0))/2 - a_AsIntergerFromByteOrder;
		if ( (b_AsIntergerFromByteOrder & signBitMask<T>::value) != 0 )//b negative
			b_AsIntergerFromByteOrder = (~static_cast<typename bitPatternType<T>::type>(0))/2 - b_AsIntergerFromByteOrder;
		return (a_AsIntergerFromByteOrder <= b_AsIntergerFromByteOrder) ?
				( numberOfAllowedFloatStepsAAndBMayDiffer  > (b_AsIntergerFromByteOrder-a_AsIntergerFromByteOrder) )
					:
				( numberOfAllowedFloatStepsAAndBMayDiffer  > (a_AsIntergerFromByteOrder-b_AsIntergerFromByteOrder) );
	};
};

template<typename T>
struct equal_up_to_significant_digits_binary_impl<std::complex<T> >{
	static bool equal_up_to_significant_digits_binary_call(
			std::complex<T> a, std::complex<T> b, int significantDigits){
		return equal_up_to_significant_digits_binary<T>( a.real() , b.real(), significantDigits ) and
					equal_up_to_significant_digits_binary<T>( a.imag() , b.imag(), significantDigits );
	}
};

template<typename T>
bool equal_up_to_significant_digits_binary(T a, T b, int significantDigits) {
	return equal_up_to_significant_digits_binary_impl<T>::equal_up_to_significant_digits_binary_call(a,b,significantDigits);
};

template<typename T>
struct equal_up_to_significant_digits_decimal_impl{
	static bool equal_up_to_significant_digits_decimal_call(T a, T b, int significantDigits) {
		const int significantDigitsBase2 = std::ceil(significantDigits/0.3010299957);

		//Determine the number of digits in binary of the mantissa that is allowed to vary
		return equal_up_to_significant_digits_binary(
				a,b,numMantissaBits<T>::value +1-significantDigitsBase2);//1 stored implicitly
	};
};

//template specialization for a complex type
template<typename T>
struct equal_up_to_significant_digits_decimal_impl< std::complex<T> >{
	static bool equal_up_to_significant_digits_decimal_call( std::complex<T> a,
			std::complex<T> b, int significantDigits){
		return equal_up_to_significant_digits_decimal<T>(a.real(),b.real(),significantDigits) and
				equal_up_to_significant_digits_decimal<T>(a.imag(),b.imag(),significantDigits);
	};
};

template<typename T>
bool equal_up_to_significant_digits_decimal(T a, T b, int significantDigits){
	return equal_up_to_significant_digits_decimal_impl<T>::equal_up_to_significant_digits_decimal_call(a,b,significantDigits);
}

template<typename T>
struct equal_up_to_significant_digits_decimal_above_threshold_impl{
	static bool equal_up_to_significant_digits_decimal_above_threshold_call(T a, T b, int significantDigits, T threshold){

		//impose threshold condition
		if ( (std::fabs( b ) < threshold) and (std::fabs(a ) < threshold) )
			return true;

		//threshold passed, check if a agrees with b on the first digits
		return equal_up_to_significant_digits_decimal(a,b,significantDigits);
	}
};

template<typename T>
struct equal_up_to_significant_digits_decimal_above_threshold_impl<std::complex<T> > {
	static bool equal_up_to_significant_digits_decimal_above_threshold_call(
			std::complex<T> a, std::complex<T> b, int significantDigits,std::complex<T> threshold){
		return equal_up_to_significant_digits_decimal_above_threshold_call(
				a.real(),b.real(),significantDigits,threshold.real()) and
				equal_up_to_significant_digits_decimal_above_threshold_call(
				a.imag(),b.imag(),significantDigits,threshold.imag());
	};
};

template<typename T>
bool equal_up_to_significant_digits_decimal_above_threshold(T a, T b, int significantDigits, T threshold){
	return equal_up_to_significant_digits_decimal_above_threshold_impl<T>::
			equal_up_to_significant_digits_decimal_above_threshold_call(a,b,significantDigits,threshold);
};

/// @end Doxygen_suppress
} /*namespace float_comparison*/
} /*namespace gslpp*/
