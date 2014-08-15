/** @file gslpp/float_comparison/FloatComparison.h
 * Contains declarations of float comparison function.
 *
 *  Created on: Jul 5, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_FLOAT_COMPARISON_H_
#define GSLPP_FLOAT_COMPARISON_H_

namespace gslpp{
namespace float_comparison{

/** Compare two float number if they do not differ more than significantDigits binary digits.
 *	Comparison done according to
 *	http://randomascii.wordpress.com/2012/01/23/stupid-float-tricks-2/
 * @param a a floating point number (T = float, double, complex<float>, complex<double>)
 * @param b a floating point number (T = float, double, complex<float>, complex<double>)
 * @param significantDigits integer with a allowed binary digits to differ
 * @return 	false if 1) a and/or b is NAN
 * 					2) a or b is infinity and the other one is not
 * 					3) a and b differ by more than significantDigits
 * 			true otherwise
 */
template<typename T>
bool equal_up_to_significant_digits_binary(T a, T b, int significantDigits);

/** Compare two float number if they do not differ more than significantDigits decimal digits.
 *
 * 	significantDigits is cast into the significant digits in base 2, rounding up.
 *  Thus, it may behave a little more strict than expected.
 *  significantDigits in binary is used in
 * 	gslpp::float_comparison::equal_up_to_significant_digits_binary(T a, T b, int significantDigits)
 */
template<typename T>
bool equal_up_to_significant_digits_decimal(T a, T b, int significantDigits);

/** Check two doubles a and b for equality with
 * 	gslpp::float_comparison::equal_up_to_significant_digits_decimal(T a, T b, int significantDigits)
 * 	 or if both are below threshold.
 */
template<typename T>
bool equal_up_to_significant_digits_decimal_above_threshold(
		T a, T b, int significantDigits, T threshold);

/**	Check a float number if it is NaN.
 *
 * @param a A float or double number.
 * @return true if the bit pattern matches NaN, false else.
 */
template<typename T>
bool is_NaN(T a);

}/*namespace float_comparison*/
}/*namespace gslpp*/

#include "gslpp/float_comparison/src/FloatComparison.hpp"
#endif /* GSLPP_FLOAT_COMPARISON_H_ */
