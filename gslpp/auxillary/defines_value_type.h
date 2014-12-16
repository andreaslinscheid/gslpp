/*
 * defines_value_type.h
 *
 *  Created on: Nov 13, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_DEFINES_VALUE_TYPE_H_
#define GSLPP_AUXILLARY_DEFINES_VALUE_TYPE_H_

#include <type_traits>

namespace gslpp{
namespace auxillary{
/** \file defines_value_type.h
    \brief Checking if an object or its parents define a type named value_type.
*/

/**
 * A Class that checks if the typename T knows a type value_type.
 *
 * Thus, gslpp::auxillary::has_iterator<T>::value is true if T::value_type exists
 * or false if it does not.
 */
template <typename T, bool isClass = std::is_class<T>::value>
class defines_value_type{
public:
    static const bool value = false;
};

template <typename T>
class defines_value_type<T,true>{
private:

	struct Fallback { typedef void value_type;};
	struct Derived : T, Fallback {};

	template<typename C, C>
	struct ChT;

	//The following will work only, if T does not define the iterator type, otherwise
	//	this will result in an ambiguity and the second version of check will be chosen
    template <typename U>
    static constexpr std::false_type check( decltype(U::value_type) *);

    template <typename U>
    static constexpr std::true_type check(...);

    typedef decltype(check<Derived>(0)) type;
public:
    static const bool value = type::value;
};

}/* namespace auxillary */
}/* namespace gslpp */

#endif /* GSLPP_AUXILLARY_DEFINES_VALUE_TYPE_H_ */
