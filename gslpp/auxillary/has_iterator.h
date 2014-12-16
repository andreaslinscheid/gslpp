/*
 * has_iterator.h
 *
 *  Created on: Nov 5, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_HAS_ITERATOR_H_
#define GSLPP_AUXILLARY_HAS_ITERATOR_H_
/** \file has_iterator.h
    \brief Checking if an object or its parents define a type named iterator.
*/

#include <type_traits>

namespace gslpp{
namespace auxillary{

/**
 * A Class that checks if the typename T knows a type iterator.
 *
 * Thus, gslpp::auxillary::has_iterator<T>::value is true if T::iterator exists
 * or false if it does not.
 */
template <typename T, bool isClass = std::is_class<T>::value>
class has_iterator{
public:
    static const bool value = false;
};

template <typename T>
class has_iterator<T,true>{
private:

	struct Fallback { typedef int iterator;};
	struct Derived : T, Fallback {};

	template<typename C, C>
	struct ChT;

	//The following will work only, if T does not define the iterator type, otherwise
	//	this will result in an ambiguity and the second version of check will be chosen
    template <typename U>
    static constexpr std::false_type check( decltype(U::iterator) *);

    template <typename U>
    static constexpr std::true_type check(...);

    typedef decltype(check<Derived>(0)) type;
public:
    static const bool value = type::value;
};

}
}

#endif /* GSLPP_AUXILLARY_HAS_ITERATOR_H_ */
