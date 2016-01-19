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

template<typename T, bool T_is_class>
struct has_iterator_impl;

/**
 * A Class that checks if the typename T knows a type iterator.
 *
 * Thus, gslpp::auxillary::has_iterator<T>::value is true if T::iterator exists
 * or false if it does not.
 */
template <typename T>
struct has_iterator : public has_iterator_impl<T,std::is_class<T>::value>
{
};

template <typename T>
struct has_iterator_impl<T,false>
{
    static const bool value = false;
};

template <typename T>
class has_iterator_impl<T,true>{
private:

	struct Fallback { int iterator;};
	struct Derived : public T, public Fallback {};

	template <typename U, U> struct ChT;

	//The following will work only, if T does not define the iterator type, otherwise
	//	this will result in an ambiguity and the second version of check will be chosen
    template <typename U>
    static constexpr std::false_type check( ChT<int Fallback::*,&U::iterator> *);

    template <typename U>
    static constexpr std::true_type check(...);

    typedef decltype(check<Derived>(0)) type;
public:
    static const bool value = type::value;
};

}
}

#endif /* GSLPP_AUXILLARY_HAS_ITERATOR_H_ */
