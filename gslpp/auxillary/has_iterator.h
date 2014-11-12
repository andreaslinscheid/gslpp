/*
 * has_iterator.h
 *
 *  Created on: Nov 5, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_HAS_ITERATOR_H_
#define GSLPP_AUXILLARY_HAS_ITERATOR_H_

#include <type_traits>

namespace gslpp{
namespace auxillary{

template <typename T>
class has_iterator{
private:
    template <typename U=T>
    static constexpr std::true_type check(typename U::iterator* x);

    template <typename U=T>
    static constexpr std::false_type check(...);

    typedef decltype(check<>(0)) type;
public:
    static const bool value = type::value;
};

}
}

#endif /* GSLPP_AUXILLARY_HAS_ITERATOR_H_ */
