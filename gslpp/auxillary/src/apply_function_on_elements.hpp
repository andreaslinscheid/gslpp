/*
 * apply_function_on_elements.hpp
 *
 *  Created on: Nov 6, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_APPLY_FUNCTION_ON_ELEMENTS_HPP_
#define GSLPP_AUXILLARY_APPLY_FUNCTION_ON_ELEMENTS_HPP_

#include "gslpp/auxillary/apply_function_on_elements.h"
#include <complex>

namespace gslpp{
namespace auxillary{

namespace delegate{
template<class method, bool resultHasMultipleElements>
struct apply_function_on_elements_impl{ };

template<class method>
struct apply_function_on_elements_impl<method,false>{

	typedef typename FunctionTraits<method>::result_type TRes;

	template <typename ...TArgs>
	static void call (TRes &result,TArgs... args){
		result = method(result,args...);
	};

//	//specialization for complex where the method should effect only the elements
//	template <typename ...TArgs>
//	static void call ( std::complex<TRes> &result,TArgs... args){
//		result.real( method( result.real() , args...) );
//		result.imag( method( result.imag() , args...) );
//	};
};

template<class method>
struct apply_function_on_elements_impl<method,true>{
	template <typename ...TArgs>
	static void call (typename FunctionTraits<method>::result_type &result,TArgs... args){
		typename FunctionTraits<method>::result_type::iterator it;
		for (it = result.begin(); it != result.end(); ++it){
			*it=method(*it, args...);
		}
	}
};

};/* namespace delegate */
};/* namespace auxillary */
};/* namespace gslpp */


#endif /* APPLY_FUNCTION_ON_ELEMENTS_HPP_ */
