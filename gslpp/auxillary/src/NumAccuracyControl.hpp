/*
 * NumAccuracyControl.hpp
 *
 *  Created on: Oct 20, 2014
 *      Author: alinsch
 */

#include "gslpp/auxillary/NumAccuracyControl.h"
#include "gslpp/auxillary/AccuracyGoal.h"
#include <algorithm>

namespace gslpp {
namespace auxillary {

template<class derived,typename T, bool THasIterator>
class NumAccuracyControl_impl { };

namespace detail{

//We delegade the comparision of lower and equality to this templates.
//	While complex does not have a strict ordering in the mathematical sense
//	these functions are intended for component-wise convergence comparision.
template<typename T>
struct cmp_lower{

	static bool cmp(T a, T b) {
		return (a < b);
	}
};

template<typename T>
struct cmp_lower<std::complex<T> >{

	static bool cmp(std::complex<T> a, std::complex<T> b) {
		return cmp_lower<T>::cmp(a.real(),b.real())
				and cmp_lower<T>::cmp(a.imag(),b.imag());
	}
};

template<typename T>
struct cmp_leq{

	static bool cmp(T a, T b) {
		return (a < b) or ( (not (a < b)) and (not (b < a)));
	}
};

template<typename T>
struct cmp_leq<std::complex<T> >{

	static bool cmp(std::complex<T> a, std::complex<T> b) {
		return cmp_leq<T>(a.real(),b.real())
				and cmp_leq<T>(a.imag(),b.imag());
	}
};

template<typename T>
struct cmp_leq_abs{

	static bool cmp(T a, T b) {
		return cmp_leq<T>::cmp(a < T(0) ? -a : a,b < T(0) ? -b : b);
	}
};

template<typename T>
struct cmp_leq_abs<std::complex<T> >{

	static bool cmp(std::complex<T> a, std::complex<T> b) {
		return cmp_leq_abs<T>::cmp(a.real(),b.real())
				and cmp_leq_abs<T>::cmp(a.imag(),b.imag());
	}
};

}; /* namespace details */

template<typename T>
class NumAccuracyControl_impl< NumAccuracyControl<T>, T,false> {
public:
	NumAccuracyControl_impl() {
		T localRelativeErrorThreshold = ( AccuracyGoal<T>::value );
		T localAbsErrorThreshold = (  AccuracyGoal<T>::value );
		T globalRelativeErrorThreshold =(   AccuracyGoal<T>::value );
		T globalAbsErrorThreshold = (  AccuracyGoal<T>::value );
		static_cast<NumAccuracyControl<T>*>(this)->set_local_error_threshold(localRelativeErrorThreshold,localAbsErrorThreshold);
		static_cast<NumAccuracyControl<T>*>(this)->set_global_error_threshold(globalRelativeErrorThreshold,globalAbsErrorThreshold);
	}

	bool abs_first_leq_than_second_impl(T first, T second) const {
		return detail::cmp_leq_abs<T>::cmp(first,second);
	}

	bool first_lower_than_second_impl(T first, T second) const {
		return detail::cmp_lower<T>::cmp(first,second);
	}

	bool first_leq_than_second_impl(T first, T second) const {
		return detail::cmp_leq<T>::cmp(first,second);
	}

	bool not_all_zero(T const& val) const{
		return val != 0;
	}

	T relative_value(T const& functionValue, T const& rel) const {
		return functionValue*rel;
	}
};

template<typename T>
class NumAccuracyControl_impl< NumAccuracyControl< std::complex<T> >, std::complex<T>,false> {
public:
	NumAccuracyControl_impl() {
		std::complex<T> localRelativeErrorThreshold = ( AccuracyGoal<T>::value );
		std::complex<T> localAbsErrorThreshold = (  AccuracyGoal<T>::value );
		std::complex<T> globalRelativeErrorThreshold =(   AccuracyGoal<T>::value );
		std::complex<T> globalAbsErrorThreshold = (  AccuracyGoal<T>::value );
		static_cast<NumAccuracyControl<std::complex<T> >*>(this)->set_local_error_threshold(localRelativeErrorThreshold,localAbsErrorThreshold);
		static_cast<NumAccuracyControl<std::complex<T> >*>(this)->set_global_error_threshold(globalRelativeErrorThreshold,globalAbsErrorThreshold);
	}

	bool abs_first_leq_than_second_impl(std::complex<T> first, std::complex<T> second) const {
		return detail::cmp_leq_abs<std::complex<T> >::cmp(first,second);
	}

	bool first_lower_than_second_impl(std::complex<T> first, std::complex<T> second) const {
		return detail::cmp_lower<std::complex<T> >::cmp(first,second);
	}

	bool first_leq_than_second_impl(std::complex<T> first, std::complex<T> second) const {
		return detail::cmp_leq<std::complex<T> >::cmp(first,second);
	}

	bool not_all_zero(std::complex<T> const& val) const{
		return (val != std::complex<T>(0));
	}

	std::complex<T> relative_value(std::complex<T> const& functionValue, std::complex<T> const& rel) const {
		return functionValue*rel;
	}
};

template<typename T>
class NumAccuracyControl_impl<NumAccuracyControl<T>, T,true> {
public:

	typedef typename std::iterator_traits<typename T::iterator>::value_type value_type;

	bool abs_first_leq_than_second_impl(T first, T second) const {
		T diff = first + second * value_type(-1.0);
		return std::any_of(diff.begin(), diff.end(),
				[] (value_type x) {return detail::cmp_leq_abs<value_type>::cmp(x,value_type(0)); });
	}

	bool first_leq_than_second_impl(T first, T second) const {
		T diff = first + second * value_type(-1.0);
		return std::any_of(diff.begin(), diff.end(),
				[] (value_type x) {return detail::cmp_leq<value_type>::cmp(x,value_type(0)); });
	}

	bool first_lower_than_second_impl(T first, T second) const {
		T diff = first + second * value_type(-1.0);
		return std::any_of(diff.begin(), diff.end(),
				[] (value_type x) {return detail::cmp_lower<value_type>::cmp(x,value_type(0)); });
	}

	bool not_all_zero(T const& val) const{
		return std::any_of(val.begin(), val.end(), [] (value_type x) {return x != value_type(0); });
	}

	T relative_value(T const& functionValue, T const& rel) const {
		T result = rel;
		auto itr = result.begin();
		std::for_each(functionValue.begin(),functionValue.end(),
				[&](value_type x){(*itr) *= x; ++itr;});
		return result;
	}
};

template<typename T>
NumAccuracyControl<T>::NumAccuracyControl() : NumAccuracyControl_impl<NumAccuracyControl<T>,T,has_iterator<T>::value >(){
	 _maxNumberOfSubdivisions =1000;
}

template<typename T>
bool NumAccuracyControl<T>::locally_sufficient(T const& localErrEstimate, T const& functionValue) const {
	bool relConv = true;
	//equality is also OK
	if ( _checkRelLocal )
		relConv = (functionValue == T(0)) or this->abs_first_leq_than_second_impl(localErrEstimate,
					this->relative_value(functionValue,_localRelativeErrorThreshold));
	bool absConv = true;
	if ( _checkAbsLocal )
		absConv = this->abs_first_leq_than_second_impl( localErrEstimate, _localAbsErrorThreshold);
	return relConv and absConv;
}

template<typename T>
bool NumAccuracyControl<T>::global_sufficient(T const& globalErrEstimate, T const& functionValue) const {
	bool relConv = true;
	if ( _checkRelGlobal )
		relConv = (functionValue == T(0)) or this->abs_first_leq_than_second_impl(globalErrEstimate,
					this->relative_value(functionValue,_globalRelativeErrorThreshold) );
	bool absConv = true;
	if ( _checkAbsGlobal )
		absConv = this->abs_first_leq_than_second_impl( globalErrEstimate, _globalAbsErrorThreshold);
	return relConv and absConv;
}

template<typename T>
bool NumAccuracyControl<T>::sub_divisions_below_max(size_t numOfSubdivisions) const {
	return _maxNumberOfSubdivisions >= numOfSubdivisions;
}

template<typename T>
bool NumAccuracyControl<T>::first_lower_than_second(T const& first, T const& second) const {
	return this->first_lower_than_second_impl(first,second);
}

template<typename T>
void NumAccuracyControl<T>::set_local_error_threshold(T const& thrRel, T const& thrAbs){
	_localRelativeErrorThreshold = thrRel;
	_checkRelLocal = this->not_all_zero(_localRelativeErrorThreshold);
	_localAbsErrorThreshold = thrAbs;
	_checkAbsLocal = this->not_all_zero(_localAbsErrorThreshold);
}

template<typename T>
void NumAccuracyControl<T>::set_max_subdiv(size_t n){
	_maxNumberOfSubdivisions = n;
}

template<typename T>
void NumAccuracyControl<T>::set_global_error_threshold(T const& thrRel, T const& thrAbs){
	_globalRelativeErrorThreshold = thrRel;
	_checkRelGlobal = this->not_all_zero(_globalRelativeErrorThreshold);
	_globalAbsErrorThreshold = thrAbs;
	_checkAbsGlobal = this->not_all_zero(_globalAbsErrorThreshold);
}

template<typename T>
size_t NumAccuracyControl<T>::get_max_num_subdiv() const {
	return _subdivisions;
}

template<typename T>
T NumAccuracyControl<T>::get_abs_error_estimate() const{
	return _errorEstimateAbs;
}

template<typename T>
T NumAccuracyControl<T>::get_rel_error_estimate() const{
	return _errorEstimateRel;
}

} /* namespace auxillary */
} /* namespace gslpp */
