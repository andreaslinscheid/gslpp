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

	bool first_lower_than_second_impl(T first, T second) const {
		return first < second;
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
		static_cast<NumAccuracyControl<T>*>(this)->set_local_error_threshold(localRelativeErrorThreshold,localAbsErrorThreshold);
		static_cast<NumAccuracyControl<T>*>(this)->set_global_error_threshold(globalRelativeErrorThreshold,globalAbsErrorThreshold);
	}

	bool first_lower_than_second_impl(std::complex<T> first, std::complex<T> second) const {
		return (first.real() < second.real())
				and (first.imag() < second.imag());
	}

	bool not_all_zero(std::complex<T> const& val) const{
		return (val.real() != 0 ) and (val.imag() != 0 );
	}

	std::complex<T> relative_value(std::complex<T> const& functionValue, std::complex<T> const& rel) const {
		return functionValue*rel;
	}
};

template<typename T>
class NumAccuracyControl_impl<NumAccuracyControl<T>, T,true> {
public:
	NumAccuracyControl_impl() {
	}

	bool first_lower_than_second_impl(T first, T second) const {
		T diff = first + second *
				static_cast< typename std::iterator_traits<typename T::iterator>::value_type >(-1.0);
		return std::any_of(diff.begin(), diff.end(), [] (
				typename std::iterator_traits<typename T::iterator>::value_type x) {return x > 0; });
	}


	bool not_all_zero(T const& val) const{
		return std::any_of(val.begin(), val.end(), [] (
				typename std::iterator_traits<typename T::iterator>::value_type x) {return x != 0; });
	}

	T relative_value(T const& functionValue, T const& rel) const {
		T result = rel;
		typename T::iterator itr = result.begin();
		std::for_each(functionValue.begin(),functionValue.end(),
				[&](typename std::iterator_traits<typename T::iterator>::value_type x)
				{(*itr) *= x; ++itr;});
		return result;
	}
private:
};

template<typename T>
NumAccuracyControl<T>::NumAccuracyControl() : NumAccuracyControl_impl<NumAccuracyControl<T>,T,has_iterator<T>::value >(){
	 _maxNumberOfSubdivisions =1000;
}

template<typename T>
bool NumAccuracyControl<T>::locally_sufficient(T const& localErrEstimate, T const& functionValue) const {
	bool relConv = true;
	if ( _checkRelLocal )
		relConv = this->first_lower_than_second_impl(localErrEstimate,
					this->relative_value(functionValue,_localRelativeErrorThreshold));
	bool absConv = true;
	if ( _checkAbsLocal )
		absConv = this->first_lower_than_second_impl( localErrEstimate, _localAbsErrorThreshold);
	return relConv and absConv;
}

template<typename T>
bool NumAccuracyControl<T>::global_sufficient(T const& globalErrEstimate, T const& functionValue) const {
	bool relConv = true;
	if ( _checkRelGlobal )
		relConv = this->first_lower_than_second_impl(globalErrEstimate,
					this->relative_value(functionValue,_globalRelativeErrorThreshold) );
	bool absConv = true;
	if ( _checkAbsGlobal )
		absConv = this->first_lower_than_second_impl( globalErrEstimate, _globalAbsErrorThreshold);
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
