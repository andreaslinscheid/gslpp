/*
 * LinearAlgebra.hpp
 *
 *  Created on: Aug 24, 2014
 *      Author: alinsch
 */
#include "gslpp/linear_algebra/LinearAlgebra.h"
#include "gslpp/error_handling/Error.h"

#ifdef LIN_ALG_EXTERNAL

#ifdef USE_MKL

#define MKL_INT int
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include <mkl_types.h>
#include <mkl_lapacke.h>
#include <mkl.h>

#else

extern "C" {
#define lapack_complex_float std::complex<float>
#define lapack_complex_double std::complex<double>
#include <cblas.h>
#include <lapacke.h>
}

#endif

#endif /* LIN_ALG_EXTERNAL */

namespace gslpp {
namespace linear_algebra {

template<typename T>
void LinearAlgebra<T>::invert_square_matrix(std::vector<T> & squareMatrix) {
	if ( squareMatrix.empty() )
		return;
	int dim = sqrt(static_cast<int>(squareMatrix.size()));
#ifdef DEBUG_BUILD
	if ( static_cast<size_t>(dim * dim) != squareMatrix.size())
		gslpp::error_handling::Error( "input appears to be not a square matrix",
				gslpp::error_handling::Error::INPUT_ERROR);
#endif
	std::vector<int> ipiv;

	// L U factorization
	ipiv.clear();
	ipiv.resize( dim );
	int info = LinearAlgebra<T>::getrf(LAPACK_ROW_MAJOR,dim,dim,squareMatrix,dim,ipiv);
	gslpp::linear_algebra::LinearAlgebra<T>::inversion_of_matrices_info_checks(info);

	//compute the inverse of the L U factorized matrix
	info = LinearAlgebra<T>::getri(LAPACK_ROW_MAJOR,dim,squareMatrix,dim,ipiv);
	gslpp::linear_algebra::LinearAlgebra<T>::inversion_of_matrices_info_checks(info);
};

#ifdef LIN_ALG_EXTERNAL
template<>
inline int LinearAlgebra< std::complex<double> >::getrf( int matrixOrder, int m, int n,
		std::vector< std::complex<double> > &a, int lda, std::vector<int> &ipiv){
	return LAPACKE_zgetrf(matrixOrder,
			static_cast<int>(m),static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

template<>
inline int LinearAlgebra< std::complex<float> >::getrf( int matrixOrder, int m, int n,
		std::vector< std::complex<float> > &a, int lda, std::vector<int> &ipiv){
	return LAPACKE_cgetrf(matrixOrder,
			static_cast<int>(m),static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

template<>
inline int LinearAlgebra<float>::getrf( int matrixOrder, int m, int n,
		std::vector<float> &a, int lda, std::vector<int> &ipiv){
	return LAPACKE_sgetrf(matrixOrder,
			static_cast<int>(m),static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

template<>
inline int LinearAlgebra<double>::getrf(int matrixOrder, int m, int n,
		std::vector<double> &a, int lda, std::vector<int> &ipiv){
	return LAPACKE_dgetrf(matrixOrder,
			static_cast<int>(m),static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

template<>
inline int LinearAlgebra<double>::getri(int matrixOrder, int n,
		std::vector<double> &a, int lda, std::vector<int> const & ipiv){
	return LAPACKE_dgetri(matrixOrder,static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

template<>
inline int LinearAlgebra<std::complex<double> >::getri(int matrixOrder, int n,
		std::vector<std::complex<double> > &a, int lda, std::vector<int> const & ipiv){
	return LAPACKE_zgetri(matrixOrder,static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

template<>
inline int LinearAlgebra<float>::getri(int matrixOrder, int n,
		std::vector<float> &a, int lda, std::vector<int> const & ipiv){
	return LAPACKE_sgetri(matrixOrder,static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

template<>
inline int LinearAlgebra<std::complex<float> >::getri(int matrixOrder, int n,
		std::vector<std::complex<float> > &a, int lda, std::vector<int> const & ipiv){
	return LAPACKE_cgetri(matrixOrder,static_cast<int>(n),
			&a[0],static_cast<int>(n),
			&ipiv[0]);
};

#else
	///@todo Implement a generic version of getrf that does not use an external BLAS and LAPACK implementation
#endif

template<typename T>
void LinearAlgebra<T>::inversion_of_matrices_info_checks(int info){
	gslpp::error_handling::Warning warning;
	if ( info < 0 )
		warning << "Argument "<< info <<" had an illegal value";
	if ( info > 0)
		warning << "U("<< info <<","<< info <<") is exactly zero. "
				<<"The factorization has been completed, but the factor U is exactly singular.";
}

} /* namespace linear_algebra */
} /* namespace gslpp */
