/*
 * LinearAlgebra.h
 *
 *  Created on: Aug 24, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_LINEAR_ALGEBRA_LINEARALGEBRA_H_
#define GSLPP_LINEAR_ALGEBRA_LINEARALGEBRA_H_

#include <vector>
#include <complex>

namespace gslpp {
namespace linear_algebra {

/** A collection of linear algebra functions.
 *
 *	We do not implement the full BLAS / LAPACK framework. This collection is intended for internal uses to limit code
 *	dependencies.
 */
template<typename T>
struct LinearAlgebra {

	/**	A function to invert a square matrix.
	 *
	 * @param squareMatrix Container with the square matrix in row major order (on input) and the inverse (on output)
	 */
	static void invert_square_matrix(std::vector<T> & squareMatrix);

	/**	A function that performs an LU factorization
	 *
	 * @param matrixOrder (input) memory layout of the %matrix
	 * @param m The number of rows of the matrix
	 * @param n The number of columns of the matrix
	 * @param a (input/output) array, dimension (m,n)
	 *          On entry, the dimRow-by-dimCol matrix to be factored.
	 *          On exit, the factors L and U from the factorization
	 *          matrix = P*L*U; the unit diagonal elements of L are not stored.
	 * @param lda  The leading dimension of array a.
	 * @param ipiv (output) integer array, dimension (min(m,n))
	 *          The pivot indices
	 * @return  = 0:  successful exit
	 *          < -i: the i-th argument had an illegal value
	 *          > i:  U(i,i) is exactly zero. The factorization
	 *                has been completed, but the factor U is exactly
	 *                singular, and division by zero will occur if it is used
	 *                to solve a system of equations.
	 */
	static int getrf( int matrixOrder, int m, int n,  std::vector<T> &a, int lda, std::vector<int> &ipiv );

	/**	A function compute the inverse of an LU factorized matrix
	 *
	 * @param matrixOrder (input) memory layout of the %matrix
	 * @param n The number of columns of the matrix
	 * @param a (input/output) array, dimension (dimRow,dimCol)
	 *          On entry, the dimRow-by-dimCol matrix to be factored.
	 *          On exit, the factors L and U from the factorization
	 *          matrix = P*L*U; the unit diagonal elements of L are not stored.
	 * @param lda  The leading dimension of array a.
	 * @param ipiv (output) integer array, dimension (min(m,n))
	 *          The pivot indices
	 *          < -i: the i-th argument had an illegal value
	 *          > i:  U(i,i) is exactly zero. The factorization
	 *                has been completed, but the factor U is exactly
	 *                singular, and division by zero will occur if it is used
	 *                to solve a system of equations.
	 */
	static int getri(int matrixOrder, int n, std::vector<T> &a, int lda, std::vector<int> const & ipiv);

	static void inversion_of_matrices_info_checks(int info);
};

} /* namespace linear_algebra */
} /* namespace gslpp */

#include "gslpp/linear_algebra/src/LinearAlgebra.hpp"
#endif /* GSLPP_LINEAR_ALGEBRA_LINEARALGEBRA_H_ */
