// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ////////////////////////////////////////////////////////////////////
// VectorStdOpsDecl.hpp

#ifndef TSFCORE_VECTOR_STD_OPS_DECL_HPP
#define TSFCORE_VECTOR_STD_OPS_DECL_HPP

#include "TSFCoreTypes.hpp"

namespace TSFCore {

/** \defgroup TSFCore_VectorStdOps_grp Collection of standard vector operations.
 */
//@{

///
/** Sum of vector elements: <tt>result = sum( v(i), i = 1...v.space()->dim() )</tt>.
 */
template<class Scalar>
Scalar sum( const Vector<Scalar>& v );

///
/** Natural norm: <tt>result = sqrt(<v,v>)</tt>.
 *
 * Returns <tt>Teuchos::ScalarTraits<Scalar>::squareroot(v.space()->scalarProd(v,v))</tt>
 */
template<class Scalar>
Scalar norm( const Vector<Scalar>& v );

///
/** One (1) norm: <tt>result = ||v||1</tt>
 */
template<class Scalar>
Scalar norm_1( const Vector<Scalar>& v );

///
/** Euclidean (2) norm: <tt>result = ||v||2</tt>
 */
template<class Scalar>
Scalar norm_2( const Vector<Scalar>& v );

///
/** Weighted Euclidean (2) norm: <tt>result = sqrt( sum( w(i)*v(i)^2) )</tt>
 */
template<class Scalar>
Scalar norm_2( const Vector<Scalar> &w, const Vector<Scalar>& v );

///
/** Infinity norm: <tt>result = ||v||inf</tt>
 */
template<class Scalar>
Scalar norm_inf( const Vector<Scalar>& v_rhs );

///
/** Dot product: <tt>result = x'*y</tt>
 */
template<class Scalar>
Scalar dot( const Vector<Scalar>& x, const Vector<Scalar>& y );

///
/** Get single element: <tt>result = v(i)</tt>
 */
template<class Scalar>
Scalar get_ele( const Vector<Scalar>& v, Index i );

///
/** Set single element: <tt>v(i) = alpha</tt>
 */
template<class Scalar>
void set_ele( Index i, Scalar alpha, Vector<Scalar>* v );

///
/** Assign all elements to a scalar: <tt>y(i) = alpha, i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void assign( Vector<Scalar>* y, const Scalar& alpha );

///
/** Vector assignment: <tt>y(i) = x(i), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void assign( Vector<Scalar>* y, const Vector<Scalar>& x );

///
/** Add a scalar to all elements: <tt>y(i) += alpha, i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void Vp_S( Vector<Scalar>* y, const Scalar& alpha );

///
/** Scale all elements by a scalar: <tt>y(i) *= alpha, i = 1...y->space()->dim()</tt>
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>y = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 */
template<class Scalar>
void Vt_S( Vector<Scalar>* y, const Scalar& alpha );

///
/** AXPY update: <tt>y(i) = alpha * x(i) + y(i), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void Vp_StV( Vector<Scalar>* y, const Scalar& alpha, const Vector<Scalar>& x );

///
/** <tt>y(i) = abs(x(i)), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void abs( Vector<Scalar>* y, const Vector<Scalar>& x );

///
/** <tt>y(i) = 1/x(i), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void reciprocal( Vector<Scalar>* y, const Vector<Scalar>& x );

///
/** <tt>y(i) += alpha * x(i) * v(i), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void ele_wise_prod( const Scalar& alpha, const Vector<Scalar>& x, const Vector<Scalar>& v, Vector<Scalar>* y );

///
/** <tt>y(i) = alpha * x(i) / v(i), i = 1...y->space()->dim()</tt>
 */
template<class Scalar>
void ele_wise_divide( const Scalar& alpha, const Vector<Scalar>& x, const Vector<Scalar>& v, Vector<Scalar>* y );

///
/** <tt>y(i) = beta*y(i) + sum( alpha[k]*x[k](i), k=0...m-1 ), i = 1...y->space()->dim()</tt>.
 *
 * @param  m          [in] Number of vectors x[]
 * @param  alpha      [in] Array (length <tt>m</tt>) of input scalars.
 * @param  x          [in] Array (length <tt>m</tt>) of input vectors.
 * @param  beta       [in] Scalar multiplier for y
 * @param  y          [in/out] Target vector that is the result of the linear combination.
 *
 * This function implements a general linear combination:
 \verbatim
 y(i) = beta*y(i) + alpha[0]*x[0](i) + alpha[1]*x[1](i) + ... + alpha[m-1]*x[m-1](i), i = 1...y->space()->dim()

 \endverbatim
 */
template<class Scalar>
void linear_combination(
	const int                m
	,const Scalar            alpha[]
	,const Vector<Scalar>*   x[]
	,const Scalar            &beta
	,Vector<Scalar>          *y
	);

///
/** Seed the random number generator used in <tt>random_vector</tt>
 */
template<class Scalar>
void seed_randomize( unsigned int );

///
/** Generate a random vector with elements uniformly distrubuted
 * elements.
 * 
 * The elements <tt>v->getEle(i)</tt> are randomly generated between
 * <tt>[l,u]</tt>.
 *
 * The seed is set using <tt>seed_randomize()</tt>
 */
template<class Scalar>
void randomize( Scalar l, Scalar u, Vector<Scalar>* v );

//@}

} // end namespace TSFCore

// ////////////////////////////
// Inline functions

template<class Scalar>
inline
Scalar TSFCore::norm( const Vector<Scalar>& v )
{
	return Teuchos::ScalarTraits<Scalar>::squareroot(v.space()->scalarProd(v,v));
}


#endif // TSFCORE_VECTOR_STD_OPS_DECL_HPP
