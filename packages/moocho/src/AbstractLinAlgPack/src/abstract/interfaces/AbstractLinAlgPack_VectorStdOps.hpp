// ////////////////////////////////////////////////////////////////////
// AbstractLinAlgPack_VectorStdOps.hpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#ifndef ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H
#define ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H

#include "AbstractLinAlgPack_VectorMutable.hpp"

namespace AbstractLinAlgPack {

/** \defgroup VectorStdOps_grp Collection of standard vector operations.
 */
//@{

/** \defgroup VectorStdOps_ROp_grp Reduction operations */
//@{

///
/** result = sum( v_rhs(i), i = 1,,,dim )
 */
value_type sum( const Vector& v_rhs );

///
/** result = v_rhs1' * v_rhs2
 */
value_type dot( const Vector& v_rhs1, const Vector& v_rhs2 );

///
/** result = v_rhs1' * sv_rhs2
 */
value_type dot( const Vector& v_rhs1, const SpVectorSlice& sv_rhs2 );

///
/** result = sv_rhs1' * v_rhs2
 */
value_type dot( const SpVectorSlice& sv_rhs1, const Vector& v_rhs2 );

///
/** Compute the maximum element in a vector.
 *
 * @param  v        [in] The vector being searched
 * @param  max_v_j  [out] The value of the element with the max abs value.
 * @param  max_j    [out] The index of the element with the max abs value.
 *
 * Returns:
 \verbatim

 max_v_j = v(max_j) s.t. |v(max_j)| <= |v(j)|, for j = 1...n
 \endverbatim
 * If there is a tie, the lowest index is returned so that the
 * result is unique no matter what order the vector elements are
 * searched.
 */
void max_abs_ele( const Vector& v, value_type* max_v_j, index_type* max_j ); 

//@}

/** \defgroup VectorStdOps_TOp_grp Transformation operations */
//@{

///
/** v_lhs += alpha
 */
void Vp_S( VectorMutable* v_lhs, const value_type& alpha );

///
/** v_lhs *= alpha
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>v_lhs = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 */
void Vt_S( VectorMutable* v_lhs, const value_type& alpha );

///
/** v_lhs = alpha * v_rhs + v_lhs
 */
void Vp_StV( VectorMutable* v_lhs, const value_type& alpha, const Vector& v_rhs );

///
/** v_lhs = alpha * sv_rhs + v_lhs
 */
void Vp_StV( VectorMutable* v_lhs, const value_type& alpha, const SpVectorSlice& sv_rhs );

///
/** v_lhs(i) += alpha * v_rhs1(i) * v_rhs2(i), i = 1,,,dim.
 */
void ele_wise_prod(
	const value_type& alpha, const Vector& v_rhs1, const Vector& v_rhs2
	,VectorMutable* v_lhs );

///
/** v_lhs(i) = alpha * v_rhs1(i) / v_rhs2(i), i = 1,,,dim.
 */
void ele_wise_divide(
	const value_type& alpha, const Vector& v_rhs1, const Vector& v_rhs2
	,VectorMutable* v_lhs );

///
/** Seed the random number generator
  */
void seed_random_vector_generator( unsigned int );

///
/** Generate a random vector with elements uniformly
  * distrubuted elements.
  * 
  * The elements are randomly generated between <tt>[l,u]</tt>.
  */
void random_vector( value_type l, value_type u, VectorMutable* v );

///
/** Compute the sign of each element in an input vector.
 *
 \verbatim

         / -1.0 : if v(i)  < 0.0
 z(i) =  |  0.0 : if v(i) == 0.0
         \ +1.0 : if v(i)  < 0.0

 , for i = 1...n
 \endverbatim
 */
void sign(
	const Vector      &v
	,VectorMutable    *z
	);

//@}

//@}

} // end namespace AbstractLinAlgPack

// /////////////////////////////////////
// Inline members

inline
AbstractLinAlgPack::value_type
AbstractLinAlgPack::dot( const SpVectorSlice& sv_rhs1, const Vector& v_rhs2 )
{
	return dot(v_rhs2,sv_rhs1);
}

#endif // ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H
