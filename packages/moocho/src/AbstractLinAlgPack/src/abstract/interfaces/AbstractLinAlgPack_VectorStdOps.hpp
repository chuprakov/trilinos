// ////////////////////////////////////////////////////////////////////
// VectorStdOps.h
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

#include "VectorWithOpMutable.h"

namespace AbstractLinAlgPack {

/** \defgroup VectorStdOps_grp Collection of standard vector operations.
 */
//@{

/** \defgroup VectorStdOps_ROp_grp Reduction operations */
//@{

///
/** result = sum( v_rhs(i), i = 1,,,dim )
 */
value_type sum( const VectorWithOp& v_rhs );

///
/** result = v_rhs1' * v_rhs2
 */
value_type dot( const VectorWithOp& v_rhs1, const VectorWithOp& v_rhs2 );

//@}

/** \defgroup VectorStdOps_TOp_grp Transformation operations */
//@{

///
/** v_lhs += alpha
 */
void Vp_S( VectorWithOpMutable* v_lhs, const value_type& alpha );

///
/** v_lhs *= alpha
 *
 * This takes care of the special cases of <tt>alpha == 0.0</tt>
 * (set <tt>v_lhs = 0.0</tt>) and <tt>alpha == 1.0</tt> (don't
 * do anything).
 */
void Vt_S( VectorWithOpMutable* v_lhs, const value_type& alpha );

///
/** v_lhs = alpha * v_rhs + v_lhs
 */
void Vp_StV( VectorWithOpMutable* v_lhs, const value_type& alpha, const VectorWithOp& v_rhs );

///
/** v_lhs = alpha * sv_rhs + v_lhs
 */
void Vp_StV( VectorWithOpMutable* v_lhs, const value_type& alpha, const SpVectorSlice& sv_rhs );

///
/** v_lhs(i) += alpha * v_rhs1(i) * v_rhs2(i), i = 1,,,dim.
 */
void ele_wise_prod(
	const value_type& alpha, const VectorWithOp& v_rhs1, const VectorWithOp& v_rhs2
	,VectorWithOpMutable* v_lhs );

///
/** v_lhs(i) = alpha * v_rhs1(i) / v_rhs2(i), i = 1,,,dim.
 */
void ele_wise_divide(
	const value_type& alpha, const VectorWithOp& v_rhs1, const VectorWithOp& v_rhs2
	,VectorWithOpMutable* v_lhs );

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
void random_vector( value_type l, value_type u, VectorWithOpMutable* v );

//@}

//@}

} // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H
