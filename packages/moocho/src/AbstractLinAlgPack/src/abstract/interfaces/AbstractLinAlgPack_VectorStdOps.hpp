// ////////////////////////////////////////////////////////////////////
// VectorStdOps.h

#ifndef ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H
#define ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H

#include <utility>

#include "VectorWithOpMutable.h"

namespace AbstractLinAlgPack {

/** \defgroup VectorStdOps_grp Collection of useful vector operations.
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

///
/** Computes the maximum positive and negative step that can be taken
  * that are within the relaxed bounds.
  *
  *	This function returns and computes the maximum (in magnitude) postive
  *	(return.first) and negative (return.second) steps u that can be taken
  *	such that the relaxed bounds:
  \verbatim
  xl - max_bnd_viol <= x + u * d <= xu - max_bnd_viol
  \endverbatim
  * are strictly satisfied.
  *
  * If return.first < 0.0 then this is a flag that x is not
  * in the relaxed bounds to begin with.  In this case
  * return.second has no meaning.
  */
std::pair<value_type,value_type>
max_near_feas_step(
	const VectorWithOp& x, const VectorWithOp& d
	,const VectorWithOp& xl, const VectorWithOp& xu
	,value_type max_bnd_viol
	); 

///
/** Count the number of finitly bounded elements in <tt>xl <= x <= xu</tt>.
 *
 * ToDo: Finish documentation!
 */
size_type num_bounded(
	const VectorWithOp& xl, const VectorWithOp& xu
	,value_type inf_bound );

//@}

/** \defgroup VectorStdOps_TOp_grp Transformation operations */
//@{

///
/** v_lhs *= alpha
 */
void Vt_S( VectorWithOpMutable* v_lhs, const value_type& alpha );

///
/** v_lhs = alpha * v_rhs + v_lhs
 */
void Vp_StV( VectorWithOpMutable* v_lhs, const value_type& alpha, const VectorWithOp& v_rhs );

///
/** v_lhs(i) = alpha * v_rhs1(i) * v_rhs2(i), i = 1,,,dim.
 */
void ele_wise_prod(
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

///
/** Force a vector in bounds.
 *
 * ToDo: Finish documentation!
 */
void force_in_bounds( const VectorWithOp& xl, const VectorWithOp& xu, VectorWithOpMutable* x );

//@}

//@}

} // end namespace AbstractLinAlgPack

// //////////////////////////////////////////////
// Inline implementations

#endif // ABSTRACT_LINALG_PACK_VECTOR_STD_OPS_H
