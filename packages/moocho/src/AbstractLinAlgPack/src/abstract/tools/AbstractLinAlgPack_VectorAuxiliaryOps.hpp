// ////////////////////////////////////////////////////////////////////
// VectorAuxiliaryOps.h
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

#ifndef ABSTRACT_LINALG_PACK_VECTOR_AUXILIARY_OPS_H
#define ABSTRACT_LINALG_PACK_VECTOR_AUXILIARY_OPS_H

#include <utility>

#include "VectorWithOpMutable.h"

namespace AbstractLinAlgPack {

/** \defgroup VectorAuxiliaryOps_grp Collection of auxiliary useful vector operations.
 */
//@{

/** \defgroup VectorAuxiliaryOps_ROp_grp Reduction operations */
//@{

///
/** Compute the maximum element in a vector.
 *
 * Returns:
  \verbatim

  max{ v(i), i = 1...n }
  \endverbatim
  */
value_type max( const VectorWithOp& v ); 

///
/** Computes the maximum positive and negative step that can be taken
  * that are within the relaxed bounds.
  *
  *	This function computes and returns the maximum (in magnitude) postive
  *	(<tt>return.first >= 0.0</tt>) and negative (<tt>return.second <= 0.0</tt>) steps
  * \c u that can be taken such that the relaxed bounds:
  \verbatim
  xl - max_bnd_viol <= x + u * d <= xu - max_bnd_viol
  \endverbatim
  * are strictly satisfied.
  *
  * If <tt>return.first < 0.0</tt> then this is a flag that \c x is not
  * in the relaxed bounds to begin with.  In this case \c return.second
  * has no meaning.
  */
std::pair<value_type,value_type>
max_near_feas_step(
	const VectorWithOp& x, const VectorWithOp& d
	,const VectorWithOp& xl, const VectorWithOp& xu
	,value_type max_bnd_viol
	); 

///
/** Computes the maximum relative step of <tt>x = x + d</tt>.
  *
  \verbatim

  return = max{ |d|/(1.0+|x(i)|), for i = 1...n }
  \endverbatim
  */
value_type max_rel_step(
	const VectorWithOp& x, const VectorWithOp& d
	);

///
/** Count the number of finitly bounded elements in <tt>xl <= x <= xu</tt>.
 *
 * ToDo: Finish documentation!
 */
size_type num_bounded(
	const VectorWithOp& xl, const VectorWithOp& xu
	,value_type inf_bound );

///
/** Computes the log barrier term:
 *
 * ToDo: Add documentation!
 */
value_type log_bound_barrier(
	const VectorWithOp    &x
	,const VectorWithOp   &xl
	,const VectorWithOp   &xu
	); 

///
/** Compute the maximum violation from a set of inequality constraints <tt>vL <= v <= vU</tt>.
 *
 * @param  v      [in] Inequality value vector.
 * @param  vL     [in] Lower inequality bounds (may be -infinity (i.e. very large negative number))
 * @param  vU     [in] Upper inequality bounds (may be +infinity (i.e. very large positive number))
 * @param  max_viol_i
 *                [out] Gives the index of the inequality with the maximum (scaled) violation.
 *                If <tt>*max_viol_i == 0</tt> on output then no inequality was violated.
 * @param  max_viol
 *                [out] The maximum (scaled violation).
 *                Only significant if <tt>*max_viol_i > 0</tt>.  
 * @param  v_i
 *                [out] Set to <tt>v.get_ele(*max_viol_i)</tt>.
 *                Only significant if <tt>*max_viol_i > 0</tt>.
 * @param  bnd_type
 *                [out] The type of bound with the maximum violation.
 *                <ul>
 *                    <li> -1 : LOWER
 *                    <li>  0 : EQUALITY
 *                    <li> +1 : UPPER
 *                </ul>
 *                Only significant if <tt>*max_viol_i > 0</tt>.
 * @param  vLU_i
 *                [out] Set to:
 *                <ul>
 *                    <li><tt>vL.get_ele(*max_viol_i)</tt> if <tt>*bnd_type <= 0</tt>
 *                    <li><tt>vU.get_ele(*max_viol_i)</tt> if <tt>*bnd_type > 0</tt>
 *                </ul>
 *                Only significant if <tt>*max_viol_i > 0</tt>.  
 *
 * @return Returns <tt>true</tt> if some constraint was violated.
 *
 * Preconditions:<ul>
 * <li> ToDo: Spell these out!
 * </ul>
 *
 * Postconditions:<ul>
 * <li> ToDo: Spell these out!
 * </ul>
 *
 * In order to make the result unique if more than one inequality
 * <tt>vL(i) <= v(i) <= vL(i)</tt> have the same maximum violation
 * then the inequality with the lowest <tt>i</tt> is returned.
 *
 */
bool max_inequ_viol(
	const AbstractLinAlgPack::VectorWithOp   &v
	,const AbstractLinAlgPack::VectorWithOp  &vL
	,const AbstractLinAlgPack::VectorWithOp  &vU
	,AbstractLinAlgPack::size_type           *max_viol_i
	,AbstractLinAlgPack::value_type          *max_viol
	,AbstractLinAlgPack::value_type          *v_i
	,int                                     *bnd_type
	,AbstractLinAlgPack::value_type          *vLU_i
	); 

//@}

/** \defgroup VectorAuxiliaryOps_TOp_grp Transformation operations */
//@{

///
/** Force a vector in bounds.
 *
 \verbatim

          / xl(i)  : if x(i) < xl(i)
  x(i) =  | x(i)   : if xl(i) <= x(i) <= xu(i)
          \ xu(i)  : if x(i) > xu(i)

  , for 1 = 1...n
 \endverbatim
 */
void force_in_bounds( const VectorWithOp& xl, const VectorWithOp& xu, VectorWithOpMutable* x );

///
/** Take the maximum value of the vector elements and a scalar.
 *
 \verbatim

 y(i) = max( y(i), min_ele ), for i = 1...n
 \endverbatim
 */
void max_vec_scalar(
	value_type              min_ele
	,VectorWithOpMutable    *y
	);

///
/** Take the maximum value of the absolute vector elements and a scalar.
 *
 \verbatim

 y(i) = max( fabs(y(i)), min_ele ), for i = 1...n
 \endverbatim
 */
void max_abs_vec_scalar(
	value_type              min_ele
	,VectorWithOpMutable    *y
	);

//@}

//@}

} // end namespace AbstractLinAlgPack

#endif // ABSTRACT_LINALG_PACK_VECTOR_AUXILIARY_OPS_H
