// ///////////////////////////////////////////////////////////
// NLPFirstOrderDirect.h
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

#ifndef NLP_FIRST_ORDER_DIRECT_H
#define NLP_FIRST_ORDER_DIRECT_H

#include "NLPObjGradient.h"

namespace NLPInterfacePack {

///
/** Interface providing only direct first order sensitivity information.
 *
 * <b>Overview:</b>
 *
 * This interface defines a basis for the equality constriants and then only
 * certain linear systems with this basis are solved for.  This interface is
 * useful in reduced space SQP-type and other related optimization algorithms.
 *
 * Specifically, the variables are partitioned into dependent and independent
 * sets <tt>x = [ x_dep'  x_indep' ]'</tt> and Jacobians of the constraints <tt>c(x)</tt>
 * and <tt>h(x)</tt> at the point <tt>x</tt> are:
 \verbatim

	del(c,x) = Gc' = [ del(c(con_decomp))   ] = [ GcD' ] = [ GcDD'  GcDI' ] = [ C  N ]
                     [ del(c(con_undecomp)) ]   [ GcU' ]   [ GcUD'  GcUI' ]   [ E  F ]

	del(h,x) = Gh' = [ GhD'  GhI' ]

    where:
      C <: R^(r x r) is nonsingular
      N <: R^(r x (n-r))
      E <: R^((m-r) x r)
      F <: R^((m-r) x (n-r))
 \endverbatim
 * This partitions the general equality constraints c(x) into two sets;
 * decomposed c(con_decomp) and undecomposed c(con_undecomp).  It is therefore
 * expected that sub-vectors and subspaces from <tt>space_x().sub_space(var_dep)</tt>,
 * <tt>space_x().sub_space(var_indep)</tt>, <tt>space_c().sub_space(con_decomp)</tt> and
 * <tt>space_c().sub_space(con_undecomp)</tt> can all be accessed.  Other sub-vectors
 * and sub-spaces may not be available (but the algorithm should not need access to other
 * sub-spaces).
 *
 * Free access to solves with the basis <tt>C</tt> is not given however and instead this interface
 * computes, for the current point \a x, the direct sensitivity matrices <tt>D = -inv(C)*N</tt>,
 * <tt>V = F - E * D</tt>, <tt>P = GhI' + GhD'* D</tt>, the auxiliary matrices
 * <tt>GcU = [ GcUD; GcUI ] = [ E';  F' ]</tt> and <tt>Gh = [ GhD;  GhI ]</tt>, and the Newton
 * step <tt>py = -inv(C)*c(con_decomp)</tt>.
 * In general, linear solves with the transpose with <tt>C</tt> are not possible and
 * therefore are not avalible.  A number of very specialized applications can only
 * provide this information but this is all that is needed by many numerical
 * optimization (and related) algorithms.
 *
 * <b>Client Usage:</b>
 *
 * The dimension of the basis matrix \c C is returned by \c r().  The ranges for the dependent and
 * independent varaibles are returned by \c var_dep() and \c var_indep().  The ranges for the
 * decomposed and undecomposed equality constraints are \c con_decomp() and \c con_undecomp().
 * Note that \c con_undecomp() may return an invalid range if there are no undecomposed equalities.
 *
 * Note that the matrix objects returned from \c space_GcU(), \c space_D(), \c space_V(),
 * \c space_P(), \c space_GcUD() and \c space_GhD() can not be expected to be usable until they are
 * passed to the calculation routines or have been intialized in some other way.
 *
 * <b>Subclass Developer's Notes:</b>
 *
 * The default implementation of this interface assumes that there are no undecomposed
 * equality constraints (i.e. <tt>this->con_decomp().size() == this->m()) and no general
 * inequality constraints (i.e. <tt>this->mI() == 0).
 *
 * ToDo: Finish Documentation!
 */
class NLPFirstOrderDirect : virtual public NLPObjGradient
{
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractLinAlgPack::MatrixSpace<MatrixWithOp> >          mat_space_ptr_t;

	/** @name Dimensionality */
	//@{

	///
	/** Returns the number of decomposed equality constraints (r <= m).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>this->con_decomp().size()</tt>.
	 * This implementation will work for all implementations.
	 */
	virtual size_type r() const;

	//@}

	/** @name Ranges for dependent and independent variables and decomposed and undecomposed equalities
	 */
	//@{

	///
	/** Return the range of dependent (i.e.\ basic) variables.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>Range1D(1,this->m())</tt>.
	 */
	virtual Range1D var_dep() const;
	///
	/** Return the range of independent (i.e.\ nonbasic) variables.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>Range1D(this->m()+1,this->n())</tt>.
	 */
	virtual Range1D var_indep() const;
	///
	/** Return the range of decomposed equality constraints.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>Range1D(1,this->m())</tt>.
	 */
	virtual Range1D con_decomp() const;
	///
	/** Return the range of undecomposed equality constraints.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>Range1D::Invalid</tt>.
	 */
	virtual Range1D con_undecomp() const;

	//@}

	/** @name MatrixSpace objects */
	//@{
	
	///
	/** Return a matrix space object for creating <tt>GcU</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation is to return <tt>return.get() == NULL</tt>.
	 * This is the proper implementation when <tt>m() == r()</tt>.
	 * When <tt>m() > r()</tt> then the subclass must override this method to
	 * return a valid matrix space object.  Moreover, the returned
	 * matrix object from <tt>this->space_GcU()->create_member()->get_sub_view(rng,Range1D())</tt>
	 * must be non-null for <tt>rng == this->var_dep()</tt> or <tt>rng == this->var_indep()</tt>.
	 * This gives access to the matrices <tt>E'</tt> and <tt>F'</tt> as shown above.
	 */
	virtual const mat_space_ptr_t& space_GcU() const;
	///
	/** Return a matrix space object for creating <tt>Gh</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation is to return <tt>return.get() == NULL</tt>.
	 * This is the proper implementation when <tt>m() == r()</tt>.  Moreover, the returned
	 * matrix object from <tt>this->space_Gh()->create_member()->get_sub_view(rng,Range1D())</tt>
	 * must be non-null for <tt>rng == this->var_dep()</tt> or <tt>rng == this->var_indep()</tt>.
	 * This gives access to the matrices <tt>GhD</tt> and <tt>GhI</tt> as shown above.
	 */
	virtual const mat_space_ptr_t& space_Gh() const;
	///
	/** Return a matrix space object for <tt>D = -inv(C)*N</tt> {abstract}.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual const mat_space_ptr_t& space_D() const = 0;
	///
	/** Return a matrix space object for <tt>V = F + E * D</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation is to return <tt>return.get() == NULL</tt>.
	 * This is the correct implementation when <tt>m() == r()</tt>.  However,
	 * when <tt>m() > r()</tt> this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_space_ptr_t& space_V() const;
	///
	/** Return a matrix space object for <tt>P = GhI' + GhD'* D</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation is to return <tt>return.get() == NULL</tt>.
	 * This is the correct implementation when <tt>m() == r()</tt>.  However,
	 * when <tt>m() > r()</tt> this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_space_ptr_t& space_P() const;
	///
	/** Return a matrix space object for a mutable matrix compatible with <tt>GcU(var_dep)</tt>.
	 *
	 * This matrix space object is designed to create mutable matrix objects compatible
	 * with <tt>GcU(var_dep)</tt>.  For example, a matrix object <tt>U</tt> created by this matrix space
	 * can be used to compute <tt>U = GcU(var_dep)' - GcU(var_indep)'*D'</tt> (this is needed
	 * by a orthogonal range/null decomposition.
	 *
	 * The default implementation is to return <tt>return.get() == NULL</tt>.
	 * This is the correct implementation when <tt>m() == r()</tt>.  However,
	 * when <tt>m() > r()</tt> this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_space_ptr_t& space_GcUD() const;
	///
	/** Return a matrix space object for a mutable matrix compatible with <tt>Gh(var_dep)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * This matrix space object is designed to create mutable matrix objects compatible
	 * with <tt>Gh(var_dep)</tt>.  For example, a matrix object <tt>Q</tt> created by this matrix space
	 * can be used to compute <tt>Q = Gh(var_dep)' - Gh(var_indep)'*D'</tt> (this is needed
	 * by a orthogonal range/null decomposition.
	 *
	 * The default implementation is to return <tt>return.get() == NULL</tt>.
	 * This is the correct implementation when <tt>mI() == 0</tt>.  However,
	 * when <tt>mI() > 0</tt> this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_space_ptr_t& space_GhD() const;

	//@}

	/** @name Calculation members */
	//@{

	///
	/** Compute all of the needed quanities for direct sensitivities.
	 *
	 *	@param	x	[in] (dim == n()) Current value of unkowns.  This vector should
	 *              have been created by <tt>this->space_x()->create_member()</tt>.
	 *	@param	f 	[out] Value of <tt>f(x)</tt>.
	 *				If f == NULL then this quantity is not computed.
	 *	@param	c 	[out] (dim == m()) Value of the equality constraints <tt>c(x)</tt>.
	 *				If c == NULL then this quantity is not computed.
	 *				If c != NULL and recalc_c == true on input then this quantity is
	 *				not recomputed and is used in the computation of
	 *				py if requested (i.e. py != NULL).  If <tt>!=NULL</tt> this this vector
	 *              should have been created by <tt>this->space_c()->create_member()</tt>.
	 *	@param  recalc_c
	 *				[in] If true then c(x) will be recomputed at x.
	 *              If false then <tt>c</tt> will not be recomputed but will be used as stated above.
	 *  @param  h   [out] (dim == mI()) Value of the general inequality constraints h(x).
	 *              If mI() == 0 then <tt>h</tt> should be set to <tt>NULL</tt> on input.
	 *              If h == NULL then this value will not be computed.  If <tt>!=NULL</tt> this
	 *              this vector should have been created by <tt>this->space_h()->create_member()</tt>.
	 *
	 *	@param	Gf	[out] (dim == n()) Gradient of <tt>f(x)</tt>.
	 *				If Gf == NULL then this quantity is not computed.  If <tt>!=NULL</tt> this this vector
	 *              should have been created by <tt>this->space_x()->create_member()</tt>.
	 *	@param	py
	 *				[out] (dim == r()) <tt>py = -inv(C)*c(con_decomp)</tt>
	 *				If py == NULL then this quantity is not computed.
	 *				If recalc_c == false on input then this c may be used
	 *				in the computation of py.  If <tt>!=NULL</tt> this this vector should have been
	 *              created by <tt>this->space_x()->sub_space(this->var_dep())->create_member()</tt>.
	 *	@param	rGf
	 *				[out] (dim == n()-r()) <tt>rGf = Gf(1,n()-r()) + D'*Gf(n()-r()+1,n())</tt>
	 *              , which is the reduced gradient of the objective function projected
	 *              into the manifold of the decomposed equality constraints.  If <tt>NULL</tt>, this
	 *              vector is not computed.  If <tt>!=NULL</tt> this this vector
	 *              should have been created by <tt>this->space_x(this->var_indep())->create_member()</tt>.
	 *  @param  GcU [out] (dim = n x (m()-r())) Auxiliary jacobian matrix <tt>del(c(con_undecomp),x)</tt>.
	 *              If m() == r() then <tt>GcU</tt> should be set to <tt>NULL</tt> on input.
	 *              If GcU == NULL then this quantitiy is not computed.  If <tt>!=NULL</tt> this this matrix
	 *              should have been created by <tt>this->space_GcU()->create_member()</tt>.
	 *  @param  Gh  [out] (dim = n x (m()-r())) Auxiliary jacobian matrix <tt>del(h,x)</tt>.
	 *              If mI() == 0 then <tt>Gh</tt> should be set to <tt>NULL</tt> on input.
	 *              If Gh == NULL then this quantitiy is not computed.  If <tt>!=NULL</tt> this this matrix
	 *              should have been created by <tt>this->space_Gh()->create_member()</tt>.
	 *	@param	D   [out] (dim = r() x (n()-r())) <tt>D = -inv(C)*N</tt>, which is the direct
	 *              sensitivity of the constraints to the independent variables.
	 *				If D == NULL then this quantity is not computed.  If <tt>!=NULL</tt> this this matrix
	 *              should have been created by <tt>this->space_D()->create_member()</tt>.
	 *	@param	V   [out] (dim = r() x (n()-r())) <tt>V = F + E * D</tt>, which is the an
	 *              auxiliary sensitivity matrix.  If m() == r() then <tt>V</tt> should be set to
	 *              <tt>NULL</tt> on input.  If V == NULL then this quantity is not computed.
	 *              If <tt>!=NULL</tt> this this matrix should have been created by
	 *              <tt>this->space_V()->create_member()</tt>.
	 *	@param	P   [out] (dim = mI() x (n()-r())) <tt>P = GhI' + GhD' * D</tt>, which is the an
	 *              auxiliary sensitivity matrix.  If mI() == 0 then <tt>P</tt> should be set to
	 *              <tt>NULL</tt> on input.  If P == NULL then this quantity is not computed.
	 *              If <tt>!=NULL</tt> this this matrix should have been created by
	 *              <tt>this->space_P()->create_member()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> Lots more!
	 * </ul>
	 */
	virtual void calc_point(
		const VectorWithOp     &x
		,value_type            *f
		,VectorWithOpMutable   *c
		,bool                  recalc_c
		,VectorWithOpMutable   *h
		,VectorWithOpMutable   *Gf
		,VectorWithOpMutable   *py
		,VectorWithOpMutable   *rGf
		,MatrixWithOp          *GcU
		,MatrixWithOp          *Gh
		,MatrixWithOp          *D
		,MatrixWithOp          *V
		,MatrixWithOp          *P
		) const = 0;

	///
	/** Calculate an approximate newton step given the Jacobian computed
	 * for the last call to calc_point(,,,).
	 *
	 * The idea behind this method is that with some applications it may be
	 * much cheaper to compute an approximate Newton step for the constraints
	 * given information computed during the last call to \Ref{calc_point}<tt>(...)</tt>.
	 * It is assumed that this approximate solution <tt>py</tt> will still be a
	 * descent direction for <tt>c(x)</tt>.  Some subclasses may have to perform an equal
	 * amount of work as <tt>calc_point(...)</tt> to perform this calculation but those
	 * are the breaks.
	 *
	 *	@param	x	[in] (dim == n()) current value of unkowns.
	 *	@param	c 	[out] (dim == m()) Value of the constraints c(x)
	 *				If c == NULL then this quantity is not computed.
	 *				If c != NULL and recalc_c == true on input then this quantity is
	 *				not recomputed and is used in the computation of
	 *				py if requested (i.e. py!=NULL).
	 *	@param  recalc_c
	 *	@param	py
	 *				[out] (size == r() on output) Approximate value of -inv(C)*c
	 *				Note that py == NULL is not allowed here.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> Lots more.
	 * </ul>
	 */
	virtual void calc_semi_newton_step(
		const VectorWithOp    &x
		,VectorWithOpMutable  *c
		,bool                 recalc_c
		,VectorWithOpMutable  *py
		) const = 0;

	//@}

	/** @name Overridden from NLP */
	//@{

	/// Returns 0
	size_type mI() const;
	/// Returns <tt>return.get() == NULL</tt>.
	vec_space_ptr_t space_h() const;
	/// Throws exception.
	const VectorWithOp& hl() const;
	/// Throws exception.
	const VectorWithOp& hu() const;

	//@}

protected:

	/** @name Overridden from NLP */
	//@{
	
	/// Does nothing (should never be called though).
	void imp_calc_h(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const;

	//@}

};	// end class NLPFirstOrderDirect

}	// end namespace NLPInterfacePack

#endif   // NLP_FIRST_ORDER_DIRECT_H
