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
/** Interface providing only directy first order sensitivity information.
  *
  * This interface defines a basis for the equality constriants and then only
  * certain linear systems with this basis are solved for.  This interface is
  * useful in reduced space SQP-type and other related algorithms.
  *
  * Specifically, the variables are partitioned into dependent and independent
  * sets #x = [ x_dep'  x_indep' ]'# and Jacobians of the constraints #c(x)#
  * and #h(x)# at the point #x# are:
  \begin{verbatim}
	del(c,x) = Gc' = [ del(c(con_decomp))   ] = [ GcD' ] = [ GcDD'  GcDI' ] = [ C  N ]
                     [ del(c(con_undecomp)) ]   [ GcU' ]   [ GcUD'  GcUI' ]   [ E  F ]

	del(h,x) = Gh' = [ GhD'  GhI' ]

    where:
      C <: R^(r x r) is nonsingular
      N <: R^(r x (n-r))
      E <: R^((m-r) x r)
      F <: R^((m-r) x (n-r))
  \end{verbatim}
  * This partitions the general equality constraints c(x) into two sets;
  * decomposed c(con_decomp) and undecomposed c(con_undecomp).  It is therefore
  * expected that sub-vectors and subspaces from #space_x().sub_space(var_dep)#,
  * #space_x().sub_space(var_indep)#, #space_c().sub_space(con_decomp)# and
  * #space_c().sub_space(con_undecomp)# can all be accessed.  Other sub-vectors
  * and sub-spaces may not be available (but the algorithm should not need
  * access to other sub-spaces).
  *
  * Free access to solves with the basis #C# is not given however and instead this interface
  * computes, for the current point x, the direct sensitivity matrices #D = -inv(C)*N#,
  * #V = F - E * D#, #P = GhI' + GhD'* D#, the auxiliary matrices
  * #GcU = [ GcUD; GcUI ] = [ E';  F' ]#
  * and #Gh = [ GhD;  GhI ]# and the Newton step #py = -inv(C)*c(con_decomp)#.
  * In general, linear solves with the transpose with #C# are not possible and
  * therefore are not avalible.  A number of very specialized applications can only
  * provide this information but this is all that is needed by many numerical
  * optimization (and related) algorithms.
  */
class NLPFirstOrderDirect : public NLPObjGradient
{
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractLinAlgPack::MatrixSpace<MatrixWithOp> >          mat_space_ptr_t;

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		const AbstractLinAlgPack::MatrixSpace<MatrixWithOpMutable> >   mat_mut_space_ptr_t;

	///
	/** Returns the number of decomposed equality constraints (r <= m).
	 *
	 * The default implementation returns #this->con_decomp().size()#.
	 * This implementation will work for all implementations.
	 */
	virtual size_type r() const;

	/** @name Return the ranges for dependent and independent variables and
	 * decomposed and undecomposed equalities.
	 *
	 * The default implementation is to return var_dep = [1,m], var_indep = [m+1,n]
	 * con_decomp = [1,m] and con_undecomp = [m+1,m+1].
	 */
	//@{

	///
	virtual Range1D var_dep() const;
	///
	virtual Range1D var_indep() const;
	///
	virtual Range1D con_decomp() const;
	///
	virtual Range1D con_undecomp() const;

	//@}

	/** @name Get access to matrix space objects for the pertinate matrix spaces.
	 */
	//@{
	
	///
	/** Return a matrix space object for creating #GcU#.
	 *
	 * The default implementation is to return #return.get() == NULL#.
	 * This is the proper implementation when #m() == r()#.
	 * When #m() > r()# then the subclass must override this method to
	 * return a valid matrix space object.  Moreover, the returned
	 * matrix object from #this->space_GcU()->create_member()->get_sub_view(rng,Range1D())#
	 * must be non-null for #rng == this->var_dep()# or #rng == this->var_indep()#.
	 * This gives access to the matrices #E'# and #F'# as shown above.
	 */
	virtual const mat_space_ptr_t& space_GcU() const;
	///
	/** Return a matrix space object for creating #Gh#.
	 *
	 * The default implementation is to return #return.get() == NULL#.
	 * This is the proper implementation when #m() == r()#.  Moreover, the returned
	 * matrix object from #this->space_Gh()->create_member()->get_sub_view(rng,Range1D())#
	 * must be non-null for #rng == this->var_dep()# or #rng == this->var_indep()#.
	 * This gives access to the matrices #GhD# and #GhI# as shown above.
	 */
	virtual const mat_space_ptr_t& space_Gh() const;
	///
	/** Return a matrix space object for #D = -inv(C)*N# {abstract}.
	 */
	virtual const mat_space_ptr_t& space_D() const = 0;
	///
	/** Return a matrix space object for #V = F + E * D#.
	 *
	 * The default implementation is to return #return.get() == NULL#.
	 * This is the correct implementation when #m() == r()#.  However,
	 * when #m() > r()# this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_space_ptr_t& space_V() const;
	///
	/** Return a matrix space object for #P = GhI' + GhD'* D#.
	 *
	 * The default implementation is to return #return.get() == NULL#.
	 * This is the correct implementation when #m() == r()#.  However,
	 * when #m() > r()# this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_space_ptr_t& space_P() const;
	///
	/** Return a matrix space object for a mutable matrix compatible with #GcU(var_dep)#.
	 *
	 * This matrix space object is designed to create mutable matrix objects compatible
	 * with #GcU(var_dep)#.  For example, a matrix object #U# created by this matrix space
	 * can be used to compute #U = GcU(var_dep)' - GcU(var_indep)'*D'# (this is needed
	 * by a orthogonal range/null decomposition.
	 *
	 * The default implementation is to return #return.get() == NULL#.
	 * This is the correct implementation when #m() == r()#.  However,
	 * when #m() > r()# this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_mut_space_ptr_t& space_GcUD() const;
	///
	/** Return a matrix space object for a mutable matrix compatible with #Gh(var_dep)#.
	 *
	 * This matrix space object is designed to create mutable matrix objects compatible
	 * with #Gh(var_dep)#.  For example, a matrix object #Q# created by this matrix space
	 * can be used to compute #Q = Gh(var_dep)' - Gh(var_indep)'*D'# (this is needed
	 * by a orthogonal range/null decomposition.
	 *
	 * The default implementation is to return #return.get() == NULL#.
	 * This is the correct implementation when #mI() == 0#.  However,
	 * when #mI() > 0# this method must be overridden to return a
	 * non-null matrix space object.
	 */
	virtual const mat_mut_space_ptr_t& space_GhD() const;

	//@}

	///
	/** Compute all of the needed quanities for direct sensitivities.
	 *
	 *	@param	x	[in] (dim == n()) Current value of unkowns.  This vector should
	 *              have been created by #this->space_x()->create_member()#.
	 *	@param	f 	[out] Value of #f(x)#.
	 *				If f == NULL then this quantity is not computed.
	 *	@param	c 	[out] (dim == m()) Value of the equality constraints #c(x)#.
	 *				If c == NULL then this quantity is not computed.
	 *				If c != NULL and recalc_c == true on input then this quantity is
	 *				not recomputed and is used in the computation of
	 *				py if requested (i.e. py != NULL).  If #!=NULL# this this vector
	 *              should have been created by #this->space_c()->create_member()#.
	 *	@param  recalc_c
	 *				[in] If true then c(x) will be recomputed at x.
	 *              If false then #c# will not be recomputed but will be used as stated above.
	 *  @param  h   [out] (dim == mI()) Value of the general inequality constraints h(x).
	 *              If mI() == 0 then #h# should be set to #NULL# on input.
	 *              If h == NULL then this value will not be computed.  If #!=NULL# this
	 *              this vector should have been created by #this->space_h()->create_member()#.
	 *
	 *	@param	Gf	[out] (dim == n()) Gradient of #f(x)#.
	 *				If Gf == NULL then this quantity is not computed.  If #!=NULL# this this vector
	 *              should have been created by #this->space_x()->create_member()#.
	 *	@param	py
	 *				[out] (dim == r()) #py = -inv(C)*c(con_decomp)#
	 *				If py == NULL then this quantity is not computed.
	 *				If recalc_c == false on input then this c may be used
	 *				in the computation of py.  If #!=NULL# this this vector should have been
	 *              created by #this->space_x()->sub_space(this->var_dep())->create_member()#.
	 *	@param	rGf
	 *				[out] (dim == n()-r()) #rGf = Gf(1,n()-r()) + D'*Gf(n()-r()+1,n())#
	 *              , which is the reduced gradient of the objective function projected
	 *              into the manifold of the decomposed equality constraints.  If #NULL#, this
	 *              vector is not computed.  If #!=NULL# this this vector
	 *              should have been created by #this->space_x(this->var_indep())->create_member()#.
	 *  @param  GcU [out] (dim = n x (m()-r())) Auxiliary jacobian matrix #del(c(con_undecomp),x)#.
	 *              If m() == r() then #GcU# should be set to #NULL# on input.
	 *              If GcU == NULL then this quantitiy is not computed.  If #!=NULL# this this matrix
	 *              should have been created by #this->space_GcU()->create_member()#.
	 *  @param  Gh  [out] (dim = n x (m()-r())) Auxiliary jacobian matrix #del(h,x)#.
	 *              If mI() == 0 then #Gh# should be set to #NULL# on input.
	 *              If Gh == NULL then this quantitiy is not computed.  If #!=NULL# this this matrix
	 *              should have been created by #this->space_Gh()->create_member()#.
	 *	@param	D   [out] (dim = r() x (n()-r())) #D = -inv(C)*N#, which is the direct
	 *              sensitivity of the constraints to the independent variables.
	 *				If D == NULL then this quantity is not computed.  If #!=NULL# this this matrix
	 *              should have been created by #this->space_D()->create_member()#.
	 *	@param	V   [out] (dim = r() x (n()-r())) #V = F + E * D#, which is the an
	 *              auxiliary sensitivity matrix.  If m() == r() then #V# should be set to
	 *              #NULL# on input.  If V == NULL then this quantity is not computed.
	 *              If #!=NULL# this this matrix should have been created by
	 *              #this->space_V()->create_member()#.
	 *	@param	P   [out] (dim = mI() x (n()-r())) #P = GhI' + GhD' * D#, which is the an
	 *              auxiliary sensitivity matrix.  If mI() == 0 then #P# should be set to
	 *              #NULL# on input.  If P == NULL then this quantity is not computed.
	 *              If #!=NULL# this this matrix should have been created by
	 *              #this->space_P()->create_member()#.
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
	 * The Idea behind this method is that with some applications it may be
	 * much cheaper to compute an approximate Newton step for the constraints
	 * given information computed during the last call to \Ref{calc_point}#(...)#.
	 * It is assumed that this approximate solution #py# will still be a
	 * descent direction for #c(x)#.  Some subclasses may have to perform an equal
	 * amount of work as #calc_point(...)# to perform this calculation but those
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
	 */
	virtual void calc_semi_newton_step(
		const VectorWithOp    &x
		,VectorWithOpMutable  *c
		,bool                 recalc_c
		,VectorWithOpMutable  *py
		) const = 0;

	/** @name Overridden from NLP.
	 *
	 * These methods are overridden for the case where the bassis #C# is
	 * full rank (#m() == r()#) and there are no general inequality constraints
	 * (#mI() == 0#).
	 */
	//@{

	/// Returns 0
	size_type mI() const;
	/// Returns #return.get() == NULL#.
	vec_space_ptr_t space_h() const;
	/// Throws exception.
	const VectorWithOp& hl() const;
	/// Throws exception.
	const VectorWithOp& hu() const;

	//@}

protected:

	/** @name Overridden from NLP.
	 *
	 * These methods are overridden for the case where the bassis #C# is
	 * full rank (#m() == r()#) and there are no general inequality constraints
	 * (#mI() == 0#).
	 */
	//@{
	
	/// Does nothing (should never be called though).
	void imp_calc_h(const VectorWithOp& x, bool newx, const ZeroOrderInfo& zero_order_info) const;

	//@}

};	// end class NLPFirstOrderDirect

}	// end namespace NLPInterfacePack

#endif   // NLP_FIRST_ORDER_DIRECT_H
