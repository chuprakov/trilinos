// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrderInfo.h
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

#ifndef NLP_FIRST_ORDER_INFO_H
#define NLP_FIRST_ORDER_INFO_H

#include "NLPObjGradient.h"

namespace NLPInterfacePack {
///
/** NLP first order information interface class {abstract}.
 *
 * <b>Overview:</b>
 *
 * This class adds Jacobian information for the constraints.  This augments the
 * information provided by the \c NLP and \c NLPObjGradient interfaces.  This interface
 * includes access to matrix space objects which must be used by the client
 * to create the matrix objects that are used with this interface.  This
 * totally decouples the client from the implementation of these matrix
 * objects.
 *
 * <b>Client Usage:</b>
 *
 * As with the <tt>NLP</tt> base interface, the <tt>initialize()</tt> method must be called before
 * the %NLP object can be used.
 *
 * The matrix space objects returned from \c factory_Gc() and \c factory_Gh() must be used to create
 * the matrix objects for \c Gc and \c Gh used with this interface.  Note that the matrix objects
 * returned from these matrix space objects can not be expected to be usable until they are
 * passed to the calculation routines.
 *
 * The methods \c set_Gc() and \c set_Gh are used to set a pointers to matrix objects to be updated
 * when \c Gc and \c Gh are computed using \c calc_Gc() and \c calc_Gh().
 *
 * The number of evaluations of \c Gc and \c Gh using calc_Gc() and calc_Gh() are returned by
 * \c num_Gc_evals() and \c num_Gh_evals().
 * 
 * <b>Subclass developer's notes:</b>
 *
 * In addition to the methods that must be overridden by the \c NLPObjGradient interface
 * (<A HREF="classNLPInterfacePack_1_1NLPObjGradient.html#must_override">see</A>) the following methods
 * must also be overridden: \c factory_Gc(), \c factory_Gh(), \c imp_calc_Gc(), \c imp_calc_Gc().
 *
 * In addition to the methods that should be overridden from <tt>%NLPObjGradient</tt> by most subclasses
 * (<A HREF="classNLPInterfacePack_1_1NLPObjGradient.html#should_override">see</A>), the following
 * additional methods should be overridden: \c initialize().
 *
 * The following methods should never have to be overridden by most subclasses except in some very
 * specialized situations: \c set_Gc(), \c get_Gc(), \c Gc(), \c set_Gh(), \c get_Gh(), \c Gh(),
 * \c num_Gc_evals(), \c num_Gh_evals().
 */
class NLPFirstOrderInfo : virtual public NLPObjGradient {
public:

	///
	typedef MemMngPack::ref_count_ptr<
		const MemMngPack::AbstractFactory<MatrixWithOp> >        mat_fcty_ptr_t;
	///
	typedef MemMngPack::ref_count_ptr<BasisSystem>               basis_sys_ptr_t;

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPFirstOrderInfo();

	//@}

	/** @name NLP initialization */
	//@{

	///
	/** Initialize the NLP for its first use.
	 *
	 * This function implementation should be called by subclass implementations
	 * in order to reset counts for \c f(x), \c c(x), \c h(x), \c Gf(x), \c Gc(x)
	 * and \c Gh(x) evaluations.  This implementation calls
	 * <tt>this->NLPObjGradient::initialize()</tt>
	 *
	 * Postconditions:<ul>
	 * <li> See <tt>NLPObjGradient::initialize()</tt>
	 * <li> <tt>this->num_Gc_evals() == 0</tt>
	 * <li> <tt>this->num_Gh_evals() == 0</tt>
	 * </ul>
	 */
	void initialize(bool test_setup);

	//@}

	/** @name Matrix factory objects */
	//@{
	
	///
	/** Return a matrix factory object for creating <tt>Gc</tt>.
	 *
	 * This method may return <tt>return.get() == NULL</tt> if <tt>m() == 0</tt>.
	 * Otherwise, it must must return a valid matrix factory object.
	 */
	virtual const mat_fcty_ptr_t factory_Gc() const = 0;

	///
	/** Return a matrix factory object for creating <tt>Gh</tt>.
	 *
	 * This method may return <tt>return.get() == NULL</tt> if <tt>mI() == 0</tt>.
	 * Otherwise, it must must return a valid matrix factory object.
	 */
	virtual const mat_fcty_ptr_t factory_Gh() const = 0;

	//@}

	/** @name BasisSystem */
	//@{

	///
	/** Return a <tt>BasisSystem</tt> object compatible with <tt>Gc</tt> and <tt>Gh</tt>.
	 *
	 * Note that multiple calls to this method may return the same <tt>return.get()</tt>
	 * value so the client must not assume that they are unique.
	 *
	 * The default implementation returns <tt>return.get() == NULL</tt>.
	 */
	virtual const basis_sys_ptr_t basis_sys() const;

	//@}

	/** @name <<std aggr>> members for the gradient of the objective function Gc(x) */
	//@{

	///
	/** Set a pointer to a matrix object to be updated when <tt>this->calc_Gc()</tt> is called.
	 *
	 * @param  Gc  [in] Pointer to matrix of gradients.  May be \c NULL.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>Gc != NULL</tt>] <tt>Gc->space().is_compatible(*this->factory_Gc(),no_trans) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_Gc() == Gc</tt>
	 * </ul>
	 */
	virtual void set_Gc(MatrixWithOp* Gc);
	///
	/** Return pointer passed to <tt>this->set_Gc()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual MatrixWithOp* get_Gc();
	///
	/** Returns non-<tt>const</tt> <tt>*this->get_Gc()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gc() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual MatrixWithOp& Gc();
	///
	/** Returns <tt>const</tt> <tt>*this->get_Gc()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gc() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual const MatrixWithOp& Gc() const;

	//@}

	/** @name <<std aggr>> members for the gradient of the objective function Gh(x) */
	//@{

	///
	/** Set a pointer to a matrix object to be updated when <tt>this->calc_Gh()</tt> is called.
	 *
	 * @param  Gh  [in] Pointer to matrix of gradients.  May be \c NULL.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>Gh != NULL</tt>] <tt>Gh->space().is_compatible(*this->factory_Gh(),no_trans) == true</tt>
	 *      (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_Gh() == Gh</tt>
	 * </ul>
	 */
	virtual void set_Gh(MatrixWithOp* Gh);
	///
	/** Return pointer passed to <tt>this->set_Gh()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual MatrixWithOp* get_Gh();
	///
	/** Returns non-<tt>const</tt> <tt>*this->get_Gh()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gh() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual MatrixWithOp& Gh();
	///
	/** Returns <tt>const</tt> <tt>*this->get_Gh()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>this->get_Gh() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 */
	virtual const MatrixWithOp& Gh() const;

	//@}


	/** @name Calculation Members. */
	//@{

	///
	/** Update the matrix for \c Gc at the point \c x and put it in the stored reference.
	 *
	 * @param  x     [in] Point at which to calculate the matrix of gradients <tt>Gc(x)</tt>.
	 * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
	 *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
	 *               If \c false, the values in \c x are not the same as the last call to a
	 *               <tt>this->calc_*(x,newx)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>this->get_Gc() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->Gc()</tt> is updated to \c Gc(x)
	 * </ul>
	 *
	 * If \c set_multi_calc(true) was called then referenced storage for \c f, \c c, \c Gf and \c Gh may also be changed
	 * but are not guarantied to be.  But no other quanities from possible subclasses are allowed
	 * to be updated as a side effect.
	 */ 
	virtual void calc_Gc(const VectorWithOp& x, bool newx = true) const;

	///
	/** Update the matrix for \c Gh at the point \c x and put it in the stored reference.
	 *
	 * @param  x     [in] Point at which to calculate the matrix of gradients <tt>Gh(x)</tt>.
	 * @param  newx  [in] (default \c true) If \c true, the values in \c x are the same as
	 *               the last call to a <tt>this->calc_*(x,newx)</tt> member.
	 *               If \c false, the values in \c x are not the same as the last call to a
	 *               <tt>this->calc_*(x,newx)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->is_initialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> <tt>x.space().is_compatible(*this->space_x()) == true</tt> (throw <tt>VectorSpace::IncompatibleVectorSpaces</tt>)
	 * <li> <tt>this->get_Gh() != NULL</tt> (throw <tt>NoRefSet</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->Gh()</tt> is updated to \c Gh(x)
	 * </ul>
	 *
	 * If \c set_multi_calc(true) was called then referenced storage for \c f, \c c, \c Gf and \c Gc may also be changed
	 * but are not guarantied to be.  But no other quanities from possible subclasses are allowed
	 * to be updated as a side effect.
	 */ 
	virtual void calc_Gh(const VectorWithOp& x, bool newx = true) const;

	//@}

	/** @name Function evaluation counts. */
	//@{

	///
	/** Gradient of constraints matrix \c Gc evaluations count.
	 *
	 * This function can be called to find out how many evaluations
	 * \c this->calc_Gc() the client requested since \c this->initialize() was called.
	 */
	virtual size_type num_Gc_evals() const;

	///
	/** Gradient of constraints matrix \c Gh evaluations count.
	 *
	 * This function can be called to find out how many evaluations
	 * \c this->calc_Gh() the client requested since \c this->initialize() was called.
	 */
	virtual size_type num_Gh_evals() const;

	//@}

protected:

	///
	/** Struct for zero and first order quantities (pointers)
	 */
	struct FirstOrderInfo {
		///
		FirstOrderInfo()
			: Gc(NULL), Gf(NULL), f(NULL), c(NULL)
		{}
		///
		FirstOrderInfo( MatrixWithOp* Gc_in, MatrixWithOp* Gh_in, const ObjGradInfo& obj_grad )
			: Gc(Gc_in), Gh(Gh_in), Gf(obj_grad.Gf), f(obj_grad.f), c(obj_grad.c), h(obj_grad.h)
		{}
		/// Pointer to Jacobian of equality constraints <tt>Gc</tt> (may be NULL if not set)
		MatrixWithOp*           Gc;
		/// Pointer to Jacobian of inequality constraints <tt>Gh</tt> (may be NULL if not set)
		MatrixWithOp*           Gh;
		/// Pointer to gradient of objective function <tt>Gf</tt> (may be NULL if not set)
		VectorWithOpMutable*    Gf;
		/// Pointer to objective function <tt>f</tt> (may be NULL if not set)
		value_type*             f;
		/// Pointer to equality constraints residule <tt>c</tt> (may be NULL if not set)
		VectorWithOpMutable*    c;
		/// Pointer to inequality constraints residule <tt>h</tt> (may be NULL if not set)
		VectorWithOpMutable*    h;
	}; // end struct FirstOrderInfo

	/// Return objective gradient and zero order information.
	const FirstOrderInfo first_order_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute \a Gc(x) and perhaps \a Gh(x), \a Gf(x), \a f(x) and \a c(x).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
	 * <li> <tt>obj_grad_info.Gc != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*obj_grad_info.Gc</tt> is updated to \a Gc(x).
	 * </ul>
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param obj_grad_info
	 *                [out] Pointers to \c f, \c c, \c h, \c Gf, \c Gc and \c Gh.
	 *                On output <tt>*obj_grad_info.Gc</tt> is updated to \a Gc(x).
	 *                Any of the other objects pointed to in
	 *                \c obj_grad_info may be set if <tt>this->multi_calc() == true</tt> but are
	 *                now guaranteed to be.
	 */
	virtual void imp_calc_Gc(const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const = 0;

	///
	/** Overridden to compute \a Gh(x) and perhaps \a Gh(x), \a Gf(x), \a f(x) and \a c(x).
	 *
	 * @param x       [in]  Unknown vector (size n).
	 * @param newx    [in]  True if is a new point.
	 * @param obj_grad_info
	 *                [out] Pointers to \c f, \c c, \c h, \c Gf, \c Gc and \c Gh.
	 *                On output <tt>*obj_grad_info.Gh</tt> is updated to \a Gh(x).
	 *                Any of the other objects pointed to in
	 *                \c obj_grad_info may be set if <tt>this->multi_calc() == true</tt> but are
	 *                now guaranteed to be.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>x.space().is_compatible(*this->space_x())</tt> (throw <tt>IncompatibleType</tt>)
	 * <li> <tt>obj_grad_info.Gh != NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*obj_grad_info.Gh</tt> is updated to \a Gh(x).
	 * </ul>
	 */
	virtual void imp_calc_Gh(const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const = 0;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	MemMngPack::AbstractFactory<AbstractLinAlgPack::MatrixWithOp>  *factory_Gc;
	MemMngPack::AbstractFactory<AbstractLinAlgPack::MatrixWithOp>  *factory_Gh;
#endif
	mutable MatrixWithOp      *Gc_;
	mutable MatrixWithOp      *Gh_;
	mutable size_type         num_Gc_evals_;
	mutable size_type         num_Gh_evals_;

};	// end class NLPFirstOrderInfo

// /////////////////////
// Inline members

inline
const NLPFirstOrderInfo::FirstOrderInfo NLPFirstOrderInfo::first_order_info() const
{
	return FirstOrderInfo(Gc_,Gh_,obj_grad_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_FIRST_ORDER_INFO_H
