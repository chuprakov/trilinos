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
  * This class adds Jacobian information for the constraints.  This augments the
  * information provided by the \Ref{NLPObjGradient} interface.  This interface
  * includes access matrix space objects which must be used by the client
  * to create the matrix objects that are used with this interface.  This
  * totally decouples the client from the implementation of these matrix
  * objects.
  */
class NLPFirstOrderInfo : public NLPObjGradient {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		AbstractLinAlgPack::MatrixSpace<const MatrixWithOp> >    mat_space_ptr_t;

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPFirstOrderInfo();

	//@}

	///
	/** Initialize the NLP for its first use.
	  */
	void initialize();

	/** @name Get access to matrix space objects for the pertinent matrix spaces.
	 *
	 * The matrix space objects returned from these methods act as "Factories"
	 * that create the needed matrix objects.
	 */
	//@{
	
	///
	/** Return a matrix space object for creating #Gc#.
	 *
	 * This method may return #return.get() == NULL# if #m() == 0#.
	 * Otherwise, it must must return a valid matrix space object.
	 * The returned matrix object may not support the creation of any
	 * sub-matrix spaces (i.e. #return->sub_space(rrng,crng).get() == NULL#
	 * for all #rrng# and #crng#).
	 */
	virtual const mat_space_ptr_t& space_Gc() const = 0;

	///
	/** Return a matrix space object for creating #Gh#.
	 *
	 * This method may return #return.get() == NULL# if #mI() == 0#.
	 * Otherwise, it must must return a valid matrix space object.
	 * The returned matrix object may not support the creation of any
	 * sub-matrix spaces (i.e. #return->sub_space(rrng,crng).get() == NULL#
	 * for all #rrng# and #crng#).
	 */
	virtual const mat_space_ptr_t& space_Gh() const = 0;

	//@}

	/** @name <<std aggr>> members for the matrix of the gradient of equality constaints vector, Gc.
	  *
	  * The is a reference to a vector storage location that is updated when
	  * #calc_Gc(...)# is called or when #Gc# is updated as a side effect
	  * when another "calc" member funciton is called.
	  */
	//@{

	///
	virtual void set_Gc(MatrixWithOp* Gc);
	///
	virtual MatrixWithOp* get_Gc();
	///
	virtual MatrixWithOp& Gc();
	///
	virtual const MatrixWithOp& Gc() const;

	//@}

	/** @name <<std aggr>> members for the matrix of the gradient of equality constaints vector, Gh.
	  *
	  * The is a reference to a vector storage location that is updated when
	  * #calc_Gh(...)# is called or when #Gh# is updated as a side effect
	  * when another "calc" member funciton is called.
	  */
	//@{

	///
	virtual void set_Gh(MatrixWithOp* Gh);
	///
	virtual MatrixWithOp* get_Gh();
	///
	virtual MatrixWithOp& Gh();
	///
	virtual const MatrixWithOp& Gh() const;

	//@}

	/** @name Calculation Members.
	  *
	  * See \Ref{NLP} for a description of the behavior of these member functions.
	  */

	//@{

	///
	/** Update the matrix for #Gc# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c#, #Gf# and #Gh# may also be changed
	  * but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_Gc(const VectorWithOp& x, bool newx = true) const;

	///
	/** Update the matrix for #Gh# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c#, #Gf# and #Gc# may also be changed
	  * but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_Gh(const VectorWithOp& x, bool newx = true) const;

	//@}

	///
	/** Number of equality constraint function Jacobian evaluations.
	  *
	  * This function can be called to find out how many evaluations
	  * the client requested since initialize() was called.
	  */
	virtual size_type num_Gc_evals() const;

	///
	/** Number of inequality constraint function Jacobian evaluations.
	  *
	  * This function can be called to find out how many evaluations
	  * the client requested since initialize() was called.
	  */
	virtual size_type num_Gh_evals() const;

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
		/// Pointer to Jacobian of equality constraints #Gc# (may be NULL if not set)
		MatrixWithOp*           Gc;
		/// Pointer to Jacobian of inequality constraints #Gh# (may be NULL if not set)
		MatrixWithOp*           Gh;
		/// Pointer to gradient of objective function #Gf# (may be NULL if not set)
		VectorWithOpMutable*    Gf;
		/// Pointer to objective function #f# (may be NULL if not set)
		value_type*             f;
		/// Pointer to equality constraints residule #c# (may be NULL if not set)
		VectorWithOpMutable*    c;
		/// Pointer to inequality constraints residule #h# (may be NULL if not set)
		VectorWithOpMutable*    h;
	}; // end struct FirstOrderInfo

	/// Return objective gradient and zero order information.
	const FirstOrderInfo first_order_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute Gc(x) and perhaps Gh(x), Gf(x), f(x) and c(x) (if multiple calculaiton = true).
	 *
	 * @param x           Unknown vector (size n).
	 * @param newx        True if is a new point.
	 * @param obj_grad    Pointers to Gf, f, c and h
	 */
	virtual void imp_calc_Gc(const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const = 0;

	///
	/** Overridden to compute Gh(x) and perhaps Gh(x), Gf(x), f(x) and c(x) (if multiple calculaiton = true).
	 *
	 * @param x           Unknown vector (size n).
	 * @param newx        True if is a new point.
	 * @param obj_grad    Pointers to Gf, f, c and h
	 */
	virtual void imp_calc_Gh(const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const = 0;

	//@}

private:
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
