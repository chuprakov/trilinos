// ///////////////////////////////////////////////////////////////////////
// NLPSecondOrderInfo.h
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

#ifndef NLP_SECOND_ORDER_INFO_H
#define NLP_SECOND_ORDER_INFO_H

#include "NLPFirstOrderInfo.h"

namespace NLPInterfacePack {

///
/** NLP second order information interface class {abstract}.
  *
  * This class adds second order inforamtion to the first order information
  * and basic information given in the \Ref{NLPFirstOrderInfo} interface class.
  *
  * Specifically the Hesssian of the Lagrangian is defined as:
  *
  * HL = Hf + sum( Hcj * lambda(j), j = 1...m ) + sub( Hhj * lambdaI(j), j = 1...mI )
  *
  * Where:\\
  * Hf is the hessian of the objective function\\
  * Hcj is the hessian of the jth equality constriant cj\\
  * Hhj is the hessian of the jth inequality constriant hj\\
  * lambda is the vector of lagrange multipliers for the equality
  * constraints c(x) \\
  * lambdaI is the vector of lagrange multipliers for the inequality
  * constraints h(x) \\
  */
class NLPSecondOrderInfo : public NLPFirstOrderInfo {
public:

	///
	typedef ReferenceCountingPack::ref_count_ptr<
		AbstractLinAlgPack::MatrixSpace<const MatrixSymWithOp> >    mat_sym_space_ptr_t;

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPSecondOrderInfo();

	//@}
	
	///
	void initialize();

	///
	/** Return a matrix space object for creating #HL#.
	 *
	 * The returned matrix object may not support the creation of any
	 * sub-matrix spaces (i.e. #return->sub_space(rrng,crng).get() == NULL#
	 * for all #rrng# and #crng#).
	 */
	virtual const mat_sym_space_ptr_t& space_HL() const = 0;

	/** @name <<std aggr>> stereotype member functions
	  */
	//@{

	/** @name <<std aggr>> members for the Hessian if the objective Lagrangian function
	  * value, HL.
	  *
	  * The is a reference to a matrix storage location that is updated when
	  * #calc_HL(...)# is called.
	  */
	//@{

	///
	virtual void set_HL(MatrixSymWithOp* HL);
	///
	virtual MatrixSymWithOp* get_HL();
	///
	virtual MatrixSymWithOp& HL();
	///
	virtual const MatrixSymWithOp& HL() const;

	//@}

	//@}

	/** @name Calculation Members.
	  *
	  * See \Ref{NLP} for a description of the behavior of these member functions.
	  */

	//@{

	///
	/** Update the matrix for #HL# at the point #x#, #lambda#, #lambdaI# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c#, #h#, #Gf#, #Gc# and #Gh# may also be
	  * changed but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  *
	  * @param  x        [in] Unknown primal variables
	  * @param  lambda   [in] Lagrange muitipliers for equality constriants.
	  *                  If #m() == 0# then #lambda# can be #NULL#.
	  * @param  lambdaI  [in] Lagrange muitipliers for inequality constriants.
	  *                  If #mI() == 0# then #lambdaI# can be #NULL#.
	  * @param  newpoint [in] If true then this is the same x that was passed to the last
	  *                  call to a calc_info(...) method.
	  *
	  */ 
	virtual void calc_HL(
		const VectorWithOp& x, const VectorWithOp* lambda, const VectorWithOp* lambdaI, bool newpoint = true) const;

	//@}

	///
	/** Number of Hessian evaluations.
	  *
	  * This function can be called to find out how many evaluations
	  * the client requested since initialize() was called.
	  */
	virtual size_type num_HL_evals() const;

protected:

	///
	/** Struct for zero, first and second order quantities (pointers)
	 */
	struct SecondOrderInfo {
		///
		SecondOrderInfo()
			: HL(NULL), Gc(NULL), Gh(NULL), Gf(NULL), f(NULL), c(NULL), h(NULL)
			{}
		///
		SecondOrderInfo( MatrixSymWithOp* HL_in, const FirstOrderInfo& first_order_info )
			: HL(HL_in), Gc(first_order_info.Gc), Gh(first_order_info.Gh), Gf(first_order_info.Gf)
			, f(first_order_info.f), c(first_order_info.c), h(first_order_info.h)
			{}
		/// Pointer to Hessiand of the Lagrangian #HL) (may be NULL is not set)
		MatrixSymWithOp*        HL;
		/// Pointer to Hessian of the equality constraints #Gc# (may be NULL if not set)
		MatrixWithOp*           Gc;
		/// Pointer to Hessian of the inequality constraints #Gh# (may be NULL if not set)
		MatrixWithOp*           Gh;
		/// Pointer to gradient of objective function #Gf# (may be NULL if not set)
		VectorWithOpMutable*    Gf;
		/// Pointer to objective function #f# (may be NULL if not set)
		value_type*             f;
		/// Pointer to equality constraints residule #c# (may be NULL if not set)
		VectorWithOpMutable*    c;
		/// Pointer to inequality constraints residule #h# (may be NULL if not set)
		VectorWithOpMutable*    h;
	}; // end struct SecondOrderInfo

	/// Return objective gradient and zero order information.
	const SecondOrderInfo second_order_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute Gc(x) and perhaps G(x), f(x) and c(x) (if multiple calculaiton = true).
	 *
	 * @param x                     [in] Unknown vector (size n).
	 * @param lambda                [in] Lagrange multipliers for equality constraints c(x).
	 *                              My be #NULL# if #m() == 0#.
	 * @param lambdaI               [in] Lagrange multipliers for inequality constraints h(x).
	 *                              My be #NULL# if #mI() == 0#.
	 * @param newpoint              [in] True if is a new point.
	 * @param second_order_info     [out] Pointers to HL, Gc, Gh, Gf, f, c and h
	 */
	virtual void imp_calc_HL(
		const VectorWithOp& x, const VectorWithOp* lambda, const VectorWithOp* lambdaI, bool newpoint
		, const SecondOrderInfo& second_order_info) const = 0;

	//@}

private:
	mutable MatrixSymWithOp   *HL_;
	mutable bool              num_HL_evals_;

};	// end class NLPSecondOrderInfo

// //////////////////
// Inline members

inline
const NLPSecondOrderInfo::SecondOrderInfo NLPSecondOrderInfo::second_order_info() const
{
	return SecondOrderInfo(HL_,first_order_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_SECOND_ORDER_INFO_H
