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
#include "SparseLinAlgPack/include/MatrixSymWithOp.h"

namespace NLPInterfacePack {

///
/** NLP second order information interface class {abstract}.
  *
  * This class adds second order inforamtion to the first order information
  * and basic information given in the \Ref{NLPFirstOrderInfo} interface class.
  *
  * Specifically the Hesssian of the Lagrangian is defined as:
  *
  * HL = Hf + sum( Hcj * lambda(j), j = 1...m )
  *
  * Where:\\
  * Hf is the hessian of the objective function\\
  * Hcj is the hessian of the jth equality constriant cj\\
  * lambda is the vector of lagrange multipliers for the equality
  * constraints c(x)
  */
class NLPSecondOrderInfo : public NLPFirstOrderInfo {
public:

	/** @name Public types */
	//@{

	///
	typedef	NLPFirstOrderInfo::InvalidMatrixType	InvalidMatrixType;

	//@}


	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPSecondOrderInfo();

	//@}
	
	///
	void initialize();

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
	/** Update the matrix for #HL# at the point #x#, #lambda# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c#, #Gf#, and #Gc# may also be
	  * changed but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_HL(const VectorSlice& x, const VectorSlice& lambda, bool newx = true) const;

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
			: Gc(NULL), Gf(NULL), f(NULL), c(NULL)
		{}
		///
		SecondOrderInfo( MatrixSymWithOp* HL_in, const FirstOrderInfo& first_order_info )
			: HL(HL_in), Gc(first_order_info.Gc), Gf(first_order_info.Gf)
			, f(first_order_info.f), c(first_order_info.c)
		{}
		/// Pointer to Hessiand of the Lagrangian #HL) (may be NULL is not set)
		MatrixSymWithOp* HL;
		/// Pointer to gradient of objective function #Gf# (may be NULL if not set)
		MatrixWithOp*    Gc;
		/// Pointer to gradient of objective function #Gf# (may be NULL if not set)
		Vector*          Gf;
		/// Pointer to objective function #f# (may be NULL if not set)
		value_type*      f;
		/// Pointer to constraints residule #c# (may be NULL if not set)
		Vector*          c;
	}; // end struct SecondOrderInfo

	/// Return objective gradient and zero order information.
	const SecondOrderInfo second_order_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute Gc(x) and perhaps G(x), f(x) and c(x) (if multiple calculaiton = true).
	 *
	 * @param x                     [in] Unknown vector (size n).
	 * @param lambda                [in] Lagrange multipliers for equality constraints
	 * @param newx                  [in] True if is a new point.
	 * @param second_order_info     [out] Pointers to HL, Gc, Gf, f and c
	 */
	virtual void imp_calc_HL(const VectorSlice& x, const VectorSlice& lambda, bool newx
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
