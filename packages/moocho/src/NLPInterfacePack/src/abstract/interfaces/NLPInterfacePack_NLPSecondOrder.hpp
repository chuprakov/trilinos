// ///////////////////////////////////////////////////////////////////////
// NLPSecondOrderInfo.h

#ifndef NLP_SECOND_ORDER_INFO_H
#define NLP_SECOND_ORDER_INFO_H

#include "NLPFirstOrderInfo.h"

namespace NLPInterfacePack {

///
/** NLP Second Order Information Interface Class.
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

	/// Clean up all of the memory #this# owns
	~NLPSecondOrderInfo();

	//@}
	
	///
	void initialize();

	/** @name <<std comp>> stereotype member functions
	  */
	//@{

	/** @name <<std comp>> members for the Hessian if the objective Lagrangian function
	  * value, HL.
	  *
	  * The is a reference to a matrix storage location that is updated when
	  * #calc_HL(...)# is called.
	  */
	//@{

	///
	virtual void set_HL(MatrixWithOp* HL, bool owns_HL = false);
	///
	virtual MatrixWithOp* get_HL();
	///
	virtual void set_owns_HL(bool owns_HL);
	///
	virtual bool owns_HL();
	///
	virtual MatrixWithOp& HL();
	///
	virtual const MatrixWithOp& HL() const;

	//@}

	//@}

	/** @name Calculation Members.
	  *
	  * See \Ref{NLP} for a description of the behavior of these member functions.
	  */

	//@{

	/// Return if #HL# is calcuated analytically (true) or by some numerical technique (false)
	virtual bool analytic_HL() const = 0;

	///
	/** Update the matrix for #HL# at the point #x#, #lambda# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c#, #Gf#, and #Gc# may also be
	  * changed but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_HL(const VectorSlice& x, const VectorSlice& lambda
		, bool newx = true) const;

	//@}

protected:

	/** @name Protected methods to be overridden by subclasses */
	//@{

	/// Override to update the protected member #HL_#
	virtual void imp_calc_HL(const VectorSlice& x, const VectorSlice& lambda
		, bool newx) const = 0;

	//@}

	/** @name Protected members to be updated by subclasses */
	//@{

	/// Updated by subclass that implements #imp_calc_HL(x,lambda,newx)#.
	mutable MatrixWithOp*			HL_;

	//@}

private:
	mutable bool					owns_HL_;
	static const char				name_HL_[];

};	// end class NLPSecondOrderInfo

}	// end namespace NLPInterfacePack 

#endif // NLP_SECOND_ORDER_INFO_H
