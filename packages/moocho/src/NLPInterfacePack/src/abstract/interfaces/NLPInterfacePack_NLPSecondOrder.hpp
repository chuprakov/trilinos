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

	/** @name <<std comp>> members for the Hessian if the objective function value, Hf.
	  *
	  * The is a reference to a matrix storage location that is updated when
	  * #calc_Hf(...)# is called or when #Hf# is updated as a side effect
	  * when another "calc" member funciton is called.  The #same_struct# argument
	  * is ment to be used by subclasses that override #set_Gc(...)# to exploit the
	  * structure of the currently set Gc matrix.
	  */
	//@{

	///
	virtual void set_Hf(MatrixWithOp* Hf, bool owns_Hf = false, bool same_struct = true);
	///
	virtual MatrixWithOp* get_Hf();
	///
	virtual void set_owns_Hf(bool owns_Hf);
	///
	virtual bool owns_Hf();
	///
	virtual MatrixWithOp& Hf();
	///
	virtual const MatrixWithOp& Hf() const;

	//@}

	/** @name <<std comp>> members for the Hessian if the jth equality constaint function, Hcj.
	  *
	  * The is a reference to a matrix storage location that is updated when
	  * #calc_Hcj(...)# is called or when #Hcj# is updated as a side effect
	  * when another "calc" member funciton is called.  The #same_struct# argument
	  * is ment to be used by subclasses that override #set_Hcj(...)# to exploit the
	  * structure of the currently set Hcj matrix.
	  */
	//@{

	///
	virtual void set_Hcj(MatrixWithOp* Hcj, bool owns_Hcj = false, bool same_struct = true);
	///
	virtual MatrixWithOp* get_Hcj();
	///
	virtual void set_owns_Hcj(bool owns_Hcj);
	///
	virtual bool owns_Hcj();
	///
	virtual MatrixWithOp& Hcj();
	///
	virtual const MatrixWithOp& Hcj() const;

	//@}

	//@}

	/** @name Calculation Members.
	  *
	  * See \Ref{NLP} for a description of the behavior of these member functions.
	  */

	//@{

	/// Return if #Hf# is calcuated analytically (true) or by some numerical technique (false)
	virtual bool analytic_Hf() const = 0;
	/// Return if #Hcj# is calcuated analytically (true) or by some numerical technique (false)
	virtual bool analytic_Hcj() const = 0;

	///
	/** Update the matrix for #Hf# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c#, #Gf#, and #Gc# may also be
	  * changed but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_Hf(const Vector& x, bool newx = true) const;
	///
	/** Update the matrix for #Hcj# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c#, #Gf#, and #Gc# may also be
	  * changed but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_Hcj(const Vector& x, size_type j, bool newx = true) const;

	//@}

protected:

	/** @name Protected methods to be overridden by subclasses */
	//@{

	/// Override to update the protected member #Gf_#
	virtual void imp_calc_Hf(const Vector& x, bool newx) const = 0;
	/// Override to update the protected member #Gc_#
	virtual void imp_calc_Hcj(const Vector& x, size_type j, bool newx) const = 0;

	//@}

	/** @name Protected members to be updated by subclasses */
	//@{

	/// Updated by subclass that implements #imp_calc_Hf(x,newx)#.
	mutable MatrixWithOp*			Hf_;
	/// Updated by subclass that implements #imp_calc_Hcj(x,j,newx)#.
	mutable MatrixWithOp*			Hcj_;

	//@}

private:
	mutable bool					owns_Hf_;
	static const char				name_Hf_[];
	mutable bool					owns_Hcj_;
	static const char				name_Hcj_[];

};	// end class NLPSecondOrderInfo

}	// end namespace NLPInterfacePack 

#endif // NLP_SECOND_ORDER_INFO_H