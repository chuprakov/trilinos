// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrderInfo.h

#ifndef NLP_FIRST_ORDER_INFO_H
#define NLP_FIRST_ORDER_INFO_H

#include "NLPObjGradient.h"

namespace NLPInterfacePack {
///
/** NLP first order information interface class {abstract}.
  *
  * This class adds Jacobian information for the constraints 
  * #Gf#, #f# and #c# info from \Ref{NLPObjGradient}.
  */
class NLPFirstOrderInfo : public NLPObjGradient {
public:

	/** @name Public Types */
	//@{

	/// Thrown if an invalid matrix type is passed in.
	class InvalidMatrixType : public std::logic_error
	{public: InvalidMatrixType(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPFirstOrderInfo();

	//@}

	///
	/** Initialize the NLP for its first use.
	  */
	void initialize();

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

	/** @name Calculation Members.
	  *
	  * See \Ref{NLP} for a description of the behavior of these member functions.
	  */

	//@{

	///
	/** Update the matrix for #Gc# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c# and #Gf# may also be changed
	  * but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_Gc(const VectorSlice& x, bool newx = true) const;

	//@}

	///
	/** Number of constraint function Jacobian evaluations.
	  *
	  * This function can be called to find out how many evaluations
	  * the client requested since initialize() was called.
	  */
	virtual size_type num_Gc_evals() const;

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
		FirstOrderInfo( MatrixWithOp* Gc_in, const ObjGradInfo& obj_grad )
			: Gc(Gc_in), Gf(obj_grad.Gf), f(obj_grad.f), c(obj_grad.c)
		{}
		/// Pointer to gradient of objective function #Gf# (may be NULL if not set)
		MatrixWithOp* Gc;
		/// Pointer to gradient of objective function #Gf# (may be NULL if not set)
		Vector*       Gf;
		/// Pointer to objective function #f# (may be NULL if not set)
		value_type*   f;
		/// Pointer to constraints residule #c# (may be NULL if not set)
		Vector*       c;
	}; // end struct FirstOrderInfo

	/// Return objective gradient and zero order information.
	const FirstOrderInfo first_order_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute Gc(x) and perhaps G(x), f(x) and c(x) (if multiple calculaiton = true).
	 *
	 * @param x           Unknown vector (size n).
	 * @param newx        True if is a new point.
	 * @param obj_grad    Pointers to Gf, f and c
	 */
	virtual void imp_calc_Gc(const VectorSlice& x, bool newx, const FirstOrderInfo& first_order_info) const = 0;

	//@}

private:
	mutable MatrixWithOp      *Gc_;
	mutable size_type         num_Gc_evals_;

};	// end class NLPFirstOrderInfo

// /////////////////////
// Inline members

inline
const NLPFirstOrderInfo::FirstOrderInfo NLPFirstOrderInfo::first_order_info() const
{
	return FirstOrderInfo(Gc_,obj_grad_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_FIRST_ORDER_INFO_H
