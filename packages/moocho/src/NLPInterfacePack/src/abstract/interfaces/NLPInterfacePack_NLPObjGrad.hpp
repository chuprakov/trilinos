// ///////////////////////////////////////////////////////////////////////
// NLPObjGradient.h

#ifndef NLP_OBJ_GRADIENT_H
#define NLP_OBJ_GRADIENT_H

#include "NLP.h"

namespace NLPInterfacePack {
///
/** NLP interface class that adds gradient information for the objective function {abstract}.
  *
  * This class adds the ability to compute the gradient of the objective function
  * Gf(x) to the basic information given in the \Ref{NLP} interface class.
  */
class NLPObjGradient : public NLP {
public:

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLPObjGradient();

	//@}

	///
	/** Initialize the NLP for its first use.
	  */
	void initialize();

	/** @name <<std aggr>> members for the gradient if the objective function value, Gf.
	  *
	  * The is a reference to a vector storage location that is updated when
	  * #calc_Gf(...)# is called or when #Gf# is updated as a side effect
	  * when another "calc" member funciton is called.
	  */
	//@{

	///
	virtual void set_Gf(Vector* Gf);
	///
	virtual Vector* get_Gf();
	///
	virtual Vector& Gf();
	///
	virtual const Vector& Gf() const;

	//@}

	/** @name Calculation Members.
	  *
	  * See \Ref{NLP} for a description of the behavior of these member functions.
	  */

	//@{

	///
	/** Update the vector for #Gf# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f# and #c#may also be updated
	  * but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect (i.e. no higher order derivatives).
	  */ 
	virtual void calc_Gf(const VectorSlice& x, bool newx = true) const;

	//@}

	///
	/** Objective gradients evaluation count.
	  *
	  * This function can be called to find out how many evaluations
	  * #calc_Gf(...)# the client requested since initialize() was called.
	  */
	virtual size_type num_Gf_evals() const;

	///
	/** Struct for gradient (objective), objective and constriants (pointers)
	 */
	struct ObjGradInfo {
		///
		ObjGradInfo()
			: Gf(NULL), f(NULL), c(NULL)
		{}
		///
		ObjGradInfo( Vector* Gf_in, const ZeroOrderInfo& f_c_in )
			: Gf(Gf_in), f(f_c_in.f), c(f_c_in.c)
		{}
		/// Pointer to gradient of objective function #Gf# (may be NULL if not set)
		Vector*       Gf;
		/// Pointer to objective function #f# (may be NULL if not set)
		value_type*   f;
		/// Pointer to constraints residule #c# (may be NULL if not set)
		Vector*       c;
	}; // end struct ObjGradInfo

protected:

	/// Return objective gradient and zero order information.
	const ObjGradInfo obj_grad_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute f(x) and perhaps c(x) (if multiple calculaiton = true).
	 *
	 * @param x           [in]  Unknown vector (size n).
	 * @param newx        [in]  True if is a new point.
	 * @param obj_grad    [out] Pointers to Gf, f and c
	 */
	virtual void imp_calc_Gf(const VectorSlice& x, bool newx, const ObjGradInfo& obj_grad_info) const = 0;

private:
	mutable Vector                  *Gf_;
	mutable size_type				num_Gf_evals_;

};	// end class NLPObjGradient

// //////////////////
// Inline members

inline
const NLPObjGradient::ObjGradInfo NLPObjGradient::obj_grad_info() const
{
	return ObjGradInfo(Gf_,zero_order_info());
}

}	// end namespace NLPInterfacePack 

#endif // NLP_OBJ_GRADIENT_H
