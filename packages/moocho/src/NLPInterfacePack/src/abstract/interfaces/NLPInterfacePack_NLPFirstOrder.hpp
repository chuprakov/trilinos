// ///////////////////////////////////////////////////////////////////////
// NLPFirstOrderInfo.h

#ifndef NLP_FIRST_ORDER_INFO_H
#define NLP_FIRST_ORDER_INFO_H

#include "NLP.h"

namespace NLPInterfacePack {
///
/** NLP First Order Information Interface Class.
  *
  * This class adds first order inforamtion to the basic information
  * given in the \Ref{NLP} interface class.
  *
  * Given first order information it is possible for the NLP
  * solver to compute Lagrange multipliers and then report
  * these back to the NLP.
  *
  * The Lagrangian for this problem is defined by:\\
  \begin{verbatim}
	L = f(x) + lambda' * c(x) + nul * ( xl - x ) + nuu * ( x - xu )
  \end{verbatim}
  *
  * The optimality conditions are given by:
  \begin{verbatim}
	del(L,x)      = del(f,x) + del(c,x) * lambda + nu = 0
	del(L,lambda) = c(x) = 0
	  where:
		nu = nuu - nul
		nul(i) * ( xl(i) - x(i) ) = 0,      for i = 1...n
		nuu(i) * ( x(i) - xu(i) ) = 0,      for i = 1...n
  \end{verbatim}

  */
class NLPFirstOrderInfo : public NLP {
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

	/// Clean up any memory #this# owns
	~NLPFirstOrderInfo();

	//@}

	///
	/** Initialize the NLP for its first use.
	  */
	void initialize();
	
	/** @name <<std comp>> stereotype member functions
	  */
	//@{

	/** @name <<std comp>> members for the gradient if the objective function value, Gf.
	  *
	  * The is a reference to a vector storage location that is updated when
	  * #calc_Gf(...)# is called or when #Gf# is updated as a side effect
	  * when another "calc" member funciton is called.
	  */
	//@{

	///
	virtual void set_Gf(Vector* Gf, bool owns_Gf = false);
	///
	virtual Vector* get_Gf();
	///
	virtual void set_owns_Gf(bool owns_Gf);
	///
	virtual bool owns_Gf();
	///
	virtual Vector& Gf();
	///
	virtual const Vector& Gf() const;

	//@}

	/** @name <<std comp>> members for the matrix of the gradient of equality constaints vector, Gc.
	  *
	  * The is a reference to a vector storage location that is updated when
	  * #calc_Gc(...)# is called or when #Gc# is updated as a side effect
	  * when another "calc" member funciton is called.  The #same_struct# argument
	  * is ment to be used by subclasses that override #set_Gc(...)# to exploit the
	  * structure of the currently set Gc matrix.
	  */
	//@{

	///
	virtual void set_Gc(MatrixWithOp* Gc, bool owns_Gc = false, bool same_struct = true);
	///
	virtual MatrixWithOp* get_Gc();
	///
	virtual void set_owns_Gc(bool owns_Gc);
	///
	virtual bool owns_Gc();
	///
	virtual MatrixWithOp& Gc();
	///
	virtual const MatrixWithOp& Gc() const;

	//@}

	//@}

	/** @name Calculation Members.
	  *
	  * See \Ref{NLP} for a description of the behavior of these member functions.
	  */

	//@{

	/// Return if #Gf# is calcuated analytically (true) or by some numerical technique (false)
	virtual bool analytic_Gf() const = 0;
	/// Return if #Gc# is calcuated analytically (true) or by some numerical technique (false)
	virtual bool analytic_Gc() const = 0;
	///
	/** Update the vector for #Gf# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c# and #Gc# may also be changed
	  * but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_Gf(const Vector& x, bool newx = true) const;
	///
	/** Update the matrix for #Gc# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then referenced storage for #f#, #c# and #Gf# may also be changed
	  * but are not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_Gc(const Vector& x, bool newx = true) const;

	//@}

	///
	/** Report the final solution and multipliers.
	  *
	  * Call this function report the final solution of the
	  * unknows x and the Lagrange multipliers for the
	  * equality constriants #lambda# and the varaible bounds
	  * #nu#.  If either of the multipliers
	  * are not known then you can pass null in for them.
	  * The default action is to call report_final_x(x,optimal)
	  * on the NLP interface and then to ignore the multipliers.
	  */
	virtual void report_final_solution(
		  const VectorSlice&	x
		, const VectorSlice*	lambda
		, const SpVectorSlice*	nu
		, bool					optimal		) const;

	/** @name Objective and constraint function gradients.
	  *
	  * These functions can be called to find out how many evaluations
	  * the client requested since initialize() was called.
	  *
	  * Also, if the client calls calc_inf( x, false ) several times
	  * then the count will be incremented for each call even though
	  * the actual quantity may not actually be
	  * recalculated each time.  This information is not known here
	  * in this base class but the subclasses can overide this behavior
	  * if they want to.
	  */
	//@{

	///
	/** Gives the number of object function f(x) evaluations called by the solver
	  * since initialize() was called..
	  */
	virtual size_type num_Gf_evals() const;

	///
	/** Gives the number of constraint function c(x) evaluations called by the solver
	  * since initialize() was called..
	  */
	virtual size_type num_Gc_evals() const;

	//@}

protected:

	/** @name Protected methods to be overridden by subclasses */
	//@{

	/// Override to update the protected member #Gf_#
	virtual void imp_calc_Gf(const Vector& x, bool newx) const = 0;
	/// Override to update the protected member #Gc_#
	virtual void imp_calc_Gc(const Vector& x, bool newx) const = 0;

	//@}

	/** @name Protected members to be updated by subclasses */
	//@{

	/// Updated by subclass that implements #imp_calc_Gf(x,newx)#.
	mutable Vector*					Gf_;
	/// Updated by subclass that implements #imp_calc_Gc(x,newx)#.
	mutable MatrixWithOp*			Gc_;

	//@}

private:
	mutable bool					owns_Gf_;
	static const char				name_Gf_[];
	mutable size_type				num_Gf_evals_;
	mutable bool					owns_Gc_;
	static const char				name_Gc_[];
	mutable size_type				num_Gc_evals_;

};	// end class NLPFirstOrderInfo

}	// end namespace NLPInterfacePack 

#endif // NLP_FIRST_ORDER_INFO_H