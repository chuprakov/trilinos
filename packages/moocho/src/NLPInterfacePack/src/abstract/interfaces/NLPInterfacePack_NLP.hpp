// ///////////////////////////////////////////////////////////////////////
// NLP.h

#ifndef NLP_H
#define NLP_H

#include <stdexcept>
#include <string>

#include "NLPInterfacePackTypes.h"
#include "Misc/include/StandardCompositionRelationshipsPack.h"

namespace NLPInterfacePack {

///
/** NLP interface class {abstract}.
  *
  * This class represents an abstract interface to a general nonlinear
  * programming problem of the form: \\
  * \\
  \begin{verbatim}
	min     f(x)
	s.t.    c(x) = 0
	        xl <= x <= xu
	where:
	        x    <: R^n
	        c(x) <: R^n -> R^m 
  \end{verbatim}
  * In the above form, none of the variables are fixed between bounds (strictly
  * xl < xu).
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
class NLP {
public:
	/** @name exceptions */
	//@{

	/// Thrown if any member functions are called before initialize() has been called.
	class UnInitialized : public std::logic_error
	{public: UnInitialized(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown from #initialize()# if some logical error occured
	class InvalidInitialization : public std::logic_error
	{public: InvalidInitialization(const std::string& what_arg) : std::logic_error(what_arg) {}};

	/// Thrown if #xl()# or #xu()# are called and #has_bounds() == false# 
	class NoBoundsOnVariables : public std::logic_error
	{public: NoBoundsOnVariables(const std::string& what_arg) : std::logic_error(what_arg) {}};

	//@}

	/** @name Constructors */
	//@{

	/// Initialize to no reference set to calculation quanities
	NLP();

	/// Destructor that cleans all the memory it owns
	virtual ~NLP();

	//@}
	
	///
	/** Initialize the NLP for its first use.
	  *
	  * Postconditions:  If this operation returns without throwing an exception then all of
	  * the member functions can be accessed.
	  *
	  * Note that subclasses must call this function to reset what needs to be
	  * reset here.
	  */
	virtual void initialize();

	///
	/** Return if this is initialized.
	  */
	virtual bool is_initialized() const = 0;
	
	/** @name Dimensionality. */
	//@{

	/// Return the number of variables
	virtual size_type n() const = 0;
	/// Return the number of equality constraints
	virtual size_type m() const = 0;
	/// Return the number of possibly decomposed equality constraints
	virtual size_type r() const = 0;	

	//@}

	/** @name Access Variable Bounds and the Initial Point */
	//@{

	///
	/** Returns true if there are bounds on the variables or false if not
	  * (conceptually equivalent to xl = -infinity, xu = infinity).
	  */
	virtual bool has_bounds() const = 0;
	///
	/** Set if the initial point must be within the bounds.
	  */
	virtual void set_force_xinit_in_bounds(bool force_xinit_in_bounds) = 0;
	///
	/** Returns if the initial point must be within the bounds.
	  */
	virtual bool force_xinit_in_bounds() const = 0;
	///
	/** Returns a reference to the vector of the initial guess for the solution #x#
	  */
	virtual const VectorSlice xinit() const = 0;
	///
	/** Returns a reference to the vector of lower bounds on the variables #x#.
	  *
	  * Preconditions: \begin{itemize}
	  * \item this->has_bounds() == true (throw NoBoundsOnVariables)
	  * \end{itemize}
	  */
	virtual const SpVectorSlice xl() const = 0;
	///
	/** Returns a reference to the vector of upper bounds on the variables #x#.
	  *
	  * Preconditions: \begin{itemize}
	  * \item this->has_bounds() == true (throw NoBoundsOnVariables)
	  * \end{itemize}
	  */
	virtual const SpVectorSlice xu() const = 0;
	///
	/** Get the initial value of the Lagrange multipliers lambda.
	  *
	  * By default this function just sets them to zero.
	  *
	  * Postconditions: \begin{itemize}
	  * \item [lambda != NULL] lambda->size() == this->n()
	  * \item [nu != NULL] nu->size() == this->n()
	  * \end{itemize}
	  *
	  * @param lambda  [out] Pointer to lagrange multipliers for equalities.
	  *                lambda == NULL is allowd in which case it will not
	  *                be set.
	  * @param nu      [out] Pointer to lagrange multipliers for bounds.
	  *                nu == NULL is allowd in which case it will not
	  *                be set.
	  */
	virtual void get_init_lagrange_mult( Vector* lambda, SpVector* nu ) const;

	//@}

	/** @name <<std aggr>> stereotype member functions
	  */
	//@{

	/** @name <<std aggr>> members for the objective function value, f.
	  *
	  * The is a reference to a scalar storage location that is updated when
	  * #calc_f(...)# is called or when #f# is updated as a side effect
	  * when another "calc" member funciton is called.
	  */
	//@{

	///
	virtual void set_f(value_type* f);
	///
	virtual value_type* get_f();
	///
	virtual value_type& f();
	///
	virtual const value_type& f() const;

	//@}

	/** @name <<std aggr>> members for the equality constaints vector, c.
	  *
	  * The is a reference to a vector storage location that is updated when
	  * #calc_c(...)# is called or when #c# is updated as a side effect
	  * when another "calc" member funciton is called.
	  */
	//@{

	///
	virtual void set_c(Vector* c);
	///
	virtual Vector* get_c();
	///
	virtual Vector& c();
	///
	virtual const Vector& c() const;

	//@}

	//@}

	/** @name Calculation Members.
	  *
	  * These members coordinate the calculation of quanaties at various points #x#.
	  * The #set_info(info)#, where #info# = #f#, or #c#, are used to set references to
	  * storage to quanities that are updated when #calc_info(x,newx)# memebers are called.
	  * See the <<std comp>> stereotype definition.
	  *
	  * Preconditions: \begin{itemize}
	  * \item #this->set_info()# has been called before this member (throw NoRefSet)
	  * \end{itemize}
	  *
	  * Postconditions: \begin{itemize}
	  * \item #this->info()# returns a const reference to the value to the updated
	  *		quantity at #x#
	  * \end{itemize}
	  *
	  * @param	x		[I] Point at which to calculate the object function #f#.
	  * @param	newx	[I] true: The values of #x# are the same as the last call to a 
	  *						#this->calc_info(x,newx)# member.
	  *						false: The values of #x# where not the same as the last call to a
	  *						#this->calc_info(x,newx)# member.
	  *						Default value = true
	  */
	//@{

	///
	/** Set whether the subclass can update multiple quanities to save recalculations.
	  *
	  * @param	set_mult_calc	[I] true: The subclass is allowed to update multiple quantities if it is
	  *							more efficient to do so.  For example, calling calc_f(...) to update 'f'
	  *							may result in an update of 'c' if it is more efficent to do so.
	  *							false: The subclass is not allowed to update multiple quantities.
	  */
	virtual void set_mult_calc(bool mult_calc) const = 0;
	///
	/** Set the scaling of the objective function.
	  */
	virtual void scale_f( value_type scale_f ) = 0;
	///
	/** Get the scaling being used for the objective function.
	  */
	virtual value_type scale_f() const = 0;
	///
	/** Update the value for #f# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then storage reference for #c# may also be changed
	  * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_f(const VectorSlice& x, bool newx = true) const;
	///
	/** Update the vector for #c# at the point #x# and put it in the stored reference.
	  *
	  * If #set_mult_calc(true)# was called then storage reference for #f# may also be changed
	  * but is not guarentied to be.  But no other quanities from possible subclasses are allowed
	  * to be updated as a side effect.
	  */ 
	virtual void calc_c(const VectorSlice& x, bool newx = true) const;

	//@}

	///
	/** Used by the solver to report the final solution and multipliers.
	  *
	  * Call this function to report the final solution of the
	  * unknows x and the Lagrange multipliers for the
	  * equality constriants #lambda# and the varaible bounds
	  * #nu#.  If either of the multipliers
	  * are not known then you can pass null in for them.
	  *
	  * The default behavior is to just ignore this.
	  */
	virtual void report_final_solution(
		  const VectorSlice&	x
		, const VectorSlice*	lambda
		, const SpVectorSlice*	nu
		, bool					optimal		) const;

	/** @name Objective and constraint function counts.
	  *
	  * These functions can be called to find out how many evaluations
	  * the client requested since initialize() was called.
	  *
	  * These do not include any function evaluations that may have
	  * been used for finite difference evaluations or anything
	  * like that.
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
	virtual size_type num_f_evals() const;

	///
	/** Gives the number of constraint function c(x) evaluations called by the solver
	  * since initialize() was called..
	  */
	virtual size_type num_c_evals() const;

	//@}

protected:

	///
	/** Struct for objective and constriants (pointer)
	 */
	struct ZeroOrderInfo {
	public:
		///
		ZeroOrderInfo() : f(NULL), c(NULL)
		{}
		///
		ZeroOrderInfo( value_type* f_in, Vector* c_in )
			: f(f_in), c(c_in)
		{}
		/// Pointer to objective function #f# (may be NULL if not set)
		value_type*   f;
		/// Pointer to constraints residule #c# (may be NULL if not set)
		Vector*       c;
	}; // end struct ZeroOrderInfo

	/// Return pointer to set quantities
	const ZeroOrderInfo zero_order_info() const;

	/** @name Protected methods to be overridden by subclasses */
	//@{

	///
	/** Overridden to compute f(x) and perhaps c(x) (if multiple calculaiton = true).
	 *
	 * @param x           [in]  Unknown vector (size n).
	 * @param newx        [in]  True if is a new point.
	 * @param zero_order  [out] Pointers to f and c
	 */
	virtual void imp_calc_f(const VectorSlice& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;
	///
	/** Overridden to compute c(x) and perhaps f(x) (if multiple calculaiton = true).
	 *
	 * @param x           [in]  Unknown vector (size n).
	 * @param newx        [in]  True if is a new point.
	 * @param zero_order  [out] Pointers to f and c
	 */
	virtual void imp_calc_c(const VectorSlice& x, bool newx, const ZeroOrderInfo& zero_order_info) const = 0;

	//@}

	/// Assert referece has been set for a quanity
	template<class T>
	void assert_ref_set(T* p, std::string info) const {
		StandardCompositionRelationshipsPack::assert_role_name_set(p, false, info);
	}

private:
	mutable ZeroOrderInfo           f_c_;
	mutable size_type				num_f_evals_;
	mutable size_type				num_c_evals_;
	
};	// end class NLP

// /////////////////
// Inline members

inline
const NLP::ZeroOrderInfo NLP::zero_order_info() const
{
	return f_c_;
}

}	// end namespace NLPInterfacePack 

#endif // NLP_H
