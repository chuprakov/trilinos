// ////////////////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblemDecl.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_DECL_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_DECL_HPP

#include "TSFCoreNonlinTypes.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Base interface to nonlinear problem composed out of a set of constraints
 * and response functions.
 *
 * The mathematical form (in Matlab-like notation) for the functions
 * that this object represents are:
 \verbatim

             c(y,{u(l)})  = 0
    gL    <= g(y,{u(l)}) <= gU
    yL    <= y           <= yU
    uL(l) <= u(l)        <= uU(l), for l=1...Nu
 \endverbatim
 * where:<ul>
 * <li><tt>y</tt> is the vector of state variables in the space <tt>space_y</tt>
 * <li><tt>u(l)</tt> is the <tt>lth</tt> subvector of auxiliary variables in
 *     the space <tt>space_u(l)</tt>
 * <li><tt>{u(l)}</tt> is the set of auxiliary subvectors <tt>{u(1),u(2),...,u(Nu)}</tt>
 * <li><tt>c(y,{u(l)})</tt> is the vector function of the state constraints that lies in the
 *      space <tt>space_c</tt>
 * <li><tt>g(y,{u(l)})</tt> is the vector function of auxiliary responses that lies in the
 *      space <tt>space_g</tt>
 * <li><tt>gL</tt> and <tt>gU</tt> are the lower and upper bound vectors on the auxiliary
 *     response functions which lie in the space <tt>space_g</tt>.
 * <li><tt>yL</tt> and <tt>yU</tt> are the lower and upper bound vectors on the state
 *     variables <tt>y</tt> which lie in the space <tt>space_y</tt>.
 * <li><tt>uL(l)</tt> and <tt>uU(l)</tt> are the lower and upper bound vectors on the auxiliary
 *     variables <tt>u(l)</tt> which lie in the space <tt>space_u(l)</tt>.
 * </ul>
 *
 * This interface provides access basic problem definition and the
 * ability to compute the zero-order vector functions
 * <tt>c(y,{u(l)})</tt> and <tt>g(y,{u(l)})</tt>.
 *
 *
 * <b>Multiple simultaneous calculations</b>
 *
 * This interface is designed to allow for multiple side-effect
 * calculations so that when one quantity is computed that other
 * quantities may also be computed.  For example, when
 * <tt>calc_c()</tt> is called to compute <tt>c</tt> this may also
 * simultaneously result in the calcuation of <tt>g</tt> (if
 * <tt>get_g() != NULL</tt>).  This allows for more efficient
 * implementations.
 *
 * <b>Notes to subclass developers</b>
 *
 * By default, <tt>Nu()</tt> returns 0 and <tt>space_g()</tt> returns
 * <tt>return.get()==NULL</tt> and all of the methods associated with
 * auxiliary variables and response functions have default
 * implementations that do nothing (but throw exceptions).
 *
 * Only the following methods must be overridden to create a concrete
 * subclass for a simulation-only nonlinear problem:
 * <tt>initialize()</tt>, <tt>isInitialized()</tt>,
 * <tt>space_y()</tt>, <tt>space_c()</tt>, <tt>yL()</tt>,
 * <tt>yU()</tt>, <tt>y0()</tt>, <tt>set_c()</tt>, non-<tt>const</tt>
 * <tt>get_c()</tt>, <tt>unsetQuantities()</tt> and <tt>calc_c()</tt>.
 *
 * In order to add auxiliary variables and response functions, the
 * following methods must be overridden also: <tt>Nu()</tt>,
 * <tt>space_u()</tt>, <tt>space_g()</tt>, <tt>uL()</tt>, <tt>uU</tt>,
 * <tt>gL()</tt>, <tt>gU</tt>, <tt>u0()</tt>, <tt>set_g()</tt>,
 * non-<tt>const</tt> <tt>get_g()</tt> and <tt>calc_g()</tt>.
 *
 * The method <tt>numResponseFunctions()</tt> and the <tt>const</tt>
 * methods <tt>get_c()</tt> and <tt>get_g()</tt> never need to be
 * overridden.
 */
template<class Scalar>
class NonlinearProblem {
public:
	
	///
	virtual ~NonlinearProblem() {}
	
	/// Value for an infinite bound.
	static Scalar infiniteBound();
	
	/** @name Initialization */
	//@{
	
	///
	/** Initialize object.
	 *
	 * @param  testSetup [in] If <tt>true</tt>, then internal tests
	 *                   will be performed the validate the setup of <tt>*this</tt>.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt>
	 * </ul>
	 *
	 * Note that subclasses must call this function to reset what needs to be
	 * reset in this base object.
	 */
	virtual void initialize( bool testSetup = false ) = 0;
	///
	/** Return if <tt>this</tt> is initialized.
	 */
	virtual bool isInitialized() const = 0;
	
	//@}
	
	/** @name Basic information */
	//@{
	
	///
	/** Number of sets of auxiliary variables.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * If this function returns 0, then there are no auxiliary variables.
	 *
	 * The default implementation returns 0.
	 */
	virtual int Nu() const;
	
	///
	/** Number of auxiliary response functions.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * If this function returns 0, then there are no auxiliary response functions.
	 *
	 * The default implementation returns <tt>this->space_g().get() ?
	 * this->space_g()->dim() : 0</tt>.
	 */
	virtual int numResponseFunctions() const;

	//@}

	/** @name VectorSpaces */
	//@{

	///
	/** VectorSpace for the state variables.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get()!=NULL</tt>
	 * </ul>
	 */
	virtual Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_y() const = 0;
	///
	/** VectorSpace for the auxiliary variables <tt>u(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get()!=NULL</tt>
	 * </ul>
	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_u(int l) const;
	///
	/** VectorSpace for the state constraints.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get()!=NULL</tt>
	 * <li> <tt>return->dim() == this->space_y()->dim()</tt>
	 * </ul>
	 */
	virtual Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_c() const = 0;
	///
	/** VectorSpace for the auxiliary response functions.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * If <tt>return.get() == NULL</tt>, then there are no auxiliary
	 * response functions.
 	 *
	 * The default implementation returns <tt>return.get() == NULL</tt>.
	 */
	virtual Teuchos::RefCountPtr<const VectorSpace<Scalar> >  space_g() const;

	//@}

	/** @name Bounds */
	//@{

	///
	/** Lower bounds on state variables.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_y()) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will give
	 * <tt>get_ele(return,i) == -NonlinearProblem::infiniteBound()</tt>.
	 */
	virtual const Vector<Scalar>& yL() const = 0;
	///
	/** Upper bounds on state variabes.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_y()) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will give
	 * <tt>get_ele(return,i) == +NonlinearProblem::infiniteBound()</tt>.
	 */
	virtual const Vector<Scalar>& yU() const = 0;
	///
	/** Lower bounds on axuliary variables <tt>u(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_u(l)) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will give
	 * <tt>get_ele(return,i) == +NonlinearProblem::infiniteBound()</tt>.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual const Vector<Scalar>& uL(int l) const;
	///
	/** Upper bounds on axuliary variables <tt>u(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_u(l)) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will give
	 * <tt>get_ele(return,i) == +NonlinearProblem::infiniteBound()</tt>.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual const Vector<Scalar>& uU(int l) const;
	///
	/** Lower bounds on response functions.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->space_g()().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_g()) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will give
	 * <tt>get_ele(return,i) == -NonlinearProblem::infiniteBound()</tt>.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->numResponseFunctions()==0</tt>.
	 */
	virtual const Vector<Scalar>& gL() const;
	///
	/** Upper bounds on response functions.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->space_g()().get() != NULL</tt> (throw <tt>std::logic_error</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_g()) == true</tt>
	 * </ul>
	 *
	 * Any bounds that are non-existant will give
	 * <tt>get_ele(return,i) == +NonlinearProblem::infiniteBound()</tt>.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->numResponseFunctions()==0</tt>.
	 */
	virtual const Vector<Scalar>& gU() const;

	//@}

	/** @name Initial values (guesses) for state and auxiliary variables */
	//@{

	///
	/** Initial values (guess) for state variables
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_y()) == true</tt>
	 * </ul>
	 */
	virtual const Vector<Scalar>& y0() const = 0;
	///
	/** Initial values (guess) for auxiliary variables <tt>u(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.space()->isCompatible(*this->space_u(l)) == true</tt>
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual const Vector<Scalar>& u0(int l) const;
	//@}

	/** @name Set and access calculation storage */
	//@{

	///
	/** Set a pointer to a vector to be updated when <tt>this->calc_c()</tt> is called.
	 *
	 * @param  c  [in] Pointer to constraint residual <tt>c(y,{u(l)})</tt>.
	 *            May be <tt>NULL</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>c!=NULL</tt>] <tt>c->space()->isCompatible(*this->space_c())==true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_c() == c</tt>
	 * </ul>
	 */
	virtual void set_c(Vector<Scalar>* c) = 0;
	///
	/** Return the non-const pointer passed to <tt>this->set_c()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual Vector<Scalar>* get_c() = 0;
	///
	/** Return the const pointer passed to <tt>this->set_c()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The function has a default implementation that returns
	 * <tt>const_cast<NonlinearProblem*>(this)->get_c()</tt>.
	 */
	virtual const Vector<Scalar>* get_c() const;

	///
	/** Set a pointer to a vector to be updated when <tt>this->calc_g()</tt> is called.
	 *
	 * @param  g  [in] Pointer to response functions <tt>g(y,{u(l)})</tt>.
	 *            May be <tt>NULL</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>g!=NULL</tt>] <tt>g->space()->isCompatible(*this->space_g())==true</tt>
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_g() == g</tt>
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual void set_g(Vector<Scalar>* g);
	///
	/** Return the non-const pointer passed to <tt>this->set_g()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual Vector<Scalar>* get_g();
	///
	/** Return the const pointer passed to <tt>this->set_g()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The function has a default implementation that returns
	 * <tt>const_cast<NonlinearProblem*>(this)->get_g()</tt>.
	 */
	virtual const Vector<Scalar>* get_g() const;
	///
	/** Call to unset all storage quantities (both in this class and all subclasses).
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->get_c()==NULL</tt>
	 * <li><tt>this->get_g()==NULL</tt>
	 * <li>Unsets quantities that may be set in subclasses as well.
	 * </ul>
	 */
	virtual void unsetQuantities() = 0;

	//@}

	/** @name Calculation methods */
	//@{

	///
	/** Update the constraint residual vector for <tt>c</tt> at the point <tt>y,{u(l)}</tt>
	 * and put it in the stored reference.
	 *
	 * @param  y     [in] The current value of the state variables
	 * @param  u     [in] Array (size <tt>this->Nu()</tt>) of the current values of the
	 *               auxiliary variables.  It is allowed for <tt>u==NULL</tt> in which
	 *               interpreted as <tt>u[l-1] == &u0(l)</tt>, for <tt>l=1...this->Nu()</tt>.
	 * @param  newPoint
	 *               [in] (default <tt>true</tt>) If <tt>false</tt>, the values in <tt>y</tt>
	 *               and <tt>u[]</tt> are assumed to be the same as the last call to
	 *               a <tt>this->calc_*(y,u,newPoint)</tt> member.  If <tt>true</tt>, the
	 *               values in <tt>y</tt> and/or <tt>u[]</tt> are not the same as the
	 *               last call to a <tt>this->calc_*(y,u,newPoint)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized()==true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>y.space()->is_compatible(*this->space_y())==true</tt>
	 *     (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li><tt>this->get_c()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li>[<tt>this->Nu()==0</tt>] <tt>u==NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*this->get_c()</tt> is updated to <tt>c(x,{u(l)})</tt>.
	 * </ul>
	 *
	 * The storage reference for <tt>g</tt> may also be updated at
	 * this point (if <tt>this->get_g()!=NULL</tt>) but is not
	 * guaranteed to be.
	 */ 
	virtual void calc_c(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]      = NULL
		,bool                    newPoint = true
		) const = 0;

	///
	/** Update the response functions vector for <tt>g</tt> at the point <tt>y,{u(l)}</tt>
	 * and put it in the stored reference.
	 *
	 * @param  y     [in] The current value of the state variables
	 * @param  u     [in] Array (size <tt>this->Nu()</tt>) of the current values of the
	 *               auxiliary variables.  It is allowed for <tt>u==NULL</tt> in which
	 *               interpreted as <tt>u[l-1] == &u0(l)</tt>, for <tt>l=1...this->Nu()</tt>.
	 * @param  newPoint
	 *               [in] (default <tt>true</tt>) If <tt>false</tt>, the values in <tt>y</tt>
	 *               and <tt>u[]</tt> are assumed to be the same as the last call to
	 *               a <tt>this->calc_*(y,u,newPoint)</tt> member.  If <tt>true</tt>, the
	 *               values in <tt>y</tt> and/or <tt>u[]</tt> are not the same as the
	 *               last call to a <tt>this->calc_*(y,u,newPoint)</tt> member.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized()==true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>y.space()->is_compatible(*this->space_y())==true</tt>
	 *     (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li><tt>this->get_g()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li>[<tt>this->Nu()==0</tt>] <tt>u==NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*this->get_g()</tt> is updated to <tt>g(x,{u(l)})</tt>.
	 * </ul>
	 *
	 * The storage reference for <tt>c</tt> may also be updated at
	 * this point (if <tt>this->get_g()!=NULL</tt>) but is not
	 * guaranteed to be.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */ 
	virtual void calc_g(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]      = NULL
		,bool                    newPoint = true
		) const;

	//@}

}; // class NonlinearProblem

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_DECL_HPP
