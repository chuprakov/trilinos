// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ////////////////////////////////////////////////////////////////////////
// TSFCoreNonlinNonlinearProblemFirstOrderDecl.hpp

#ifndef TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_DECL_HPP
#define TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_DECL_HPP

#include "TSFCoreNonlinNonlinearProblem.hpp"
#include "Teuchos_AbstractFactory.hpp"

namespace TSFCore {
namespace Nonlin {

///
/** Specializes the <tt>NonlinearProblem</tt> interface by adding
 * first derivative quantities (i.e.~Jacobians).
 *
 * Here the Jacobian objects are represented as:
 * <ul>
 * <li><tt>D(c)/D(y) = op(DcDy)</tt> (where <tt>DcDy</tt> is a <tt>LinearOpWithSolve</tt> object and
 *      <tt>op(DcDy) == DcDy</tt> if <tt>opDcDy()==NOTRANS</tt>,
 *      and <tt>op(DcDy) == DcDy'</tt> if <tt>opDcDy()==TRANS</tt> 
 * <li><tt>D(c)/D(u(l)) = op(DcDu(l))</tt> (where <tt>DcDu(l)</tt> is a <tt>LinearOp</tt> object and
 *      <tt>op(DcDu(l)) == DcDu(l)'</tt> if <tt>opDcDu(l)()==NOTRANS</tt>,
 *      and <tt>op(DcDu(l)) == DcDu(l)</tt> if <tt>opDcDu(l)()==TRANS</tt>
 * <li><tt>D(g)/D(y) = DgDy'</tt> where <tt>DgDy</tt> is a <tt>MultiVector</tt> object.
 * <li><tt>D(g)/D(u(l)) = DgDu(l)'</tt> where <tt>DgDu(l)</tt> is a <tt>MultiVector</tt> object.
 * </ul>
 *
 * In the above notation <tt>D(f)/D(x)</tt> is a Jacobian matrix of
 * dimmension <tt>dim(f) x dim(x)</tt> where <tt>(D(f)/D(x))(i,j)</tt>
 * is the derivative of the <tt>ith</tt> function <tt>f(i)</tt> with
 * respect to the <tt>jth</tt> variable <tt>x(j)</tt>.
 *
 * <b>Warning!</b> Note that the <tt>MultiVector</tt> objects that are
 * computed for <tt>DcDy</tt> and <tt>DgDu</tt> are actually
 * transposes of the Jacobians <tt>D(g)/D(y)</tt> and
 * <tt>D(g)/D(u)</tt>.  This is a little confusing but make sense in
 * the context of things.
 * 
 * The reason fo the methods <tt>opDcDy()</tt> and <tt>opDcDu()</tt>
 * is to allow the subclass to represent <tt>DcDy</tt> and/or
 * <tt>DcDu(l)</tt> as the transpose or not.  This relieves the need
 * to define adapter subclasses which would destroy the ability
 * to perform dynamic casting.
 *
 * <b>Multiple simultaneous calculations</b>
 *
 * This interface is designed to allow for multiple side-effect
 * calculations so that when one quantity is computed that other
 * quantities may also be computed.  However, the calculation of
 * zero-order quantities (i.e. using <tt>calc_c()</tt> and
 * <tt>calc_g()</tt>) may not result in the side effect calculation of
 * any first-order quantities.  However, the caluation of any
 * first-order quantity may result in the calculation of any other
 * first-order or zero-order quantity (e.g. a call to
 * <tt>calc_DcDy()</tt> may result in calculations of <tt>c</tt>,
 * <tt>g</tt>, <tt>DgDy</tt>, <tt>DcDu(l)</tt> and <tt>DgDu(l)</tt>,
 * for <tt>l=1...Nu()</tt>).  The reason that this is important is so
 * that the subclass can take advantage of opportunities to share
 * computations an other overhead when computing function values and
 * derivatives when possible.  For example, when automatic
 * differentiation is used to compute first-order derivatives, the
 * zero-oerder function values are automatically computed at the same
 * time.
 *
 * <b>Notes to subclass developers</b>
 *
 * If adjoints (i.e. transposed) operations can not be supported for
 * both <tt>DcDy</tt> and <tt>DcDu_l</tt>, then the method
 * <tt>adjointSupported()</tt> should be overridden to return
 * <tt>false</tt>.
 *
 * In order to implement a concrete subclass for a simple nonlinear
 * problem with no auxiilary variables or response functions, only
 * the following methods must be overridden (in addition to the
 * required methods in <tt>NonlinearProblem</tt>): <tt>factory_DcDy()</tt>,
 * <tt>set_DcDy()</tt>, non-<tt>const</tt> <tt>get_DcDy()</tt>
 * and <tt>calc_DcDy()<tt>.
 *
 * In order for a subclass to support auxiliary variables and response
 * functions, the following methods must also be overridden:
 * <tt>factory_DcDu()</tt>, <tt>set_DcDu()</tt>, non-<tt>const</tt>
 * <tt>get_DcDu()</tt>, <tt>set_DgDy()</tt>, non-<tt>const</tt>
 * <tt>get_DgDy()</tt>, <tt>set_DgDu()</tt>, non-<tt>const</tt>
 * <tt>get_DgDu()</tt>, <tt>calc_DcDu()</tt>, <tt>calc_DgDy()</tt> and
 * <tt>calc_DgDu()</tt>.
 *
 * The <tt>const</tt> methods <tt>get_DcDy(), <tt>get_DcDu(),
 * <tt>get_DgDy() and <tt>get_DgDu() never need to be overridden.
 */
template<class Scalar>
class NonlinearProblemFirstOrder : virtual public NonlinearProblem<Scalar> {
public:

	/** @name Adjoints supported? */
	//@{

	///
	/** Return if adjoints are supported on the Jacobian objects
	 * <tt>DcDy</tt> and <tt>DcDu</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The default implementation returns <tt>true</tt>.
	 */
	virtual bool adjointSupported() const;

	//@}

	/** @name Factories for linear operators */
	//@{

	///
	/** Return a factory object for creating objects for <tt>DcDy</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>return.get()!=NULL</tt>
	 * </ul>
	 */
	virtual Teuchos::RefCountPtr< const Teuchos::AbstractFactory< LinearOpWithSolve<Scalar> > > factory_DcDy() const = 0;

	///
	/** Return a factory object for creating objects for <tt>DcDu(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
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
	virtual Teuchos::RefCountPtr< const Teuchos::AbstractFactory<LinearOp<Scalar > > > factory_DcDu(int l) const;

	//@}

	/** @name Transpose arguments */
	//@{

	///
	/** Return the transpose argument for <tt>DcDy</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual ETransp opDcDy() const = 0;

	///
	/** Return th transpose argument for <tt>DcDu(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual ETransp opDcDu(int l) const;

	//@}

	/** @name Set and access calculation storage */
	//@{

	///
	/** Set a pointer to a <tt>LinearOpInvertiable</tt> object to be
	 * updated when <tt>this->calc_DcDy()</tt> is called.
	 *
	 * @param  DcDy  [in] Pointer to state constraint jacobian <tt>D(c)/D(y)</tt>.
	 *                May be <tt>NULL</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>DcDy!=NULL</tt>] Blah blah blah ...
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*this->get_DcDy()==DcDy</tt>
	 * </ul>
	 */
	virtual void set_DcDy(LinearOpWithSolve<Scalar>* DcDy) = 0;
	///
	/** Return the non-const pointer passed to <tt>this->set_DcDy()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 */
	virtual LinearOpWithSolve<Scalar>* get_DcDy() = 0;
	///
	/** Return the const pointer passed to <tt>this->set_DcDy()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The function has a default implementation that returns
	 * <tt>const_DcDyast<NonlinearProblemFirstOrder*>(this)->get_DcDy()</tt>.
	 */
	virtual const LinearOpWithSolve<Scalar>* get_DcDy() const;

	///
	/** Set a pointer to a <tt>LinearOp</tt> object to be
	 * updated when <tt>this->calc_DcDu()</tt> is called.
	 *
	 * @param  l       [in] Index to select the set of auxiliary variables <tt>u(l)</tt>.
	 * @param  DcDu_l  [in] Pointer to state constraint jacobian <tt>D(c)/D(y)</tt>.
	 *                 May be <tt>NULL</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= l <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>DcDu_l!=NULL</tt>] Blah blah blah ...
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_DcDu(l)==DcDu_l</tt>
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual void set_DcDu(int l, LinearOp<Scalar>* DcDu_l);
	///
	/** Return the non-const pointer passed to <tt>this->set_DcDu(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */
	virtual LinearOp<Scalar>* get_DcDu(int l);
	///
	/** Return the const pointer passed to <tt>this->set_DcDu(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * The function has a default implementation that returns
	 * <tt>const_cast<NonlinearProblemFirstOrder*>(this)->get_DcDu(l)</tt>.
	 */
	virtual const LinearOp<Scalar>* get_DcDu(int l) const;

	///
	/** Set a pointer to a <tt>MultiVector</tt> object to be
	 * updated when <tt>this->calc_DgDy()</tt> is called.
	 *
	 * @param  DgDy  [in] Pointer to Jacobian <tt>D(g)/D(y)</tt>.
	 *                May be <tt>NULL</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li> [<tt>DgDy!=NULL</tt>] Blah blah blah ...
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_DgDy()==DgDy</tt>
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->numResponseFunctions()==0</tt>.
	 */
	virtual void set_DgDy(MultiVector<Scalar>* DgDy);
	///
	/** Return the non-const pointer passed to <tt>this->set_DgDy()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->numResponseFunctions()==0</tt>.
	 */
	virtual MultiVector<Scalar>* get_DgDy();
	///
	/** Return the const pointer passed to <tt>this->set_DgDy()</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * </ul>
	 *
	 * The function has a default implementation that returns
	 * <tt>const_DgDyast<NonlinearProblemFirstOrder*>(this)->get_DgDy()</tt>.
	 */
	virtual const MultiVector<Scalar>* get_DgDy() const;

	///
	/** Set a pointer to a <tt>LinearOp</tt> object to be
	 * updated when <tt>this->calc_DgDu()</tt> is called.
	 *
	 * @param  l      [in] Index to select the set of auxiliary variables <tt>u(l)</tt>.
	 * @param  DgDu_l [in] Pointer to state constraint jacobian <tt>D(c)/D(y)</tt>.
	 *                May be <tt>NULL</tt>.
	 *
	 * Preconditions:<ul>
	 * <li><tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= l <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> [<tt>DgDu_l!=NULL</tt>] Blah blah blah ...
	 *      (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->get_DgDu(l)==DgDu_l</tt>
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu() == 0</tt> and <tt>this->numResponseFunctions()==0</tt>.
	 */
	virtual void set_DgDu(int l, MultiVector<Scalar>* DgDu_l);
	///
	/** Return the non-const pointer passed to <tt>this->set_DgDu(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu() == 0</tt> and <tt>this->numResponseFunctions()==0</tt>.
	 */
	virtual MultiVector<Scalar>* get_DgDu(int l);
	///
	/** Return the const pointer passed to <tt>this->set_DgDu(l)</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>this->isInitialized() == true</tt> (throw <tt>NotInitialized</tt>)
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * The function has a default implementation that returns
	 * <tt>const_cast<NonlinearProblemFirstOrder*>(this)->get_DgDy()</tt>.
	 */
	virtual const MultiVector<Scalar>* get_DgDu(int l) const;

	//@}

	/** @name Calculation methods */
	//@{

	///
	/** Update the Jacobian <tt>D(c)/D(y)</tt> at the point <tt>y,{u(l)}</tt>
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
	 * <li><tt>this->get_DcDy()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li>[<tt>this->Nu()==0</tt>] <tt>u==NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*this->get_DcDy()</tt> is updated to <tt>D(c)/D(y)</tt>.
	 * <li> [<tt>this->adjointSupported()</tt>] <tt>this->get_DcDy()->opSupported(trans_trans(this->opDcDy(),TRANS))==true</tt>
	 * </ul>
	 *
	 * The set storage references for the other Jacobian matrices and
	 * the vector functions <tt>c</tt> and/or <tt>g</tt> may also be
	 * updated at this point but is not guaranteed to be.
	 */ 
	virtual void calc_DcDy(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]      = NULL
		,bool                    newPoint = true
		) const = 0;

	///
	/** Update the Jacobian <tt>D(c)/D(u(l))</tt> at the point <tt>y,{u(l)}</tt>
	 * and put it in the stored reference.
	 *
	 * @param  l     [in] Specifies the <tt>lth</tt> Jacobian <tt>D(c)/D(u(l))</tt> to compute.
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
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li><tt>y.space()->is_compatible(*this->space_y())==true</tt>
	 *     (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li><tt>this->get_DcDu(l)()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li>[<tt>this->Nu()==0</tt>] <tt>u==NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*this->get_DcDu(l)</tt> is updated to <tt>D(c)/D(u(l))</tt>.
	 * <li> [<tt>this->adjointSupported()</tt>] <tt>this->get_DcDu()->opSupported(trans_trans(this->opDcDu(),TRANS))==true</tt>
	 * </ul>
	 *
	 * The set storage references for the other Jacobian matrices and
	 * the vector functions <tt>c</tt> and/or <tt>g</tt> may also be
	 * updated at this point but is not guaranteed to be.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu()==0</tt>.
	 */ 
	virtual void calc_DcDu(
		int                      l
		,const Vector<Scalar>    &y
		,const Vector<Scalar>*   u[]      = NULL
		,bool                    newPoint = true
		) const;


	///
	/** Update the Jacobian <tt>D(g)/D(y)</tt> at the point <tt>y,{u(l)}</tt>
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
	 * <li><tt>this->get_DgDy()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li>[<tt>this->Nu()==0</tt>] <tt>u==NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*this->get_DgDy()</tt> is updated to <tt>D(g)/D(y)'</tt>.
	 * </ul>
	 *
	 * The set storage references for the other Jacobian matrices and
	 * the vector functions <tt>c</tt> and/or <tt>g</tt> may also be
	 * updated at this point but is not guaranteed to be.
	 *
	 * <b>Warning!</b> Note that <tt>*this->get_DgDy()</tt> is
	 * actually the transpose of <tt>D(g)/D(y)</tt> as stated in the
	 * class overview and the above postconditions.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->numResponseFunctions()==0</tt>.
	 */ 
	virtual void calc_DgDy(
		const Vector<Scalar>     &y
		,const Vector<Scalar>*   u[]      = NULL
		,bool                    newPoint = true
		) const;

	///
	/** Update the Jacobian <tt>D(g)/D(u(l))</tt> at the point <tt>y,{u(l)}</tt>
	 * and put it in the stored reference.
	 *
	 * @param  l     [in] Specifies the <tt>lth</tt> Jacobian <tt>D(g)/D(u(l))</tt> to compute.
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
	 * <li><tt>this->Nu() > 0</tt> (throw <tt>std::logic_error</tt>)
	 * <li><tt>1 <= u <= this->Nu()</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li><tt>y.space()->is_compatible(*this->space_y())==true</tt>
	 *     (throw <tt>Exceptions::IncompatibleVectorSpaces</tt>)
	 * <li><tt>this->get_DgDu(l)()!=NULL</tt> (throw <tt>std::logic_error</tt>)
	 * <li>[<tt>this->Nu()==0</tt>] <tt>u==NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>*this->get_DgDu(l)</tt> is updated to <tt>D(g)/D(u(l))'</tt>.
	 * </ul>
	 *
	 * The set storage references for the other Jacobian matrices and
	 * the vector functions <tt>c</tt> and/or <tt>g</tt> may also be
	 * updated at this point but is not guaranteed to be.
	 *
	 * <b>Warning!</b> Note that <tt>*this->get_DgDu(l)</tt> is
	 * actually the transpose of <tt>D(g)/D(u(l))</tt> as stated in
	 * the class overview and the above postconditions.
 	 *
	 * The default implementation throws an exception since by default
	 * <tt>this->Nu() == 0</tt> and <tt>this->numResponseFunctions()==0</tt>.
	 */
	virtual void calc_DgDu(
		int                      l
		,const Vector<Scalar>    &y
		,const Vector<Scalar>*   u[]      = NULL
		,bool                    newPoint = true
		) const;

	//@}

}; // class NonlinearProblemFirstOrder

} // namespace Nonlin
} // namespace TSFCore

#endif // TSFCORE_NONLIN_NONLINEAR_PROBLEM_FIRST_ORDER_DECL_HPP
