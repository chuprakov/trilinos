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

// ///////////////////////////////////////////////////////////////
// TSFCoreSolversConvergenceTesterDecl.hpp

#ifndef TSFCORE_SOLVERS_CONVERGENCE_TESTER_DECL_HPP
#define TSFCORE_SOLVERS_CONVERGENCE_TESTER_DECL_HPP

#include "TSFCoreSolversTypes.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Abstract interface for convergence tests for systems of equations.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class ConvergenceTester {
public:

  ///
  virtual ~ConvergenceTester() {}

	/** @name Pure virtual functions that must be overridden */
	//@{

	///
	/** Reset the convergence tester for a new set of iterations for solving
	 * the system of equations.
	 *
	 * Calling this function before each set of iterative solvers on a
	 * block of right-hand sides allows <tt>*this</tt> convergence
	 * tester object to keep track of what is going on.
	 */
	virtual void reset() = 0;

	///
	/** Determing the convergence status for the currently active linear systems.
	 *
	 * @param  solver  [in] The state of the solver object.
	 * @param  currNumSystems
	 *                 [in] The current number of systems actively being solved.
	 * @param  isConverged
	 *                 [out] Array (length <tt>currNumSystems</tt>) that on ouput states
	 *                 which linear systems are converged.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>currNumSystems == solver.SolverState::currNumSystems()</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>isConverged[k]==true</tt> if the current linear system <tt>k</tt> is converged.
	 * </ul>
	 *
	 * This method is declared non-constant since it is assumed that the
	 * internal implementation will record these calls and will
	 * therefore modify the state of the object.
	 */
	virtual void convStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		) = 0;

	///
	/** Attach another convergence tester object.
	 *
	 * @param  convTester  [in] Smart pointer to another convergence tester that
	 *                     can be considered in some way.  It is allowed for
	 *                     <tt>convTester.get()==NULL</tt> in which case
	 *                     the current convergence tester (if one is currently
	 *                     attached) will be unattached.
	 *
	 * Postconditions:<ul>
	 * <li><tt>this->getAttachedConvTester().get() == convTester.get()</tt>
	 * </ul>
	 *
	 * The convergence tester <tt>convTester</tt> being attached can be
	 * dealt any way that <tt>*this</tt> chooses.
	 */
	virtual void attach( const Teuchos::RefCountPtr<ConvergenceTester<Scalar> > &convTester ) = 0;

	///
	/** Get a smart pointer to non-<tt>const</tt> attached convergence
	 * tester.
	 */
	virtual Teuchos::RefCountPtr<ConvergenceTester<Scalar> > getAttachedConvTester() = 0;

	//@}

	/** @name Virtual functions with default implementation */
	//@{

	///
	/** Get a smart pointer to <tt>const</tt> attached convergence
	 * tester.
	 *
	 * The default implementation returns
	 \code

	   const_cast<ConvergenceTester<Scalar>*>(this)->getAttachedConvTester()

   \endcode
	 *
	 * No override of this function should be needed.
	 */
	virtual Teuchos::RefCountPtr<const ConvergenceTester<Scalar> > getAttachedConvTester() const;

	///
	/** Clone the convergence test if supported.
	 *
	 * Default returns <tt>return.get()==NULL</tt>
	 */
	virtual Teuchos::RefCountPtr<ConvergenceTester<Scalar> > clone();

	//@}

}; // class ConvergenceTester

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_CONVERGENCE_TESTER_DECL_HPP
