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

	/** @name Pure virtual functions that must be overridden */
	//@{

	///
	/** Return the definition of the norm to be used in computing relative errors.
	 *
	 * The default implementation returns <tt>return.get()!=NULL</tt>
	 * and uses the default definition implementation <tt>Norm</tt>.
	 */
	virtual Teuchos::RefCountPtr<const Norm<Scalar> > norm() const;

	///
	/** Reset the convergence tester for a new set of iterations for solving
	 * the system of equations.
	 *
	 * ToDo: Finish documentation!
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
	 * <tt><tt>isConverged[k]==true</tt> if the linear system <tt>j</tt> is converged where
	 *     <tt>j = activeSystems[k]</tt> and where <tt>activeSystems[]</tt> is returned from
	 *     the <tt>solver.SolverState::currActiveSystems()</tt>.
	 * </ul>
	 *
	 * This method is declared conconstant since it is assumed that the internal implementation
	 * will record these calls and will therefore modify the state of the object.
	 */
	virtual void convStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		) = 0;

	///
	virtual Teuchos::RefCountPtr<ConvergenceTester<Scalar> > clone();

	//@}

}; // class ConvergenceTester

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_CONVERGENCE_TESTER_DECL_HPP
