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
// TSFCoreSolversSolverStateDecl.hpp

#ifndef TSFCORE_SOLVERS_SOLVER_STATE_DECL_HPP
#define TSFCORE_SOLVERS_SOLVER_STATE_DECL_HPP

#include "TSFCoreSolversTypes.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Abstract interface for equation solvers that is designed for use by <tt>ConvergenceTester</tt>.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class SolverState {
public:

	///
	virtual ~SolverState() {}

	/** @name Pure virtual functions that must be overridden */
	//@{

	///
	/** Return the total number of systems of equations being solved.
	 */
	virtual Index totalNumSystems() const = 0;

	///
	/** Return the current number of systems actively being solved.
	 *
	 * Postconditions:<ul>
	 * <li><tt>return <= this->totalNumSystems()</tt>
	 * </ul>
	 */
	virtual Index currNumSystems() const = 0;

	///
	/** Return the current iteration count (1-based).
	 */
	virtual int currIteration() const = 0;

	///
	/** Return the indexes of the current active systems currently being solved.
	 *
	 * @param  activeSystems  [out] Array (length <tt>this->currNumSystems()</tt>) of the indexes
	 *                        of the current active systems being solved and are not yet converged.
	 *                        The values <tt>activeSystems[k]</tt>, for <tt>k=1..this->currNumSystems()</tt>
	 *                        determine the indexes <tt>j</tt> of the currently active systems.
	 *
	 * Postconditions:<ul>
	 * <li> The array <tt>activeSystems[]</tt> must be sorted in assending order (i.e.
	 *      <tt>activeSystems[k] < activeSystems[k+1]</tt>, for <tt>k = 0...this->currNumSystemsJ()-2</tt>).
	 * </ul>
	 */
	virtual void currActiveSystems( Index activeSystems[] ) const = 0;

	///
	/** Return norms of the estimates of the relative residual errors in the currently active systems.
	 *
	 * @param  norms  [out] Array (length <tt>this->currNumSystems()</tt>) of the norms of the
	 *                estimates of the relative residual errors of the current active systems being
	 *                solved that are not yet converged.  The values <tt>norm[k]</tt>, for <tt>k=1..this->currNumSystems()</tt>
	 *                give the norms of the extimates of the current active systems specified by
	 *                <tt>activeSystems[k]</tt>, for <tt>k=1..this->currNumSystems()</tt> where <tt>activeSystems[]</tt>
	 *                is returned from <tt>this->currActiveSystems()</tt>.
	 */
	virtual void currEstRelResidualNorms( Scalar norms[] ) const = 0;

	//@}

}; // class SolverState

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_SOLVER_STATE_DECL_HPP
