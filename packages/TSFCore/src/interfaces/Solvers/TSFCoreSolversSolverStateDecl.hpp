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
