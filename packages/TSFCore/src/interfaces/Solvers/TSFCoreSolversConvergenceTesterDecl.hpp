// ///////////////////////////////////////////////////////////////
// TSFCoreSolversConvergenceTesterDecl.hpp

#ifndef TSFCORE_SOLVERS_CONVERGENCE_TESTER_DECL_HPP
#define TSFCORE_SOLVERS_CONVERGENCE_TESTER_DECL_HPP

#include "TSFCoreSolversTypes.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Abstract interface for convergence tests ofr systems of equations.
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
	virtual MemMngPack::ref_count_ptr<const Norm<Scalar> > norm() const;

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
	virtual MemMngPack::ref_count_ptr<ConvergenceTester<Scalar> > clone();

	//@}

}; // class ConvergenceTester

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_CONVERGENCE_TESTER_DECL_HPP
