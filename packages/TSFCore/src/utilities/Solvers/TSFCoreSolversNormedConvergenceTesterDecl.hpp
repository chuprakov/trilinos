// ///////////////////////////////////////////////////////////////
// NormedConvergenceTesterDecl.hpp

#ifndef TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
#define TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP

#include <vector>

#include "TSFCoreSolversConvergenceTester.hpp"

namespace TSFCore {
namespace Solvers {

///
/** Convergence test based on relative errors.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class NormedConvergenceTester : public ConvergenceTester<Scalar> {
public:

	/** @name Constructors / initializers */
	//@{
	
	/// Calls <tt>initialize()</tt>
	NormedConvergenceTester(
		const Scalar                                           tol             = 1e-12
		,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm           = MemMngPack::null
		);

	/// Calls <tt>initialize()</tt>
	NormedConvergenceTester(
		const Index                                            totalNumSystems
		,const Scalar                                          tols[]
		,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm           = MemMngPack::null
		);

	///
	/** Initialize with one set of tolerances for multiple linear systems
	 */
	void initialize(
		const Scalar                                           tol             = 1e-12
		,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm           = MemMngPack::null
		);

	///
	/** Initialize with different tolerances for multiple linear systems.
	 */
	void initialize(
		const Index                                            totalNumSystems
		,const Scalar                                          tols[]
		,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm           = MemMngPack::null
		);

	//@}

	/** @name Overridden from ConvergenceTester */
	//@{

	///
	MemMngPack::ref_count_ptr<const Norm<Scalar> > norm() const;
	///
	void convStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		);

	//@}

private:

	MemMngPack::ref_count_ptr<const Norm<Scalar> >   norm_;
	std::vector<Scalar>                              tols_;           // if size()==1 then works for multiple systems
	std::vector<Index>                               activeSystems_;  // cache
	std::vector<Scalar>                              norms_;          // cache

}; // class NormedConvergenceTester

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
