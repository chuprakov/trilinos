// ///////////////////////////////////////////////////////////////
// NormedConvergenceTesterDecl.hpp

#ifndef TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
#define TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP

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
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	/// Calls <tt>initialize()</tt>
	NormedConvergenceTester(
		const Index                                            totalNumSystems
		,const Scalar                                          tols[]
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	///
	/** Initialize with one set of tolerances for multiple linear systems
	 */
	void initialize(
		const Scalar                                           tol             = 1e-12
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	///
	/** Initialize with different tolerances for multiple linear systems.
	 */
	void initialize(
		const Index                                            totalNumSystems
		,const Scalar                                          tols[]
		,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm           = Teuchos::null
		);

	/// Return the minimum (over all iterations) maximum (over all right-hand-sides) tolerance seen
	Scalar minMaxError() const;

	//@}

	/** @name Overridden from ConvergenceTester */
	//@{

	///
	Teuchos::RefCountPtr<const Norm<Scalar> > norm() const;
	///
	void reset();
	///
	void convStatus(
		const SolverState<Scalar>     &solver
		,const Index                  currNumSystems
		,bool                         isConverged[]
		);

	//@}

private:

	Teuchos::RefCountPtr<const Norm<Scalar> >   norm_;
	std::vector<Scalar>                              tols_;           // if size()==1 then works for multiple systems
	std::vector<Index>                               activeSystems_;  // cache
	std::vector<Scalar>                              norms_;          // cache

	Scalar                                           minMaxErr_;

}; // class NormedConvergenceTester

// ///////////////////////////
// Inline members

template<class Scalar>
inline
Scalar NormedConvergenceTester<Scalar>::minMaxError() const
{
	return minMaxErr_;
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_DECL_HPP
