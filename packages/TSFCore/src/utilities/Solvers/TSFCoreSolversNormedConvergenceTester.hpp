// ///////////////////////////////////////////////////////////////
// TSFCoreSolversNormedConvergenceTester.hpp

#ifndef TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_HPP
#define TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreSolversNormedConvergenceTesterDecl.hpp"
#include "TSFCoreSolversNorm.hpp"
#include "TSFCoreSolversSolverState.hpp"
#include "ThrowException.hpp"

namespace TSFCore {
namespace Solvers{

// Constructors / initializers

template<class Scalar>
NormedConvergenceTester<Scalar>::NormedConvergenceTester(
	const Scalar                                           tol
	,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm
	)
{
	initialize(tol,norm);
}

template<class Scalar>
NormedConvergenceTester<Scalar>::NormedConvergenceTester(
	const Index                                            totalNumSystems
	,const Scalar                                          tols[]
	,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm
	)
{
	initialize(totalNumSystems,tols,norm);
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::initialize(
	const Scalar                                           tol
	,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm
	)
{
	tols_.resize(1);
	tols_[0] = tol;
	if(norm.get())   norm_ = norm;
	else             norm_ = MemMngPack::rcp(new Norm<Scalar>());
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::initialize(
	const Index                                            totalNumSystems
	,const Scalar                                          tols[]
	,const MemMngPack::ref_count_ptr<const Norm<Scalar> >  &norm
	)
{
	tols_.resize(totalNumSystems);
	std::copy( tols, tols + totalNumSystems, tols_.begin() );
	if(norm.get())   norm_ = norm;
	else             norm_ = MemMngPack::rcp(new Norm<Scalar>());
}

// Overridden from ConvergenceTester

template<class Scalar>
MemMngPack::ref_count_ptr<const Norm<Scalar> >
NormedConvergenceTester<Scalar>::norm() const
{
	return norm_;
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::convStatus(
	const SolverState<Scalar>     &solver
	,const Index                  currNumSystems
	,bool                         isConverged[]
	)
{
	const Index  totalNumSystems = solver.totalNumSystems();
	const int    tols_size       = tols_.size();
#ifdef _DEBUG
	THROW_EXCEPTION(
		tols_size == 0, std::logic_error
		,"NormedConvergenceTester<Scalar>::convStatus(...): Error!"
		);
	THROW_EXCEPTION(
		tols_size > 1 && totalNumSystems != tols_size, std::logic_error
		,"NormedConvergenceTester<Scalar>::convStatus(...): Error!"
		);
	THROW_EXCEPTION(
		currNumSystems > totalNumSystems, std::logic_error
		,"NormedConvergenceTester<Scalar>::convStatus(...): Error!"
		);
#endif
	if(activeSystems_.size() != totalNumSystems) {
		activeSystems_.resize(totalNumSystems);
		norms_.resize(totalNumSystems);
	}
	solver.currActiveSystems(&activeSystems_[0]);
	solver.currEstRelResidualNorms(&norms_[0]);
	for(int k = 0; k < currNumSystems; ++k )
		isConverged[k] = ( norms_[k] <= ( tols_size > 1 ? tols_[activeSystems_[k]] : tols_[0] ) );
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_HPP
