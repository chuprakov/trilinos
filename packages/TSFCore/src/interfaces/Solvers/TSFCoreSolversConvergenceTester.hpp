// ///////////////////////////////////////////////////////////////
// TSFCoreSolversConvergenceTester.hpp

#ifndef TSFCORE_SOLVERS_CONVERGENCE_TESTER_HPP
#define TSFCORE_SOLVERS_CONVERGENCE_TESTER_HPP

#include "TSFCoreSolversConvergenceTesterDecl.hpp"

namespace TSFCore {
namespace Solvers{

template<class Scalar>
MemMngPack::ref_count_ptr<const Norm<Scalar> >
ConvergenceTester<Scalar>::norm() const
{
	return MemMngPack::rcp( new Norm<Scalar>() );
}

template<class Scalar>
MemMngPack::ref_count_ptr<ConvergenceTester<Scalar> >
ConvergenceTester<Scalar>::clone()
{
	return MemMngPack::null;
}


} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_CONVERGENCE_TESTER_HPP
