// ///////////////////////////////////////////////////////////////
// TSFCoreSolversConvergenceTester.hpp

#ifndef TSFCORE_SOLVERS_CONVERGENCE_TESTER_HPP
#define TSFCORE_SOLVERS_CONVERGENCE_TESTER_HPP

#include "TSFCoreSolversConvergenceTesterDecl.hpp"

namespace TSFCore {
namespace Solvers{

template<class Scalar>
Teuchos::RefCountPtr<const Norm<Scalar> >
ConvergenceTester<Scalar>::norm() const
{
	return Teuchos::rcp( new Norm<Scalar>() );
}

template<class Scalar>
Teuchos::RefCountPtr<ConvergenceTester<Scalar> >
ConvergenceTester<Scalar>::clone()
{
	return Teuchos::null;
}


} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_CONVERGENCE_TESTER_HPP
