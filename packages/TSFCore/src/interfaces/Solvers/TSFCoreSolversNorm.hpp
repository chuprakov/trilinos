// ///////////////////////////////////////////////////////////////
// TSFCoreSolversNorm.hpp

#ifndef TSFCORE_SOLVERS_NORM_HPP
#define TSFCORE_SOLVERS_NORM_HPP

#include "TSFCoreSolversNormDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreMultiVectorCols.hpp"

namespace TSFCore {
namespace Solvers {

template<class Scalar>
Scalar Norm<Scalar>::norm(const Vector<Scalar>& x) const
{
	const MultiVectorCols<Scalar> X( Teuchos::rcp( const_cast<Vector<Scalar>*>(&x), false ) );
	Scalar norms[1];
	this->norms(X,norms);
	return norms[0];
}

template<class Scalar>
void Norm<Scalar>::norms( const MultiVector<Scalar>& X, Scalar norms[] ) const
{
	const int m = X.domain()->dim();
	X.range()->scalarProds(X,X,norms);
	for( int j = 0; j < m; ++j )
		norms[j] = sqrt(norms[j]);
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_NORM_HPP
