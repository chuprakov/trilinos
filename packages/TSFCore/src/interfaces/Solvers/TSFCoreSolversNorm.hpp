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
// TSFCoreSolversNorm.hpp

#ifndef TSFCORE_SOLVERS_NORM_HPP
#define TSFCORE_SOLVERS_NORM_HPP

#include "TSFCoreSolversNormDecl.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "TSFCoreVector.hpp"
#include "TSFCoreMultiVectorCols.hpp"
#include "TSFCoreMultiVectorStdOps.hpp"

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
	TSFCore::norms(X,norms);
}

} // namespace Solvers
} // namespace TSFCore

#endif // TSFCORE_SOLVERS_NORM_HPP
