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
// TSFCoreSolversNormedConvergenceTester.hpp

#ifndef TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_HPP
#define TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_HPP

#include "TSFCoreSolversNormedConvergenceTesterDecl.hpp"
#include "TSFCoreSolversSolverState.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Solvers{

// Constructors / initializers

template<class Scalar>
NormedConvergenceTester<Scalar>::NormedConvergenceTester(
	const ScalarMagnitude tol
	,const EAttachmentMode attachmentMode
	)
	:AttachConvergenceTesterBase<Scalar>(attachmentMode)
{
	this->tol(tol);
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::tol(
	const ScalarMagnitude tol
	)
{
	tol_ = tol;
	minMaxErr_ = ScalarMagnitude(1e+50);
}

// Overridden from ConvergenceTester

template<class Scalar>
void NormedConvergenceTester<Scalar>::protectedReset()
{
	return; // There is nothing to reset!
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::protectedConvStatus(
	const SolverState<Scalar>     &solver
	,const Index                  currNumSystems
	,bool                         isConverged[]
	)
{
	if( static_cast<Index>(norms_.size()) < currNumSystems )
		norms_.resize(currNumSystems);
	solver.currEstRelResidualNorms(&norms_[0]);
	ScalarMagnitude maxErr = 0.0;
	for(int k = 0; k < currNumSystems; ++k ) {
		isConverged[k] = ( norms_[k] <= tol_ );
		if( maxErr < norms_[k] ) maxErr = norms_[k]; 
	}
	if( minMaxErr_ > maxErr ) minMaxErr_ = maxErr; 
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_HPP
