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
#include "TSFCoreSolversNorm.hpp"
#include "TSFCoreSolversSolverState.hpp"
#include "Teuchos_TestForException.hpp"

namespace TSFCore {
namespace Solvers{

// Constructors / initializers

template<class Scalar>
NormedConvergenceTester<Scalar>::NormedConvergenceTester(
	const Scalar                                      tol
	,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm
	)
{
	initialize(tol,norm);
}

template<class Scalar>
NormedConvergenceTester<Scalar>::NormedConvergenceTester(
	const Index                                       totalNumSystems
	,const Scalar                                     tols[]
	,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm
	)
{
	initialize(totalNumSystems,tols,norm);
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::initialize(
	const Scalar                                      tol
	,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm
	)
{
	tols_.resize(1);
	tols_[0] = tol;
	if(norm.get())   norm_ = norm;
	else             norm_ = Teuchos::rcp(new Norm<Scalar>());
	minTol_    = tol;
	minMaxErr_ = 1e+50;
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::initialize(
	const Index                                       totalNumSystems
	,const Scalar                                     tols[]
	,const Teuchos::RefCountPtr<const Norm<Scalar> >  &norm
	)
{
	tols_.resize(totalNumSystems);
	std::copy( tols, tols + totalNumSystems, tols_.begin() );
	if(norm.get())   norm_ = norm;
	else             norm_ = Teuchos::rcp(new Norm<Scalar>());
	minTol_    = *std::min_element( tols, tols + totalNumSystems );
	minMaxErr_ = 1e+50;
}

// Overridden from ConvergenceTester

template<class Scalar>
Teuchos::RefCountPtr<const Norm<Scalar> >
NormedConvergenceTester<Scalar>::norm() const
{
	return norm_;
}

template<class Scalar>
void NormedConvergenceTester<Scalar>::reset()
{
	return; // There is nothing to reset!
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
	TEST_FOR_EXCEPTION(
		tols_size == 0, std::logic_error
		,"NormedConvergenceTester<Scalar>::convStatus(...): Error!"
		);
	TEST_FOR_EXCEPTION(
		tols_size > 1 && totalNumSystems != tols_size, std::logic_error
		,"NormedConvergenceTester<Scalar>::convStatus(...): Error!"
		);
	TEST_FOR_EXCEPTION(
		currNumSystems > totalNumSystems, std::logic_error
		,"NormedConvergenceTester<Scalar>::convStatus(...): Error!"
		);
#endif
	if(static_cast<Index>(activeSystems_.size()) != totalNumSystems) {
		activeSystems_.resize(totalNumSystems);
		norms_.resize(totalNumSystems);
	}
	solver.currActiveSystems(&activeSystems_[0]);
	solver.currEstRelResidualNorms(&norms_[0]);
	Scalar maxErr = 0.0;
	for(int k = 0; k < currNumSystems; ++k ) {
		isConverged[k] = ( norms_[k] <= ( tols_size > 1 ? tols_[activeSystems_[k]] : tols_[0] ) );
		if( maxErr < norms_[k] ) maxErr = norms_[k]; 
	}
	if( minMaxErr_ > maxErr ) minMaxErr_ = maxErr; 
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_NORMED_CONVERGENCE_TESTER_HPP
