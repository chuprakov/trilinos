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
// TSFCoreSolversSummaryOutputter.hpp

#ifndef TSFCORE_SOLVERS_SUMMARY_OUTPUTTER_HPP
#define TSFCORE_SOLVERS_SUMMARY_OUTPUTTER_HPP

#include "TSFCoreSolversSummaryOutputterDecl.hpp"

namespace TSFCore {
namespace Solvers {

// Constructors / initializers

template<class Scalar>
SummaryOutputter<Scalar>::SummaryOutputter()
	:resetCalled_(true)
{}

template<class Scalar>
SummaryOutputter<Scalar>::SummaryOutputter(
	const Teuchos::RefCountPtr<std::ostream>       &out
	,const std::string                             &leadingOutputStr
	)
	:out_(out)
	,leadingOutputStr_(leadingOutputStr)
	,resetCalled_(true)
{
	TEST_FOR_EXCEPT(out.get()==NULL);
}

// Overridden from ConvergenceTester

template<class Scalar>
void SummaryOutputter<Scalar>::reset()
{
	Teuchos::RefCountPtr<ConvergenceTester<Scalar> > attachedConvTester = getAttachedConvTester();
	if(attachedConvTester.get()) attachedConvTester->reset();
	resetCalled_ = true;
}

template<class Scalar>
void SummaryOutputter<Scalar>::convStatus(
	const SolverState<Scalar>     &solver
	,const Index                  currNumSystems
	,bool                         isConverged[]
	)
{
	Teuchos::RefCountPtr<ConvergenceTester<Scalar> > attachedConvTester = getAttachedConvTester();
	if(attachedConvTester.get()) attachedConvTester->convStatus(solver,currNumSystems,isConverged);
	else std::fill_n(isConverged,currNumSystems,false);
	if(out_.get()) {
		const std::string &leadstr = leadingOutputStr();
		if(resetCalled_) {
			*out_ << leadstr << "Linear solve using type \'"<<typeid(solver).name()<<"\':\n"; 
			resetCalled_ = false;
		}
		if( static_cast<int>(norms_.size()) < currNumSystems ) norms_.resize(currNumSystems);
		solver.currEstRelResidualNorms(&norms_[0]);
		const Scalar maxNorm = *std::max_element(&norms_[0],&norms_[0]+currNumSystems);
		*out_ << leadstr << "iter="<<solver.currIteration()<<", max{||R||}="<<maxNorm<<std::endl;
	}
	// ToDo: Print the status
}

// Overridden from AttachedConvergenceTesterBase

template<class Scalar>
void SummaryOutputter<Scalar>::protectedReset()
{
	TEST_FOR_EXCEPT(true); // Never called but just in case ...
}

template<class Scalar>
void SummaryOutputter<Scalar>::protectedConvStatus(
	const SolverState<Scalar>     &solver
	,const Index                  currNumSystems
	,bool                         isConverged[]
	)
{
	TEST_FOR_EXCEPT(true); // Never called but just in case ...
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_SUMMARY_OUTPUTTER_HPP
