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
// TSFCoreSolversAttachConvergenceTesterBase.hpp

#ifndef TSFCORE_SOLVERS_ATTACH_CONVERGENCE_TESTER_BASE_HPP
#define TSFCORE_SOLVERS_ATTACH_CONVERGENCE_TESTER_BASE_HPP

#include "TSFCoreSolversAttachConvergenceTesterBaseDecl.hpp"
#include "Teuchos_Workspace.hpp"

namespace TSFCore {
namespace Solvers {

// Constructors / initializers

template<class Scalar>
AttachConvergenceTesterBase<Scalar>::AttachConvergenceTesterBase(
	const EAttachmentMode attachmentMode
	)
{
	this->attachmentMode(attachmentMode);
}

// Overridden from ConvergenceTester

template<class Scalar>
void AttachConvergenceTesterBase<Scalar>::reset()
{
	protectedReset();
	if( attachedConvTester_.get() && attachmentMode()!=ATTACHED_TEST_EXCLUDE )
		attachedConvTester_->reset();
}

template<class Scalar>
void AttachConvergenceTesterBase<Scalar>::convStatus(
	const SolverState<Scalar>     &solver
	,const Index                  currNumSystems
	,bool                         isConverged[]
	)
{
	using Teuchos::Workspace;
	Teuchos::WorkspaceStore* wss = Teuchos::get_default_workspace_store().get();

	protectedConvStatus(solver,currNumSystems,isConverged);

	if( attachedConvTester_.get() && attachmentMode()!=ATTACHED_TEST_EXCLUDE ) {
		Workspace<bool> isConverged2(wss,currNumSystems,false);
		attachedConvTester_->convStatus(solver,currNumSystems,&isConverged2[0]);
		switch(attachmentMode()) {
			case ATTACHED_TEST_IGNORE: {
				break; // Just ignore what is in isConverged2
			}
			case ATTACHED_TEST_INSTEAD: {
				for( int k = 0; k < currNumSystems; ++k ) isConverged[k] = isConverged2[k];
				break;
			}
			case ATTACHED_TEST_AND: {
				for( int k = 0; k < currNumSystems; ++k ) isConverged[k] = (isConverged[k] && isConverged2[k]);
				break;
			}
			case ATTACHED_TEST_OR: {
				for( int k = 0; k < currNumSystems; ++k ) isConverged[k] = (isConverged[k] || isConverged2[k]);
				break;
			}
			default:
				TEST_FOR_EXCEPT(true);
		}
	}
}

template<class Scalar>
void AttachConvergenceTesterBase<Scalar>::attach(
	const Teuchos::RefCountPtr<ConvergenceTester<Scalar> > &convTester
	)
{
	attachedConvTester_ = convTester;
}

template<class Scalar>
Teuchos::RefCountPtr<ConvergenceTester<Scalar> >
AttachConvergenceTesterBase<Scalar>::getAttachedConvTester()
{
	return attachedConvTester_;
}

} // namespace Solvers
} // namespace TSFCore

#endif  // TSFCORE_SOLVERS_ATTACH_CONVERGENCE_TESTER_BASE_HPP
