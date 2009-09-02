// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_TestForException.hpp"
#include "FEApp_ProblemFactory.hpp"
#include "FEApp_BrusselatorProblem.hpp"
#include "FEApp_HeatNonlinearSourceProblem.hpp"

FEApp::ProblemFactory::ProblemFactory(
       const Teuchos::RCP<Teuchos::ParameterList>& problemParams_,
       const Teuchos::RCP<ParamLib>& paramLib_) :
  problemParams(problemParams_),
  paramLib(paramLib_)
{
}

Teuchos::RCP<FEApp::AbstractProblem>
FEApp::ProblemFactory::create()
{
  Teuchos::RCP<FEApp::AbstractProblem> strategy;

  std::string& method = problemParams->get("Name", "Brusselator");
  if (method == "Brusselator") {
    strategy = Teuchos::rcp(new FEApp::BrusselatorProblem(problemParams, 
                                                          paramLib));
  }
  else if (method == "Heat Nonlinear Source") {
    strategy = 
      Teuchos::rcp(new FEApp::HeatNonlinearSourceProblem(problemParams,
                                                         paramLib));
  }
  else {
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                       std::endl << 
                       "Error!  Unknown problem " << method << 
                       "!" << std::endl << "Supplied parameter list is " << 
                       std::endl << *problemParams);
  }

  return strategy;
}
