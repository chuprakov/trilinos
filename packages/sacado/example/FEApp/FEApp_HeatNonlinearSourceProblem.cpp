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

#include "FEApp_HeatNonlinearSourceProblem.hpp"
#include "FEApp_ConstantNodeBCStrategy.hpp"
#include "FEApp_BoundaryFlux1DResponseFunction.hpp"
#include "FEApp_SolutionAverageResponseFunction.hpp"
#include "FEApp_SolutionTwoNormResponseFunction.hpp"

FEApp::HeatNonlinearSourceProblem::
HeatNonlinearSourceProblem(
	            const Teuchos::RCP<Teuchos::ParameterList>& params_,
              const Teuchos::RCP<ParamLib>& paramLib_) :
  params(params_),
  paramLib(paramLib_)
{
  leftBC = params->get("Left BC", 0.0);
  rightBC = params->get("Right BC", 0.0);
}

FEApp::HeatNonlinearSourceProblem::
~HeatNonlinearSourceProblem()
{
}

unsigned int
FEApp::HeatNonlinearSourceProblem::
numEquations() const
{
  return 1;
}

void
FEApp::HeatNonlinearSourceProblem:: 
buildProblem(
      const Epetra_Map& dofMap,
      const Epetra_Map& overlapped_dofMap,
      FEApp::AbstractPDE_TemplateManager<EvalTypes>& pdeTM,
      std::vector< Teuchos::RCP<FEApp::NodeBC> >& bcs,
      std::vector< Teuchos::RCP<FEApp::AbstractResponseFunction> >& responses,
      const Teuchos::RCP<Epetra_Vector>& u)
{
  // Build PDE equations
  FEApp::HeatNonlinearSourcePDE_TemplateBuilder pdeBuilder(params, paramLib);
  pdeTM.buildObjects(pdeBuilder);

  // Build boundary conditions
  FEApp::ConstantNodeBCStrategy_TemplateBuilder leftBuilder(0, 0, leftBC, 1,
                                                            paramLib);
  FEApp::ConstantNodeBCStrategy_TemplateBuilder rightBuilder(0, 0, rightBC, 2,
                                                             paramLib);
  int left_node = dofMap.MinAllGID();
  int right_node = dofMap.MaxAllGID();
  bcs.resize(2);
  bcs[0] = Teuchos::rcp(new FEApp::NodeBC(dofMap, overlapped_dofMap,
                                          left_node, 1, leftBuilder));
  bcs[1] = Teuchos::rcp(new FEApp::NodeBC(dofMap, overlapped_dofMap,
                                          right_node, 1, rightBuilder));

  // Build response functions
  Teuchos::ParameterList& responseList = params->sublist("Response Functions");
  int num_responses = responseList.get("Number", 0);
  responses.resize(num_responses);
  for (int i=0; i<num_responses; i++) {
     std::ostringstream ss;
     ss << "Response " << i;
     std::string name = responseList.get(ss.str(), "??");

     if (name == "Boundary Flux 1D") {
       double h = 1.0 / (dofMap.NumGlobalElements() - 1);
       responses[i] =
         Teuchos::rcp(new BoundaryFlux1DResponseFunction(left_node,
                                                         right_node,
                                                         0, 1, h,
                                                         dofMap));
     }

     else if (name == "Solution Average")
       responses[i] = Teuchos::rcp(new SolutionAverageResponseFunction());

     else if (name == "Solution Two Norm")
       responses[i] = Teuchos::rcp(new SolutionTwoNormResponseFunction());

     else {
       TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                          std::endl <<
                          "Error!  Unknown response function " << name <<
                          "!" << std::endl << "Supplied parameter list is " <<
                          std::endl << responseList);
     }
  }

  // Build initial solution
  u->PutScalar(1.0);
}
