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

#include "FEApp_BrusselatorProblem.hpp"
#include "FEApp_BrusselatorNodeBCStrategy.hpp"

FEApp::BrusselatorProblem::
BrusselatorProblem(
      const Teuchos::RCP<Teuchos::ParameterList>& params,
      const Teuchos::RCP<ParamLib>& paramLib_) :
  paramLib(paramLib_)
{
  alpha = params->get("alpha", 1.0);
  beta = params->get("beta", 1.0);
  D1 = params->get("D1", 1.0);
  D2 = params->get("D2", 1.0);
}

FEApp::BrusselatorProblem::
~BrusselatorProblem()
{
}

unsigned int
FEApp::BrusselatorProblem::
numEquations() const
{
  return 2;
}

void
FEApp::BrusselatorProblem::
buildProblem(const Epetra_Map& dofMap,
             const Epetra_Map& overlapped_dofMap,
             FEApp::AbstractPDE_TemplateManager<EvalTypes>& pdeTM,
             std::vector< Teuchos::RCP<FEApp::NodeBC> >& bcs,
             const Teuchos::RCP<Epetra_Vector>& u)
{
  // Build PDE equations
  FEApp::BrusselatorPDE_TemplateBuilder pdeBuilder(alpha, beta, D1, D2, 
                                                   paramLib);
  pdeTM.buildObjects(pdeBuilder);

  // Build boundary conditions
  FEApp::BrusselatorNodeBCStrategy_TemplateBuilder bcBuilder(alpha, beta,
                                                             paramLib);
  int left_node = dofMap.MinAllGID();
  int right_node = 
    (dofMap.MaxAllGID() - dofMap.MinAllGID())/2 + dofMap.MinAllGID();
  bcs.resize(2);
  bcs[0] = Teuchos::rcp(new FEApp::NodeBC(dofMap, overlapped_dofMap,
					  left_node, 2, bcBuilder));
  bcs[1] = Teuchos::rcp(new FEApp::NodeBC(dofMap, overlapped_dofMap,
					  right_node, 2, bcBuilder));

  // Build initial solution
//   for (int i=0; i<u->MyLength()/2; i++) {
//     (*u)[2*i]   = alpha;
//     (*u)[2*i+1] = beta/alpha;
//   }
  u->PutScalar(0.0);
}

