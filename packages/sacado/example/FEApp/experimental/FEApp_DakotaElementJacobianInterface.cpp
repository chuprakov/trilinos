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

#include "Sacado_ConfigDefs.h"

#ifdef HAVE_DAKOTA

#include "FEApp_DakotaElementJacobianInterface.hpp"
#include "Teuchos_TestForException.hpp"

// Define interface class
FEApp::DakotaElementJacobianInterface::
DakotaElementJacobianInterface(
      const Dakota::ProblemDescDB& problem_db_,
		  const Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> >& pde_,
      const ParamVec& pvec_,
      const Teuchos::RCP<const FEApp::AbstractQuadrature>& quad_,
      unsigned int ndof_,
      double alpha_,
      double beta_) : 
  Dakota::DirectApplicInterface(problem_db_),
  pde(pde_),
  pvec(pvec_),
  quad(quad_),
  numParameters(pvec.size()),
  numEquations(ndof_),
  e(),
  sg_xdot(),
  sg_x(),
  xdot(numEquations),
  x(numEquations),
  f(numEquations),
  alpha(alpha_),
  beta(beta_)
{
  for (unsigned int i=0; i<numEquations; i++) {
    x[i].resize(numEquations);
    xdot[i].resize(numEquations);
    x[i].fastAccessDx(i) = beta;
    xdot[i].fastAccessDx(i) = alpha;
  }
}

void
FEApp::DakotaElementJacobianInterface::
reset(const Teuchos::RCP<const FEApp::AbstractElement>& e_,
      const Teuchos::RCP< std::vector<SGFadType> >& sg_xdot_,
      const Teuchos::RCP< std::vector<SGFadType> >& sg_x_)
{
  e = e_;
  sg_xdot = sg_xdot_;
  sg_x = sg_x_;
}

int 
FEApp::DakotaElementJacobianInterface::
derived_map_ac(const Dakota::String& ac_name)
{
  // test for consistency of problem definition between ModelEval and Dakota
  TEST_FOR_EXCEPTION(numVars != numParameters, logic_error,
                     "FEApp_Dakota Adapter Error: ");
  TEST_FOR_EXCEPTION(numADV != 0, logic_error,
                     "FEApp_Dakota Adapter Error: ");
  TEST_FOR_EXCEPTION(numFns != numEquations*numEquations, logic_error,
                     "FEApp_Dakota Adapter Error: ");
  TEST_FOR_EXCEPTION(hessFlag, logic_error,
                     "FEApp_Dakota Adapter Error: ");
  TEST_FOR_EXCEPTION(gradFlag, logic_error,
                     "FEApp_Dakota Adapter Error: ");

  // Set parameters in FEApp
  for (std::size_t i=0; i<pvec.size(); ++i)
    pvec[i].family->setRealValue<FEApp::JacobianType>(xC[i]);

  std::vector<double> xxC(xC.length());
  for (int i=0; i<xC.length(); i++)
    xxC[i] = xC[i];

  // Evaluate xdot, x on SG basis
  for (unsigned int i=0; i<numEquations; i++) {
    x[i].val() = ((*sg_x)[i]).val().evaluate(xxC);
    if (sg_xdot != Teuchos::null)
      xdot[i].val() = ((*sg_xdot)[i]).val().evaluate(xxC);
  }
  
  // Evaluate element residual
  if (sg_xdot != Teuchos::null)
    pde->evaluateElementResidual(*quad, *e, &xdot, x, f);
  else
    pde->evaluateElementResidual(*quad, *e, NULL, x, f);
  
  // Load function values back into Dakota vector
  for (std::size_t i=0; i<numEquations; i++)
    for (std::size_t j=0; j<numEquations; j++) {
      fnVals[j+i*numEquations] = f[i].fastAccessDx(j);
  }

  return 0;
}

#endif // HAVE_DAKOTA

