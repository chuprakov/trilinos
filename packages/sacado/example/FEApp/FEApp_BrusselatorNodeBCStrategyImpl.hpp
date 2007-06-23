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

#include "FEApp_BrusselatorParameters.hpp"

template <typename ScalarT>
FEApp::BrusselatorNodeBCStrategy<ScalarT>::
BrusselatorNodeBCStrategy(
       const ScalarT& alpha_, const ScalarT& beta_,
       const Teuchos::RefCountPtr<Sacado::ScalarParameterLibrary>& paramLib) :
  alpha(alpha_),
  beta(beta_),
  offsets(2),
  pl(paramLib)
{
  offsets[0] = 0;
  offsets[1] = 1;

  // Add "alpha" to parameter library
  std::string name = "Brusselator Alpha";
  if (!pl->isParameter(name))
    pl->addParameterFamily(name, true, false);
  if (!pl->template isParameterForType<ScalarT>(name)) {
    Teuchos::RefCountPtr< BrusselatorAlphaParameter<ScalarT> > tmpa = 
      Teuchos::rcp(new BrusselatorAlphaParameter<ScalarT>(alpha));
    pl->template addEntry<ScalarT>(name, tmpa);
  }

  // Add "beta" to parameter library
  name = "Brusselator Beta";
  if (!pl->isParameter(name))
    pl->addParameterFamily(name, true, false);
  if (!pl->template isParameterForType<ScalarT>(name)) {
    Teuchos::RefCountPtr< BrusselatorBetaParameter<ScalarT> > tmpb = 
      Teuchos::rcp(new BrusselatorBetaParameter<ScalarT>(beta));
    pl->template addEntry<ScalarT>(name, tmpb);
  }
}

template <typename ScalarT>
FEApp::BrusselatorNodeBCStrategy<ScalarT>::
~BrusselatorNodeBCStrategy()
{
}

template <typename ScalarT>
const std::vector<unsigned int>&
FEApp::BrusselatorNodeBCStrategy<ScalarT>::
getOffsets() const
{
  return offsets;
}

template <typename ScalarT>
void
FEApp::BrusselatorNodeBCStrategy<ScalarT>::
evaluateResidual(const std::vector<ScalarT>* dot,
		 const std::vector<ScalarT>& solution,
		 std::vector<ScalarT>& residual) const
{
  alpha = pl->template getValue<ScalarT>("Brusselator Alpha");
  beta = pl->template getValue<ScalarT>("Brusselator Beta");

  residual[0] = solution[0] - alpha;
  residual[1] = solution[1] - beta/alpha;
}
