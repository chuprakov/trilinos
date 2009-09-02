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
#include "FEApp_QuadratureFactory.hpp"
#include "FEApp_GaussianQuadrature2.hpp"

FEApp::QuadratureFactory::QuadratureFactory(
	    const Teuchos::RCP<Teuchos::ParameterList>& quadParams_) :
  quadParams(quadParams_)
{
}

Teuchos::RCP<FEApp::AbstractQuadrature>
FEApp::QuadratureFactory::create()
{
  Teuchos::RCP<FEApp::AbstractQuadrature> strategy;

  std::string& method = quadParams->get("Method", "Gaussian");
  if (method == "Gaussian") {
    int num_points = quadParams->get("Num Points", 2);

    if (num_points == 2) {
      strategy = Teuchos::rcp(new FEApp::GaussianQuadrature2);
    }
    else {
      TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                         std::endl << 
                         "Error!  Number of quadrature points = " << 
                         num_points << 
                         " is not supported for Gaussian quadrature!" << 
                         std::endl << "Supplied parameter list is " << 
                         std::endl << *quadParams);
    }
  }
  else {
    TEST_FOR_EXCEPTION(true, Teuchos::Exceptions::InvalidParameter,
                       std::endl << 
                       "Error!  Unknown quadrature method " << method << 
                       "!" << std::endl << "Supplied parameter list is " << 
                       std::endl << *quadParams);
  }
  
  return strategy;
}
