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

#include "FEApp_GaussianQuadrature2.hpp"
#include <cmath> // for sqrt

FEApp::GaussianQuadrature2::GaussianQuadrature2() :
  qp(2),
  w(2)
{
  qp[0] = -1.0 / std::sqrt(3.0);
  qp[1] = -qp[0];

  w[0] = 1.0;
  w[1] = 1.0;
}

FEApp::GaussianQuadrature2::~GaussianQuadrature2() 
{
}

unsigned int
FEApp::GaussianQuadrature2::numPoints() const
{
  return 2;
}

const std::vector<double>& 
FEApp::GaussianQuadrature2::quadPoints() const
{
  return qp;
}

const std::vector<double>& 
FEApp::GaussianQuadrature2::weights() const
{
  return w;
}
