// @HEADER
// ***********************************************************************
// 
//              Meros: Segregated Preconditioning Package
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

#include "KayLoghinSchurFactory.h"

#include "SchurFactory.h"

#include "SchurFactoryBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearSolverBase.h"
#include "TSFIdentityOperator.h"
#include "TSFOperatorSourceBase.h"
#include "KayLoghinRightOperatorSource.h"

using namespace SPP;
using namespace TSF;
using std::string;

KayLoghinSchurFactory::KayLoghinSchurFactory(TSFLinearSolver& ApSolver)
  : ApSolver_(ApSolver)
{}



TSFLinearOperator KayLoghinSchurFactory
::getSchurInvApprox(const TSFOperatorSource& OpSrc) const
{
  TSFLinearOperator Xinv; // get the right domain and range space for Xinv?
 
  // Cast OpSrc to KayLoghinRightOperatorSource type.
  const KayLoghinRightOperatorSource* klSrcPtr = 
    dynamic_cast<const KayLoghinRightOperatorSource*>(OpSrc.ptr());

  // Add check that operator source is really the right type.

  // Get saddle operator if we need it:
  //  TSFLinearOperator S = klSrcPtr->getOp();
  //  cout << "Describe S:" << endl;
  //  S.describe();

  // Get Ap operator and set up Apinv using ApSolver
  TSFLinearOperator Ap = klSrcPtr->getAp();
  TSFLinearOperator Apinv = Ap.inverse(ApSolver_);

  // Get Fp operator
  TSFLinearOperator Fp = klSrcPtr->getFp();

  // Using Mp = I for now
  TSFLinearOperator Mpinv = new TSFIdentityOperator(Fp.domain());

  Xinv = -Mpinv * Fp * Apinv;
  return Xinv;
}

string KayLoghinSchurFactory::toString() const
{
	return "KayLoghinSchurFactory toString";
}
