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

#include "BJBSchurFactory.h"

#include "SchurFactory.h"

#include "SchurFactoryBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearSolverBase.h"
#include "TSFIdentityOperator.h"
#include "TSFOperatorSourceBase.h"
#include "BJBRightOperatorSource.h"

using namespace Meros;
using namespace TSF;
using std::string;

BJBSchurFactory::BJBSchurFactory(TSFLinearSolver& ApSolver)
  : ApSolver_(ApSolver)
{}

TSFLinearOperator BJBSchurFactory
::getSchurInvApprox(const TSFOperatorSource& OpSrc) const
{
  TSFLinearOperator Xinv; // get the right domain and range space for Xinv?
  
  // Cast OpSrc to BJBRightOperatorSource type.
  const BJBRightOperatorSource* bjbSrcPtr = 
    dynamic_cast<const BJBRightOperatorSource*>(OpSrc.ptr());

  // Add check that operator source is really the right type.

  // Get saddle operator if we need it:
  TSFLinearOperator S = bjbSrcPtr->getOp();

  TSFLinearOperator F = S.getBlock(0,0);
  TSFLinearOperator Bt = S.getBlock(0,1);
  TSFLinearOperator B = S.getBlock(1,0);
  TSFLinearOperator C = S.getBlock(1,1);

  //  cout << "Describe S:" << endl;
  //  S.describe();

 if (C.isZeroOperator())
   {
  // Get Ap operator and set up Apinv using ApSolver
     cerr << "\n C is zero!";
  TSFLinearOperator Ap = bjbSrcPtr->getAp();
  TSFLinearOperator Apinv = Ap.inverse(ApSolver_);
  Xinv = -Apinv * B * F * Bt * Apinv;
   }
 else
   {
     TSFReal alpha, beta;
     bool noscal = bjbSrcPtr->getScal();
     if(noscal)
       {
   alpha = bjbSrcPtr->getAlpha();
   beta = bjbSrcPtr->getBeta();
   cerr << "\nC is not zero! creating preconditioner with beta= " << beta << " and alpha= " << alpha << endl;
       }
     else
       {
       alpha = 1.0;
       beta = 1.0;  
       cerr << "\n Alpha and beta not given, calculating ... alpha is " << alpha << " and beta is " << beta << endl;
       }
  TSFLinearOperator Dinv = bjbSrcPtr->getDinv();
  TSFLinearOperator Cinv = C.inverse(ApSolver_);
  Xinv = -beta*Dinv + alpha*alpha*Cinv*(B*F*Bt)*Cinv;
   }
 Xinv.describe();
  return Xinv;
}

string BJBSchurFactory::toString() const
{
	return "BJBSchurFactory toString";
}
