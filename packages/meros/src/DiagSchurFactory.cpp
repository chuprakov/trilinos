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

#include "DiagSchurFactory.h"

#include "SchurFactory.h"

#include "SchurFactoryBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFLinearOperatorBase.h"
#include "TSFLinearSolverBase.h"
#include "TSFIdentityOperator.h"
#include "PetraMatrix.h"
#include "TSFOperatorSourceBase.h"
#include "DiagRightOperatorSource.h"
#include "Aztec2TSF.h"

using namespace Meros;
using namespace TSF;
using std::string;

DiagSchurFactory::DiagSchurFactory(TSFLinearSolver schurSolver)
  : schurSolver_(schurSolver)
{}



TSFLinearOperator DiagSchurFactory
::getSchurInvApprox(const TSFOperatorSource& OpSrc) const
{
  cerr << "In DiagSchurFactor: getSchurInvApprox" << endl;
  TSFLinearOperator Xinv; // get the right domain and range space for Xinv?
 
  // Cast OpSrc to DiagRightOperatorSource type.
  const DiagRightOperatorSource* diagSrcPtr = 
    dynamic_cast<const DiagRightOperatorSource*>(OpSrc.ptr());

  // Add check that operator source is really the right type.

  // Get saddle operator 
  TSFLinearOperator S = diagSrcPtr->getOp();

  // Get needed blocks from saddle operator
  TSFLinearOperator Bt = S.getBlock(0,1);
  TSFLinearOperator B = S.getBlock(1,0);
  // TSFLinearOperator C = S.getBlock(1,1);

  TSFLinearOperator F = S.getBlock(0,0);

  // Check for C

  // Get negative Dinv operator
  TSFLinearOperator Dinv = diagSrcPtr->getDinv();

  // If using Multigrid,
  // do explicit matrix multiplies and adds to get Schur approx
  // Need explicit operations so we have a matrix for multigrid in the end
  // Could put a check here: if schurSolver is multigrid, use explicit mults
  // otherwise, the regular TSF implicit mult is OK and cheaper

  TSFLinearOperator schurOp;
  TSFLinearOperator tempOp1;
  TSF_MatrixMult(Dinv, Bt, tempOp1); // tempOp1 = Dinv*Bt
  TSF_MatrixMult(B, tempOp1, schurOp); // tempOp2 = B*tempOp1
  cerr << "got here after matmult" << endl;

  S.describe();
  
  Dinv.describe();
  Bt.describe();
  tempOp1.describe();

  B.describe();
  schurOp.describe();


  TSFLinearOperator tempOp2;
  TSF_MatrixMult(B, Bt, tempOp2); // tempOp2 = B*tempOp1

  tempOp2.describe();



  // If not ML, can use TSF's regular deferred multiplications

  //  schurOp = B*Dinv*Bt;

  Xinv = schurOp.inverse(schurSolver_);

  return Xinv;

}

string DiagSchurFactory::toString() const
{
	return "DiagSchurFactory toString";
}
