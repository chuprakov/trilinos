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

#include "FEApp_SolutionTwoNormResponseFunction.hpp"
#include "Epetra_Comm.h"

FEApp::SolutionTwoNormResponseFunction::
SolutionTwoNormResponseFunction()
{
}

FEApp::SolutionTwoNormResponseFunction::
~SolutionTwoNormResponseFunction()
{
}

unsigned int
FEApp::SolutionTwoNormResponseFunction::
numResponses() const 
{
  return 1;
}

void
FEApp::SolutionTwoNormResponseFunction::
evaluateResponses(const Epetra_Vector* xdot,
		  const Epetra_Vector& x,
		  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
		  Epetra_Vector& g)
{
  x.Norm2(&g[0]);
}

void
FEApp::SolutionTwoNormResponseFunction::
evaluateTangents(
	   const Epetra_Vector* xdot,
	   const Epetra_Vector& x,
	   const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
	   const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
	   const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
	   const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
	   Epetra_Vector* g,
	   const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& gt)
{
  // Evaluate response g
  if (g != NULL)
    x.Norm2(&(*g)[0]);

  // Evaluate tangent of g = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
  // dg/dx = 1/||x|| * x^T
  for (unsigned int j=0; j<gt.size(); j++)
    if (gt[j] != Teuchos::null)
      gt[j]->Multiply('T','N',1.0,x,*dx_dp[j],0.0);
}

void
FEApp::SolutionTwoNormResponseFunction::
evaluateGradients(
	  const Epetra_Vector* xdot,
	  const Epetra_Vector& x,
	  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
	  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
	  Epetra_Vector* g,
	  Epetra_MultiVector* dg_dx,
	  Epetra_MultiVector* dg_dxdot,
	  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dg_dp)
{

  // Evaluate response g
  if (g != NULL)
    x.Norm2(&(*g)[0]);

  // Evaluate dg/dx
  if (dg_dx != NULL) {
    double nrm;
    if (g != NULL)
      nrm = (*g)[0];
    else
      x.Norm2(&nrm);
    dg_dx->Scale(1.0/nrm,x);
  }

  // Evaluate dg/dxdot
  if (dg_dxdot != NULL)
    dg_dxdot->PutScalar(0.0);

  // Evaluate dg/dp
  for (unsigned int j=0; j<dg_dp.size(); j++)
    if (dg_dp[j] != Teuchos::null)
      dg_dp[j]->PutScalar(0.0);
}

#ifdef SG_ACTIVE
void
FEApp::SolutionTwoNormResponseFunction::
evaluateSGResponses(const Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_xdot,
		    const Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_x,
		    const ParamVec* p,
		    const ParamVec* sg_p,
		    const Teuchos::Array<SGType>* sg_p_vals,
		    Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_g)
{
  int sz = sg_x.size();
  int N = sg_x[0].MyLength();
  SGType nrm_local = 0.0;
  SGType x(sz);
  for (int i=0; i<N; i++) {
    for (int k=0; k<sz; k++)
      x.fastAccessCoeff(k) = sg_x[k][i];
    nrm_local += x*x;
  }
#ifdef HAVE_MPI
  SGType nrm(sz);
  sg_x[0].Map().Comm().SumAll(nrm_local.coeff(), nrm.coeff(), sz);
#else
  SGType& nrm = nrm_local;
#endif
  nrm = std::sqrt(nrm);
  for (int k=0; k<sz; k++)
    sg_g[k][0] = nrm.fastAccessCoeff(k);
}
#endif
