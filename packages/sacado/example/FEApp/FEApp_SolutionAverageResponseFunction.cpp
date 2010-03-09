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

#include "FEApp_SolutionAverageResponseFunction.hpp"

FEApp::SolutionAverageResponseFunction::
SolutionAverageResponseFunction()
{
}

FEApp::SolutionAverageResponseFunction::
~SolutionAverageResponseFunction()
{
}

unsigned int
FEApp::SolutionAverageResponseFunction::
numResponses() const 
{
  return 1;
}

void
FEApp::SolutionAverageResponseFunction::
evaluateResponses(const Epetra_Vector* xdot,
		  const Epetra_Vector& x,
		  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
		  Epetra_Vector& g)
{
  x.MeanValue(&g[0]);
}

void
FEApp::SolutionAverageResponseFunction::
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
    x.MeanValue(&(*g)[0]);

  // Evaluate tangent of g = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
  for (int j=0; j<gt.size(); j++)
    if (gt[j] != Teuchos::null)
      for (int i=0; i<dx_dp[i]->NumVectors(); i++)
	(*dx_dp[j])(i)->MeanValue(&(*gt[j])[i][0]);
}

void
FEApp::SolutionAverageResponseFunction::
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
    x.MeanValue(&(*g)[0]);

  // Evaluate dg/dx
  if (dg_dx != NULL)
    dg_dx->PutScalar(1.0 / x.GlobalLength());

  // Evaluate dg/dxdot
  if (dg_dxdot != NULL)
    dg_dxdot->PutScalar(0.0);

  // Evaluate dg/dp
  for (int j=0; j<dg_dp.size(); j++)
    if (dg_dp[j] != Teuchos::null)
      dg_dp[j]->PutScalar(0.0);
}

#if SG_ACTIVE
void
FEApp::SolutionAverageResponseFunction::
evaluateSGResponses(const Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_xdot,
		    const Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_x,
		    const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
		    const Teuchos::Array<SGType>* sg_p_vals,
		    Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_g)
{
  int sz = sg_x.size();
  for (int i=0; i<sz; i++)
    sg_x[i].MeanValue(&sg_g[i][0]);
}

void
FEApp::SolutionAverageResponseFunction::
evaluateSGTangents(
  const Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_xdot,
  const Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_x,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  const Teuchos::Array<SGType>* sg_p_vals,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
  Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_g,
  const Teuchos::Array< Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_MultiVector> > >& sg_gt)
{
  int sz = sg_x.size();
    
  // Evaluate response g
  if (sg_g != NULL)
    for (int i=0; i<sz; i++)
      sg_x[i].MeanValue(&((*sg_g)[i][0]));

  // Evaluate tangent of g = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
  for (int j=0; j<sg_gt.size(); j++) {
    if (sg_gt[j] != Teuchos::null)
      sg_gt[j]->init(0.0);
      for (int i=0; i<dx_dp[i]->NumVectors(); i++)
	(*dx_dp[j])(i)->MeanValue(&(*sg_gt[j])[0][i][0]);
  }
}

void
FEApp::SolutionAverageResponseFunction::
evaluateSGGradients(
  const Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_xdot,
  const Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_x,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  const Teuchos::Array<SGType>* sg_p_vals,
  Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_g,
  Stokhos::VectorOrthogPoly<Epetra_MultiVector>* sg_dg_dx,
  Stokhos::VectorOrthogPoly<Epetra_MultiVector>* sg_dg_dxdot,
  const Teuchos::Array< Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_MultiVector> > >& sg_dg_dp)
{
  int sz = sg_x.size();
    
  // Evaluate response g
  if (sg_g != NULL) {
    for (int i=0; i<sz; i++)
      sg_x[i].MeanValue(&((*sg_g)[i][0]));
  }

  // Evaluate dg/dx
  if (sg_dg_dx != NULL) {
    sg_dg_dx->init(0.0);
    (*sg_dg_dx)[0].PutScalar(1.0 / sg_x[0].GlobalLength());
  }

  // Evaluate dg/dxdot
  if (sg_dg_dxdot != NULL)
    sg_dg_dxdot->init(0.0);

  // Evaluate dg/dp
  for (int j=0; j<sg_dg_dp.size(); j++)
    if (sg_dg_dp[j] != Teuchos::null)
      sg_dg_dp[j]->init(0.0);
}
#endif
