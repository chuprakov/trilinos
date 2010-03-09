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

#include "FEApp_BoundaryFlux1DResponseFunction.hpp"

FEApp::BoundaryFlux1DResponseFunction::
BoundaryFlux1DResponseFunction(unsigned int left_gid,
                               unsigned int right_gid,
                               unsigned int eqn,
                               unsigned int num_eqns,
                               double grid_spacing_,
                               const Epetra_Map& dofMap) :
  grid_spacing(grid_spacing_),
  importer(NULL),
  bv(NULL)
{
  // Compute GID's of left/right DOF's
  unsigned int left_dof = num_eqns*left_gid+eqn;
  unsigned int right_dof = num_eqns*right_gid+eqn;

  // Build importer to bring in left/right DOF's to all proc's
  int gids[4];
  gids[0] = left_dof;
  gids[1] = left_dof+num_eqns;
  gids[2] = right_dof-num_eqns;
  gids[3] = right_dof;
  boundaryMap = new Epetra_Map(4, 4, gids, 0, dofMap.Comm());
  importer = new Epetra_Import(*boundaryMap, dofMap);
  bv = new Epetra_Vector(*boundaryMap);
}

FEApp::BoundaryFlux1DResponseFunction::
~BoundaryFlux1DResponseFunction()
{
  delete boundaryMap;
  delete importer;
  delete bv;
}

unsigned int
FEApp::BoundaryFlux1DResponseFunction::
numResponses() const 
{
  return 2;
}

void
FEApp::BoundaryFlux1DResponseFunction::
evaluateResponses(const Epetra_Vector* xdot,
		  const Epetra_Vector& x,
		  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
		  Epetra_Vector& g)
{
  // Import boundary values
  bv->Import(x, *importer, Insert);

  // Compute fluxes
  g[0] = ((*bv)[1] - (*bv)[0]) / grid_spacing;
  g[1] = ((*bv)[3] - (*bv)[2]) / grid_spacing;
}

void
FEApp::BoundaryFlux1DResponseFunction::
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
  if (g != NULL) {
    bv->Import(x, *importer, Insert);
    (*g)[0] = ((*bv)[1] - (*bv)[0]) / grid_spacing;
    (*g)[1] = ((*bv)[3] - (*bv)[2]) / grid_spacing;
  }

  // Evaluate tangent of g = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
  for (int j=0; j<gt.size(); j++)
    if (gt[j] != Teuchos::null) {
      Epetra_MultiVector bvt(*boundaryMap, dx_dp[j]->NumVectors());
      bvt.Import(*dx_dp[j], *importer, Insert);
      for (int i=0; i<dx_dp[j]->NumVectors(); i++) {
	(*gt[j])[i][0] = (bvt[i][1] - bvt[i][0]) / grid_spacing;
	(*gt[j])[i][1] = (bvt[i][3] - bvt[i][2]) / grid_spacing;
      }
    }
}

void
FEApp::BoundaryFlux1DResponseFunction::
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
  if (g != NULL) {
    bv->Import(x, *importer, Insert);
    (*g)[0] = ((*bv)[1] - (*bv)[0]) / grid_spacing;
    (*g)[1] = ((*bv)[3] - (*bv)[2]) / grid_spacing;
  }

  // Evaluate dg/dx
  if (dg_dx != NULL) {
    Epetra_MultiVector bv_dx(*boundaryMap, 2);
    bv_dx.PutScalar(0.0);
    bv_dx[0][0] = -1.0 / grid_spacing;
    bv_dx[0][1] =  1.0 / grid_spacing;
    bv_dx[1][2] = -1.0 / grid_spacing;
    bv_dx[1][3] =  1.0 / grid_spacing;
    dg_dx->PutScalar(0.0);
    dg_dx->Export(bv_dx, *importer, Insert);
  }

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
FEApp::BoundaryFlux1DResponseFunction::
evaluateSGResponses(const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
		    const Stokhos::EpetraVectorOrthogPoly& sg_x,
		    const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
		    const Teuchos::Array<SGType>* sg_p_vals,
		    Stokhos::EpetraVectorOrthogPoly& sg_g)
{
  int sz = sg_x.size();
  for (int i=0; i<sz; i++) {

    // Import boundary values
    bv->Import(sg_x[i], *importer, Insert);
    
    // Compute fluxes
    sg_g[i][0] = ((*bv)[1] - (*bv)[0]) / grid_spacing;
    sg_g[i][1] = ((*bv)[3] - (*bv)[2]) / grid_spacing;

  }
}

void
FEApp::BoundaryFlux1DResponseFunction::
evaluateSGTangents(
  const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
  const Stokhos::EpetraVectorOrthogPoly& sg_x,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  const Teuchos::Array<SGType>* sg_p_vals,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
  const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
  Stokhos::EpetraVectorOrthogPoly* sg_g,
  const Teuchos::Array< Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly > >& sg_gt)
{
  if (sg_g != NULL)
    evaluateSGResponses(sg_xdot, sg_x, p, sg_p_vals, *sg_g);

  // For tangent, only mean is non-zero, which is the tangent of the mean
  const Epetra_Vector *xdot_mean = NULL;
  if (sg_xdot != NULL)
    xdot_mean = sg_xdot->getCoeffPtr(0).get();
  const Epetra_Vector& x_mean = sg_x[0];
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > gt_mean(sg_gt.size());
  for (int i=0; i<sg_gt.size(); i++)
    if (sg_gt[i] != Teuchos::null) {
      sg_gt[i]->init(0.0);
      gt_mean[i] = sg_gt[i]->getCoeffPtr(0);
    }
  evaluateTangents(xdot_mean, x_mean, p, deriv_p, dxdot_dp, dx_dp, NULL, 
		   gt_mean);
}

void
FEApp::BoundaryFlux1DResponseFunction::
evaluateSGGradients(
  const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
  const Stokhos::EpetraVectorOrthogPoly& sg_x,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
  const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
  const Teuchos::Array<SGType>* sg_p_vals,
  Stokhos::EpetraVectorOrthogPoly* sg_g,
  Stokhos::EpetraMultiVectorOrthogPoly* sg_dg_dx,
  Stokhos::EpetraMultiVectorOrthogPoly* sg_dg_dxdot,
  const Teuchos::Array< Teuchos::RCP<Stokhos::EpetraMultiVectorOrthogPoly > >& sg_dg_dp)
{
  if (sg_g != NULL)
    evaluateSGResponses(sg_xdot, sg_x, p, sg_p_vals, *sg_g);

  // For gradients, only mean is non-zero, which is the gradient of the mean
  const Epetra_Vector *xdot_mean = NULL;
  if (sg_xdot != NULL)
    xdot_mean = sg_xdot->getCoeffPtr(0).get();
  const Epetra_Vector& x_mean = sg_x[0];
  Epetra_MultiVector* dg_dx_mean = NULL;
  if (sg_dg_dx != NULL) {
    sg_dg_dx->init(0.0);
    dg_dx_mean = sg_dg_dx->getCoeffPtr(0).get();
  }
  Epetra_MultiVector* dg_dxdot_mean = NULL;
  if (sg_dg_dxdot != NULL) {
    sg_dg_dxdot->init(0.0);
    dg_dxdot_mean = sg_dg_dxdot->getCoeffPtr(0).get();
  }
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > dg_dp_mean(sg_dg_dp.size());
  for (int i=0; i<sg_dg_dp.size(); i++)
    if (sg_dg_dp[i] != Teuchos::null) {
      sg_dg_dp[i]->init(0.0);
      dg_dp_mean[i] = sg_dg_dp[i]->getCoeffPtr(0);
    }
  evaluateGradients(xdot_mean, x_mean, p, deriv_p, NULL, 
		    dg_dx_mean, dg_dxdot_mean, dg_dp_mean);
}
#endif
