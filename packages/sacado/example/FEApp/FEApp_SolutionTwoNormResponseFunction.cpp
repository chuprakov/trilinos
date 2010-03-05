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
  double nrm;
  x.Norm2(&nrm);

  // Evaluate response g
  if (g != NULL)
    (*g)[0] = nrm;

  // Evaluate tangent of g = dg/dx*dx/dp + dg/dxdot*dxdot/dp + dg/dp
  // dg/dx = 1/||x|| * x^T
  for (int j=0; j<gt.size(); j++)
    if (gt[j] != Teuchos::null) {
      gt[j]->Multiply('T','N',1.0,x,*dx_dp[j],0.0);
      gt[j]->Scale(1.0/nrm);
    }
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
  for (int j=0; j<dg_dp.size(); j++)
    if (dg_dp[j] != Teuchos::null)
      dg_dp[j]->PutScalar(0.0);
}

#if SG_ACTIVE
void 
FEApp::SolutionTwoNormResponseFunction::
init_sg(
  const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& basis,
  const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& quad,
  const Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> >& exp)
{
  sg_basis = basis;
  sg_quad = quad;
}

void
FEApp::SolutionTwoNormResponseFunction::
evaluateSGResponses(const Stokhos::EpetraVectorOrthogPoly* sg_xdot,
		    const Stokhos::EpetraVectorOrthogPoly& sg_x,
		    const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
		    const Teuchos::Array<SGType>* sg_p_vals,
		    Stokhos::EpetraVectorOrthogPoly& sg_g)
{
  // Get basis data
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();

  // Get quadrature data
  const Teuchos::Array< Teuchos::Array<double> >& points = 
    sg_quad->getQuadPoints();
  const Teuchos::Array<double>& weights = sg_quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<double> >& vals = 
    sg_quad->getBasisAtQuadPoints();
  int nqp = points.size();

  // Temporaries for storing inputs, outputs evaluated at quad points
  Epetra_Vector x(sg_x[0]);
  Teuchos::RCP<Epetra_Vector> xdot;
  if (sg_xdot != NULL)
    xdot = Teuchos::rcp(new Epetra_Vector((*sg_xdot)[0]));
  Epetra_Vector g(sg_g[0]);
  Teuchos::Array<double> pvals_orig;
  if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null) {
    pvals_orig.resize(sg_p_vals->size());
    for (int i=0; i<sg_p_vals->size(); i++)
      pvals_orig[i] = (*p[1])[i].baseValue;
  }

  // Compute sg_g via quadrature
  sg_g.init(0.0);
  for (int qp=0; qp<nqp; qp++) {
    sg_x.evaluate(vals[qp], x);
    if (sg_xdot != NULL)
      sg_xdot->evaluate(vals[qp], *xdot);
    if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null) 
      for (int i=0; i<sg_p_vals->size(); i++)
	(*p[1])[i].baseValue = (*sg_p_vals)[i].evaluate(points[qp], vals[qp]);
    evaluateResponses(xdot.get(), x, p, g);
    sg_g.sumIntoAllTerms(weights[qp], vals[qp], norms, g);
  }

  // Restore original parameter values
  if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null)
    for (int i=0; i<sg_p_vals->size(); i++)
      (*p[1])[i].baseValue = pvals_orig[i];
}

void
FEApp::SolutionTwoNormResponseFunction::
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
  // Get basis data
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();

  // Get quadrature data
  const Teuchos::Array< Teuchos::Array<double> >& points = 
    sg_quad->getQuadPoints();
  const Teuchos::Array<double>& weights = sg_quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<double> >& vals = 
    sg_quad->getBasisAtQuadPoints();
  int nqp = points.size();

  // Temporaries for storing inputs, outputs evaluated at quad points
  Epetra_Vector x(sg_x[0]);
  Teuchos::RCP<Epetra_Vector> xdot;
  if (sg_xdot != NULL)
    xdot = Teuchos::rcp(new Epetra_Vector((*sg_xdot)[0]));
  Teuchos::Array<double> pvals_orig;
  if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null) {
    pvals_orig.resize(sg_p_vals->size());
    for (int i=0; i<sg_p_vals->size(); i++)
      pvals_orig[i] = (*p[1])[i].baseValue;
  }
  Teuchos::RCP<Epetra_Vector> g;
  if (sg_g != NULL) {
    g = Teuchos::rcp(new Epetra_Vector((*sg_g)[0]));
    sg_g->init(0.0);
  }
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > gt(sg_gt.size());
  for (int i=0; i<sg_gt.size(); i++)
    if (sg_gt[i] != Teuchos::null) {
      gt[i] = Teuchos::rcp(new Epetra_MultiVector((*sg_gt[i])[0]));
      sg_gt[i]->init(0.0);
    }

  // Compute sg_g, sg_gt via quadrature
  for (int qp=0; qp<nqp; qp++) {
    sg_x.evaluate(vals[qp], x);
    if (sg_xdot != NULL)
      sg_xdot->evaluate(vals[qp], *xdot);
    if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null) 
      for (int i=0; i<sg_p_vals->size(); i++)
	(*p[1])[i].baseValue = (*sg_p_vals)[i].evaluate(points[qp], vals[qp]);
    evaluateTangents(xdot.get(), x, p, deriv_p, dxdot_dp, dx_dp, g.get(), gt);
    if (sg_g != NULL)
      sg_g->sumIntoAllTerms(weights[qp], vals[qp], norms, *g);
    for (int i=0; i<sg_gt.size(); i++)
      if (sg_gt[i] != Teuchos::null)
	sg_gt[i]->sumIntoAllTerms(weights[qp], vals[qp], norms, *gt[i]);
  }

  // Restore original parameter values
  if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null)
    for (int i=0; i<sg_p_vals->size(); i++)
      (*p[1])[i].baseValue = pvals_orig[i];
}

void
FEApp::SolutionTwoNormResponseFunction::
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
  // Get basis data
  const Teuchos::Array<double>& norms = sg_basis->norm_squared();

  // Get quadrature data
  const Teuchos::Array< Teuchos::Array<double> >& points = 
    sg_quad->getQuadPoints();
  const Teuchos::Array<double>& weights = sg_quad->getQuadWeights();
  const Teuchos::Array< Teuchos::Array<double> >& vals = 
    sg_quad->getBasisAtQuadPoints();
  int nqp = points.size();

  // Temporaries for storing inputs, outputs evaluated at quad points
  Epetra_Vector x(sg_x[0]);
  Teuchos::RCP<Epetra_Vector> xdot;
  if (sg_xdot != NULL)
    xdot = Teuchos::rcp(new Epetra_Vector((*sg_xdot)[0]));
  Teuchos::Array<double> pvals_orig;
  if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null) {
    pvals_orig.resize(sg_p_vals->size());
    for (int i=0; i<sg_p_vals->size(); i++)
      pvals_orig[i] = (*p[1])[i].baseValue;
  }
  Teuchos::RCP<Epetra_Vector> g;
  if (sg_g != NULL) {
    g = Teuchos::rcp(new Epetra_Vector((*sg_g)[0]));
    sg_g->init(0.0);
  }
  Teuchos::RCP<Epetra_MultiVector> dg_dx;
  if (sg_dg_dx != NULL) {
    dg_dx = Teuchos::rcp(new Epetra_MultiVector((*sg_dg_dx)[0]));
    sg_dg_dx->init(0.0);
  }
  Teuchos::RCP<Epetra_MultiVector> dg_dxdot;
  if (sg_dg_dxdot != NULL) {
    dg_dxdot = Teuchos::rcp(new Epetra_MultiVector((*sg_dg_dxdot)[0]));
    sg_dg_dxdot->init(0.0);
  }
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > dg_dp(sg_dg_dp.size());
  for (int i=0; i<sg_dg_dp.size(); i++)
    if (sg_dg_dp[i] != Teuchos::null) {
      dg_dp[i] = Teuchos::rcp(new Epetra_MultiVector((*sg_dg_dp[i])[0]));
      sg_dg_dp[i]->init(0.0);
    }

  // Compute sg_g, sg_dg_dx, sg_dg_dxdot, sg_dg_dp via quadrature
  for (int qp=0; qp<nqp; qp++) {
    sg_x.evaluate(vals[qp], x);
    if (sg_xdot != NULL)
      sg_xdot->evaluate(vals[qp], *xdot);
    if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null) 
      for (int i=0; i<sg_p_vals->size(); i++)
	(*p[1])[i].baseValue = (*sg_p_vals)[i].evaluate(points[qp], vals[qp]);
    evaluateGradients(xdot.get(), x, p, deriv_p, g.get(), dg_dx.get(), 
		      dg_dxdot.get(), dg_dp);
    if (sg_g != NULL)
      sg_g->sumIntoAllTerms(weights[qp], vals[qp], norms, *g);
    if (sg_dg_dx != NULL)
      sg_dg_dx->sumIntoAllTerms(weights[qp], vals[qp], norms, *dg_dx);
    if (sg_dg_dxdot != NULL)
      sg_dg_dxdot->sumIntoAllTerms(weights[qp], vals[qp], norms, *dg_dxdot);
    for (int i=0; i<sg_dg_dp.size(); i++)
      if (sg_dg_dp[i] != Teuchos::null)
	sg_dg_dp[i]->sumIntoAllTerms(weights[qp], vals[qp], norms, *dg_dp[i]);
  }

  // Restore original parameter values
  if (sg_p_vals != NULL && p.size() > 1 && p[1] != Teuchos::null)
    for (int i=0; i<sg_p_vals->size(); i++)
      (*p[1])[i].baseValue = pvals_orig[i];
}
#endif
