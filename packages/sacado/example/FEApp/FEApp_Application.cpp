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

#include "FEApp_Application.hpp"
#include "FEApp_ProblemFactory.hpp"
#include "FEApp_QuadratureFactory.hpp"
#include "FEApp_DiscretizationFactory.hpp"
#include "Epetra_LocalMap.h"
#if SG_ACTIVE
#include "FEApp_SGGaussQuadResidualGlobalFill.hpp"
#include "FEApp_SGGaussQuadJacobianGlobalFill.hpp"
#include "Stokhos_MatrixFreeEpetraOp.hpp"
#include "Stokhos_MeanEpetraOp.hpp"
#endif
#include "Teuchos_TimeMonitor.hpp"

FEApp::Application::Application(
		   const std::vector<double>& coords,
		   const Teuchos::RCP<const Epetra_Comm>& comm,
		   const Teuchos::RCP<Teuchos::ParameterList>& params_,
		   bool is_transient,
		   const Epetra_Vector* initial_soln) :
  params(params_),
  transient(is_transient)
{
  // Create parameter library
  paramLib = Teuchos::rcp(new ParamLib);

  // Create problem object
  Teuchos::RCP<Teuchos::ParameterList> problemParams = 
    Teuchos::rcp(&(params->sublist("Problem")),false);
  FEApp::ProblemFactory problemFactory(problemParams, paramLib);
  Teuchos::RCP<FEApp::AbstractProblem> problem = 
    problemFactory.create();

  // Get number of equations
  unsigned int num_equations = problem->numEquations();

  // Create quadrature object
  Teuchos::RCP<Teuchos::ParameterList> quadParams = 
    Teuchos::rcp(&(params->sublist("Quadrature")),false);
  FEApp::QuadratureFactory quadFactory(quadParams);
  quad = quadFactory.create();

  // Create discretization object
  Teuchos::RCP<Teuchos::ParameterList> discParams = 
    Teuchos::rcp(&(params->sublist("Discretization")),false);
  FEApp::DiscretizationFactory discFactory(discParams);
  disc = discFactory.create(coords, num_equations, comm);
  disc->createMesh();
  disc->createMaps();
  disc->createJacobianGraphs();

  // Create Epetra objects
  importer = Teuchos::rcp(new Epetra_Import(*(disc->getOverlapMap()), 
                                            *(disc->getMap())));
  exporter = Teuchos::rcp(new Epetra_Export(*(disc->getOverlapMap()), 
                                            *(disc->getMap())));
  overlapped_x = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  if (transient)
    overlapped_xdot = 
      Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_f = Teuchos::rcp(new Epetra_Vector(*(disc->getOverlapMap())));
  overlapped_jac = 
    Teuchos::rcp(new Epetra_CrsMatrix(Copy, 
                                      *(disc->getOverlapJacobianGraph())));

  // Initialize problem
  initial_x = Teuchos::rcp(new Epetra_Vector(*(disc->getMap())));
  problem->buildProblem(*(disc->getMap()), *(disc->getOverlapMap()), 
                        pdeTM, bc, responses, initial_x);
  if (initial_soln != NULL)
    initial_x->Scale(1.0, *initial_soln);
  typedef FEApp::AbstractPDE_TemplateManager<EvalTypes>::iterator iterator;
  int nqp = quad->numPoints();
  int nn = disc->getNumNodesPerElement();
  for (iterator it = pdeTM.begin(); it != pdeTM.end(); ++it)
    it->init(nqp, nn);

  // Create response map
  unsigned int total_num_responses = 0;
  for (unsigned int i=0; i<responses.size(); i++)
    total_num_responses += responses[i]->numResponses();
  if (total_num_responses > 0)
    response_map = Teuchos::rcp(new Epetra_LocalMap(total_num_responses, 0,
                                                    *comm));

#if SG_ACTIVE
  bool enable_sg = params->get("Enable Stochastic Galerkin",false);
  if (enable_sg) {
    sg_expansion = params->get< Teuchos::RCP<Stokhos::OrthogPolyExpansion<int,double> > >("Stochastic Galerkin expansion");
    if (params->isParameter("Stochastic Galerkin quadrature"))
      sg_quad = params->get< Teuchos::RCP<const Stokhos::Quadrature<int,double> > >("Stochastic Galerkin quadrature");
    sg_basis = sg_expansion->getBasis();

    // Create Epetra orthogonal polynomial objects
    sg_overlapped_x = 
      Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Vector>(sg_basis,
							  *overlapped_x));
    if (transient)
      sg_overlapped_xdot = 
	Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Vector>(sg_basis,
							    *overlapped_xdot));
    sg_overlapped_f = 
      Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_Vector>(sg_basis,
							  *overlapped_f));
    // Delay creation of sg_overlapped_jac until needed
  }
#endif
}

FEApp::Application::~Application()
{
}

Teuchos::RCP<const Epetra_Map>
FEApp::Application::getMap() const
{
  return disc->getMap();
}

Teuchos::RCP<const Epetra_Map>
FEApp::Application::getResponseMap() const
{
  return response_map;
}

Teuchos::RCP<const Epetra_CrsGraph>
FEApp::Application::getJacobianGraph() const
{
  return disc->getJacobianGraph();
}

Teuchos::RCP<const Epetra_Vector>
FEApp::Application::getInitialSolution() const
{
  return initial_x;
}

Teuchos::RCP<ParamLib> 
FEApp::Application::getParamLib()
{
  return paramLib;
}

bool
FEApp::Application::isTransient() const
{
  return transient;
}

Teuchos::RCP<Epetra_Operator>
FEApp::Application::createW() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
					   *(disc->getJacobianGraph())));
}

Teuchos::RCP<Epetra_Operator>
FEApp::Application::createPrec() const
{
  return Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
					   *(disc->getJacobianGraph())));
}

void
FEApp::Application::computeGlobalResidual(
                          const Epetra_Vector* xdot,
			  const Epetra_Vector& x,
			  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
			  Epetra_Vector& f)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalResidual");

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (unsigned int i=0; i<p.size(); i++) {
    if (p[i] != Teuchos::null)
      for (unsigned int j=0; j<p[i]->size(); j++)
	(*(p[i]))[j].family->setRealValueForAllTypes((*(p[i]))[j].baseValue);
  }

  // Zero out overlapped residual
  overlapped_f->PutScalar(0.0);
  f.PutScalar(0.0);

  // Create residual init/post op
  Teuchos::RCP<FEApp::ResidualOp> op = 
    Teuchos::rcp(new FEApp::ResidualOp(overlapped_xdot, overlapped_x, 
                                       overlapped_f));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::ResidualType> > pde = 
    pdeTM.getAsObject<FEApp::ResidualType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::ResidualType> globalFill(disc->getMesh(), quad, 
                                                    pde, bc, transient);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  f.Export(*overlapped_f, *exporter, Add);
}

void
FEApp::Application::computeGlobalJacobian(
			    double alpha, double beta,
			    const Epetra_Vector* xdot,
			    const Epetra_Vector& x,
			    const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
			    Epetra_Vector* f,
			    Epetra_Operator& jacOp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalJacobian");

  // Cast jacOp to an Epetra_CrsMatrix
  Epetra_CrsMatrix& jac = dynamic_cast<Epetra_CrsMatrix&>(jacOp);

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (unsigned int i=0; i<p.size(); i++) {
    if (p[i] != Teuchos::null)
      for (unsigned int j=0; j<p[i]->size(); j++)
	(*(p[i]))[j].family->setRealValueForAllTypes((*(p[i]))[j].baseValue);
  }

  // Zero out overlapped residual
  Teuchos::RCP<Epetra_Vector> overlapped_ff;
  if (f != NULL) {
    overlapped_ff = overlapped_f;
    overlapped_ff->PutScalar(0.0);
    f->PutScalar(0.0);
  }

  // Zero out Jacobian
  overlapped_jac->PutScalar(0.0);
  jac.PutScalar(0.0);

  // Create Jacobian init/post op
  Teuchos::RCP<FEApp::JacobianOp> op
    = Teuchos::rcp(new FEApp::JacobianOp(alpha, beta, overlapped_xdot, 
                                         overlapped_x, overlapped_ff, 
                                         overlapped_jac));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > pde = 
    pdeTM.getAsObject<FEApp::JacobianType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::JacobianType> globalFill(disc->getMesh(), 
                                                    quad, pde, bc, 
                                                    transient);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (f != NULL)
    f->Export(*overlapped_f, *exporter, Add);

  // Assemble global Jacobian
  jac.Export(*overlapped_jac, *exporter, Add);

  jac.FillComplete(true);
}

void
FEApp::Application::computeGlobalPreconditioner(
			    double alpha, double beta,
			    const Epetra_Vector* xdot,
			    const Epetra_Vector& x,
			    const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
			    Epetra_Vector* f,
			    Epetra_Operator& jacOp)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalPreconditioner");

  // Cast jacOp to an Epetra_CrsMatrix
  Epetra_CrsMatrix& jac = dynamic_cast<Epetra_CrsMatrix&>(jacOp);

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (unsigned int i=0; i<p.size(); i++) {
    if (p[i] != Teuchos::null)
      for (unsigned int j=0; j<p[i]->size(); j++)
	(*(p[i]))[j].family->setRealValueForAllTypes((*(p[i]))[j].baseValue);
  }

  // Zero out overlapped residual
  Teuchos::RCP<Epetra_Vector> overlapped_ff;
  if (f != NULL) {
    overlapped_ff = overlapped_f;
    overlapped_ff->PutScalar(0.0);
    f->PutScalar(0.0);
  }

  // Zero out Jacobian
  overlapped_jac->PutScalar(0.0);
  jac.PutScalar(0.0);

  // Create Jacobian init/post op
  Teuchos::RCP<FEApp::JacobianOp> op = 
    Teuchos::rcp(new FEApp::JacobianOp(alpha, beta, overlapped_xdot, 
                                       overlapped_x, overlapped_ff, 
                                       overlapped_jac));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > pde = 
    pdeTM.getAsObject<FEApp::JacobianType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::JacobianType> globalFill(disc->getMesh(), 
                                                    quad, pde, bc, 
                                                    transient);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (f != NULL)
    f->Export(*overlapped_f, *exporter, Add);

  // Assemble global Jacobian
  jac.Export(*overlapped_jac, *exporter, Add);

  jac.FillComplete(true);
}

void
FEApp::Application::computeGlobalTangent(
			      double alpha, double beta,
			      bool sum_derivs,
			      const Epetra_Vector* xdot,
			      const Epetra_Vector& x,
			      const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
			      ParamVec* deriv_p,
			      const Epetra_MultiVector* Vx,
			      const Teuchos::SerialDenseMatrix<int,double>* Vp,
			      Epetra_Vector* f,
			      Epetra_MultiVector* JVx,
			      Epetra_MultiVector* fVp)
{
  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  for (unsigned int i=0; i<p.size(); i++) {
    if (p[i] != Teuchos::null)
      for (unsigned int j=0; j<p[i]->size(); j++)
	(*(p[i]))[j].family->setRealValueForAllTypes((*(p[i]))[j].baseValue);
  }

  // Zero out overlapped residual
  Teuchos::RCP<Epetra_Vector> overlapped_ff;
  if (f != NULL) {
    overlapped_ff = overlapped_f;
    overlapped_ff->PutScalar(0.0);
    f->PutScalar(0.0);
  }

  Teuchos::RCP<Epetra_MultiVector> overlapped_JVx;
  if (JVx != NULL) {
    overlapped_JVx = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          JVx->NumVectors()));
    overlapped_JVx->PutScalar(0.0);
    JVx->PutScalar(0.0);
  }
  
  Teuchos::RCP<Epetra_MultiVector> overlapped_fVp;
  if (fVp != NULL) {
    overlapped_fVp = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          fVp->NumVectors()));
    overlapped_fVp->PutScalar(0.0);
    fVp->PutScalar(0.0);
  }

  Teuchos::RCP<Epetra_MultiVector> overlapped_Vx;
  if (Vx != NULL) {
    overlapped_Vx = 
      Teuchos::rcp(new Epetra_MultiVector(*(disc->getOverlapMap()), 
                                          Vx->NumVectors()));
  }

  Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> > vp =
    Teuchos::rcp(Vp, false);
  Teuchos::RCP<ParamVec> params = 
    Teuchos::rcp(deriv_p, false);

  // Create Jacobian init/post op
  Teuchos::RCP<FEApp::TangentOp> op = 
    Teuchos::rcp(new FEApp::TangentOp(alpha, beta, sum_derivs,
                                      overlapped_xdot, 
                                      overlapped_x,
                                      params,
                                      overlapped_Vx,
                                      overlapped_Vx,
                                      vp,
                                      overlapped_ff, 
                                      overlapped_JVx,
                                      overlapped_fVp));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::TangentType> > pde = 
    pdeTM.getAsObject<FEApp::TangentType>();

  // Do global fill
  FEApp::GlobalFill<FEApp::TangentType> globalFill(disc->getMesh(), 
                                                   quad, pde, bc, 
                                                   transient);
  globalFill.computeGlobalFill(*op);

  // Assemble global residual
  if (f != NULL)
    f->Export(*overlapped_f, *exporter, Add);

  // Assemble derivatives
  if (JVx != NULL)
    JVx->Export(*overlapped_JVx, *exporter, Add);
  if (fVp != NULL)
    fVp->Export(*overlapped_fVp, *exporter, Add);
}

void
FEApp::Application::
evaluateResponses(const Epetra_Vector* xdot,
                  const Epetra_Vector& x,
                  const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
                  Epetra_Vector& g)
{
  const Epetra_Comm& comm = x.Map().Comm();
  unsigned int offset = 0;
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Epetra_LocalMap local_response_map(num_responses, 0, comm);

    // Create Epetra_Vector for response function
    Epetra_Vector local_g(local_response_map);

    // Evaluate response function
    responses[i]->evaluateResponses(xdot, x, p, local_g);

    // Copy result into combined result
    for (unsigned int j=0; j<num_responses; j++)
      g[offset+j] = local_g[j];

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::
evaluateResponseTangents(
	   const Epetra_Vector* xdot,
	   const Epetra_Vector& x,
	   const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
	   const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
	   const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dxdot_dp,
	   const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dx_dp,
	   Epetra_Vector* g,
	   const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& gt)
{
  const Epetra_Comm& comm = x.Map().Comm();
  unsigned int offset = 0;
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > local_gt(gt.size());
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Epetra_LocalMap local_response_map(num_responses, 0, comm);

    // Create Epetra_Vectors for response function
    Teuchos::RCP<Epetra_Vector> local_g;
    if (g != NULL)
      local_g = Teuchos::rcp(new Epetra_Vector(local_response_map));
    for (unsigned int j=0; j<gt.size(); j++)
      if (gt[j] != Teuchos::null)
	local_gt[j] = Teuchos::rcp(new Epetra_MultiVector(local_response_map, 
							  gt[j]->NumVectors()));

    // Evaluate response function
    responses[i]->evaluateTangents(xdot, x, p, deriv_p, dxdot_dp, dx_dp, 
				   local_g.get(), local_gt);

    // Copy results into combined result
    for (unsigned int j=0; j<num_responses; j++) {
      if (g != NULL)
        (*g)[offset+j] = (*local_g)[j];
      for (unsigned int l=0; l<gt.size(); l++)
	if (gt[l] != Teuchos::null)
	  for (int k=0; k<gt[l]->NumVectors(); k++)
	    (*gt[l])[k][offset+j] = (*local_gt[l])[k][j];
    }

    // Increment offset in combined result
    offset += num_responses;
  }
}

void
FEApp::Application::
evaluateResponseGradients(
	    const Epetra_Vector* xdot,
	    const Epetra_Vector& x,
	    const Teuchos::Array< Teuchos::RCP<ParamVec> >& p,
	    const Teuchos::Array< Teuchos::RCP<ParamVec> >& deriv_p,
	    Epetra_Vector* g,
	    Epetra_MultiVector* dg_dx,
	    Epetra_MultiVector* dg_dxdot,
	    const Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> >& dg_dp)
{
  const Epetra_Comm& comm = x.Map().Comm();
  unsigned int offset = 0;
  Teuchos::Array< Teuchos::RCP<Epetra_MultiVector> > local_dgdp(dg_dp.size());
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Epetra_LocalMap local_response_map(num_responses, 0, comm);

    // Create Epetra_Vectors for response function
    Teuchos::RCP<Epetra_Vector> local_g;
    if (g != NULL)
      local_g = Teuchos::rcp(new Epetra_Vector(local_response_map));
    Teuchos::RCP<Epetra_MultiVector> local_dgdx;
    if (dg_dx != NULL)
      local_dgdx = Teuchos::rcp(new Epetra_MultiVector(dg_dx->Map(), 
                                                       num_responses));
    Teuchos::RCP<Epetra_MultiVector> local_dgdxdot;
    if (dg_dxdot != NULL)
      local_dgdxdot = Teuchos::rcp(new Epetra_MultiVector(dg_dxdot->Map(), 
                                                          num_responses));

    for (unsigned int j=0; j<dg_dp.size(); j++)
      if (dg_dp[j] != Teuchos::null)
	local_dgdp[j] = Teuchos::rcp(new Epetra_MultiVector(local_response_map, 
							    dg_dp[j]->NumVectors()));

    // Evaluate response function
    responses[i]->evaluateGradients(xdot, x, p, deriv_p, local_g.get(), 
                                    local_dgdx.get(), local_dgdxdot.get(), 
                                    local_dgdp);

    // Copy results into combined result
    for (unsigned int j=0; j<num_responses; j++) {
      if (g != NULL)
        (*g)[offset+j] = (*local_g)[j];
      if (dg_dx != NULL)
        (*dg_dx)(offset+j)->Update(1.0, *((*local_dgdx)(j)), 0.0);
      if (dg_dxdot != NULL)
        (*dg_dxdot)(offset+j)->Update(1.0, *((*local_dgdxdot)(j)), 0.0);
      for (unsigned int l=0; l<dg_dp.size(); l++)
	if (dg_dp[l] != Teuchos::null)
	  for (int k=0; k<dg_dp[l]->NumVectors(); k++)
	    (*dg_dp[l])[k][offset+j] = (*local_dgdp[l])[k][j];
    }

    // Increment offset in combined result
    offset += num_responses;
  }
}

#if SG_ACTIVE
void
FEApp::Application::computeGlobalSGResidual(
			const Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_xdot,
			const Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_x,
			const ParamVec* p,
			const ParamVec* sg_p,
			const Teuchos::Array<SGType>* sg_p_vals,
			Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_f)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGResidual");

  for (int i=0; i<sg_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*sg_overlapped_x)[i].Import(sg_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (transient)
      (*sg_overlapped_xdot)[i].Import((*sg_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    (*sg_overlapped_f)[i].PutScalar(0.0);
    sg_f[i].PutScalar(0.0);

  }

  // Set real parameters
  if (p != NULL) {
    for (unsigned int i=0; i<p->size(); ++i) {
      (*p)[i].family->setRealValueForAllTypes((*p)[i].baseValue);
    }
  }

  // Set SG parameters
  if (sg_p != NULL && sg_p_vals != NULL) {
    for (unsigned int i=0; i<sg_p->size(); ++i) {
      (*sg_p)[i].family->setValue<FEApp::SGResidualType>((*sg_p_vals)[i]);
    }
  }

  // Create residual init/post op
  Teuchos::RCP<FEApp::SGResidualOp> sg_res_fill_op = 
    Teuchos::rcp(new FEApp::SGResidualOp(sg_expansion, 
					 sg_overlapped_xdot, 
					 sg_overlapped_x, 
					 sg_overlapped_f));
    
  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGResidualType> > pde = 
    pdeTM.getAsObject<FEApp::SGResidualType>();

  // Instantiate global fill
  if (sg_res_global_fill == Teuchos::null) {
    std::string method = params->get("SG Method", "AD");
    if (method == "AD") {
      sg_res_global_fill = 
        Teuchos::rcp(new FEApp::GlobalFill<FEApp::SGResidualType>(
	  disc->getMesh(), quad, pde, bc, transient));
    }
    else if (method == "Gauss Quadrature") {
      Teuchos::RCP< FEApp::AbstractPDE<FEApp::ResidualType> > res_pde = 
	pdeTM.getAsObject<FEApp::ResidualType>();
      sg_res_global_fill = 
        Teuchos::rcp(new FEApp::SGGaussQuadResidualGlobalFill(disc->getMesh(), 
                                                              quad, pde, bc, 
                                                              transient,
                                                              sg_basis,
							      sg_quad,
                                                              res_pde,
                                                              sg_p));
    }
  }

  // Do global fill
  sg_res_global_fill->computeGlobalFill(*sg_res_fill_op);

  // Assemble global residual
  for (int i=0; i<sg_f.size(); i++)
    sg_f[i].Export((*sg_overlapped_f)[i], *exporter, Add);
}

void
FEApp::Application::computeGlobalSGJacobian(
			double alpha, double beta,
			const Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_xdot,
			const Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_x,
			const ParamVec* p,
			const ParamVec* sg_p,
			const Teuchos::Array<SGType>* sg_p_vals,
			Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_f,
			Stokhos::VectorOrthogPoly<Epetra_Operator>& sg_jac)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGJacobian");

  for (int i=0; i<sg_x.size(); i++) {

    // Scatter x to the overlapped distrbution
    (*sg_overlapped_x)[i].Import(sg_x[i], *importer, Insert);

    // Scatter xdot to the overlapped distribution
    if (transient)
      (*sg_overlapped_xdot)[i].Import((*sg_xdot)[i], *importer, Insert);

    // Zero out overlapped residual
    if (sg_f != NULL) {
      (*sg_overlapped_f)[i].PutScalar(0.0);
      (*sg_f)[i].PutScalar(0.0);
    }

  }

  // Create, resize and initialize overlapped Jacobians
  if (sg_overlapped_jac == Teuchos::null || 
      sg_overlapped_jac->size() < sg_jac.size())
    sg_overlapped_jac = 
      Teuchos::rcp(new Stokhos::VectorOrthogPoly<Epetra_CrsMatrix>(
		     sg_basis,  *overlapped_jac, sg_jac.size()));
  else if (sg_overlapped_jac->size() > sg_jac.size())
    sg_overlapped_jac->resize(sg_jac.size());
  for (int i=0; i<sg_overlapped_jac->size(); i++)
    (*sg_overlapped_jac)[i].PutScalar(0.0);

  // Set real parameters
  if (p != NULL) {
    for (unsigned int i=0; i<p->size(); ++i) {
      (*p)[i].family->setRealValueForAllTypes((*p)[i].baseValue);
    }
  }

  // Set SG parameters
  if (sg_p != NULL && sg_p_vals != NULL) {
    for (unsigned int i=0; i<sg_p->size(); ++i) {
      (*sg_p)[i].family->setValue<FEApp::SGJacobianType>((*sg_p_vals)[i]);
    }
  }

  // Create Jacobian init/post op
  Teuchos::RCP< Stokhos::VectorOrthogPoly<Epetra_Vector> > sg_overlapped_ff;
  if (sg_f != NULL)
    sg_overlapped_ff = sg_overlapped_f;
  Teuchos::RCP<FEApp::SGJacobianOp> sg_jac_fill_op = 
    Teuchos::rcp(new FEApp::SGJacobianOp(sg_expansion,
					 alpha, beta, 
					 sg_overlapped_xdot, 
					 sg_overlapped_x, 
					 sg_overlapped_ff, 
					 sg_overlapped_jac));

  // Get template PDE instantiation
  Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGJacobianType> > pde = 
    pdeTM.getAsObject<FEApp::SGJacobianType>();

  // Instantiate global fill
  if (sg_jac_global_fill == Teuchos::null) {
    std::string method = params->get("SG Method", "AD");
    if (method == "AD") {
      sg_jac_global_fill = 
	Teuchos::rcp(new FEApp::GlobalFill<FEApp::SGJacobianType>(
	    disc->getMesh(), quad, pde, bc,transient));
    }
    else if (method == "Gauss Quadrature") {
      Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > jac_pde = 
	pdeTM.getAsObject<FEApp::JacobianType>();
      sg_jac_global_fill = 
	Teuchos::rcp(new FEApp::SGGaussQuadJacobianGlobalFill(disc->getMesh(),
							      quad, pde, bc, 
							      transient,
							      sg_basis,
							      sg_quad,
							      jac_pde,
							      sg_p, 
							      alpha, beta));
    }
  }

  // Do global fill
  sg_jac_global_fill->computeGlobalFill(*sg_jac_fill_op);
  
  // Assemble global residual
  if (sg_f != NULL)
    for (int i=0; i<sg_f->size(); i++)
      (*sg_f)[i].Export((*sg_overlapped_f)[i], *exporter, Add);
    
  // Assemble block Jacobians
  Teuchos::RCP<Epetra_CrsMatrix> jac;
  for (int i=0; i<sg_jac.size(); i++) {
    jac = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(sg_jac.getCoeffPtr(i), 
						      true);
    jac->PutScalar(0.0);
    jac->Export((*sg_overlapped_jac)[i], *exporter, Add);
    jac->FillComplete(true);
  }
}

void
FEApp::Application::
evaluateSGResponses(const Stokhos::VectorOrthogPoly<Epetra_Vector>* sg_xdot,
		    const Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_x,
		    const ParamVec* p,
		    const ParamVec* sg_p,
		    const Teuchos::Array<SGType>* sg_p_vals,
		    Stokhos::VectorOrthogPoly<Epetra_Vector>& sg_g)
{
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::evaluateSGResponses");

  const Epetra_Comm& comm = sg_x[0].Map().Comm();
  unsigned int offset = 0;
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis = 
    sg_x.basis();
  for (unsigned int i=0; i<responses.size(); i++) {

    // Create Epetra_Map for response function
    unsigned int num_responses = responses[i]->numResponses();
    Epetra_LocalMap local_response_map(num_responses, 0, comm);

    // Create Epetra_Vector for response function
    Stokhos::VectorOrthogPoly<Epetra_Vector> local_sg_g(basis,
							local_response_map);

    // Evaluate response function
    responses[i]->evaluateSGResponses(sg_xdot, sg_x, p, sg_p, sg_p_vals,
				      local_sg_g);

    // Copy result into combined result
    for (int k=0; k<sg_g.size(); k++)
      for (unsigned int j=0; j<num_responses; j++)
	sg_g[k][offset+j] = local_sg_g[k][j];

    // Increment offset in combined result
    offset += num_responses;
  }
}

#endif
