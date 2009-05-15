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
		   bool is_transient) :
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
                        pdeTM, bc, initial_x);
  typedef FEApp::AbstractPDE_TemplateManager<EvalTypes>::iterator iterator;
  int nqp = quad->numPoints();
  int nn = disc->getNumNodesPerElement();
  for (iterator it = pdeTM.begin(); it != pdeTM.end(); ++it)
    it->init(nqp, nn);

#if SG_ACTIVE
  enable_sg = params->get("Enable Stochastic Galerkin",false);
  if (enable_sg) {
    sg_solver_method = params->get("SG Solver Method", "Fully Assembled");
    sg_basis = params->get< Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > >("Stochastic Galerkin basis");
    sg_quad = params->get< Teuchos::RCP<const Stokhos::Quadrature<int,double> > >("Stochastic Galerkin quadrature");
    Cijk = params->get< Teuchos::RCP<const tp_type> >("Stochastic Galerkin triple product");

    bool makeJacobian = sg_solver_method == "Fully Assembled";
    sg_disc = Teuchos::rcp(new FEApp::BlockDiscretization(comm, disc,
                                                          sg_basis,
                                                          makeJacobian));

    // Create Epetra objects
    sg_initial_x = 
      Teuchos::rcp(new EpetraExt::BlockVector(*(disc->getMap()),
                                              *(sg_disc->getMap())));
    sg_importer = Teuchos::rcp(new Epetra_Import(*(sg_disc->getOverlapMap()), 
                                                 *(sg_disc->getMap())));
    sg_exporter = Teuchos::rcp(new Epetra_Export(*(sg_disc->getOverlapMap()), 
                                                 *(sg_disc->getMap())));
    sg_overlapped_x = 
      Teuchos::rcp(new EpetraExt::BlockVector(*(disc->getOverlapMap()),
                                              *(sg_disc->getOverlapMap())));
    if (transient)
      sg_overlapped_xdot = 
        Teuchos::rcp(new EpetraExt::BlockVector(*(disc->getOverlapMap()),
                                                *(sg_disc->getOverlapMap())));
    sg_overlapped_f = 
      Teuchos::rcp(new EpetraExt::BlockVector(*(disc->getOverlapMap()),
                                              *(sg_disc->getOverlapMap())));
    if (sg_solver_method == "Fully Assembled")
      sg_overlapped_jac = sg_disc->getOverlapJacobian();
    else if (sg_solver_method == "Matrix Free Mean Prec")
      precFactory = params->get<Teuchos::RCP<FEApp::PreconditionerFactory> >("SG Preconditioner Factory");
  }
#endif
}

FEApp::Application::~Application()
{
}

Teuchos::RCP<const Epetra_Map>
FEApp::Application::getMap() const
{
#if SG_ACTIVE
  if (enable_sg)
    return sg_disc->getMap();
  else
#endif
    return disc->getMap();
}

Teuchos::RCP<const Epetra_CrsGraph>
FEApp::Application::getJacobianGraph() const
{
#if SG_ACTIVE
  if (enable_sg)
    return sg_disc->getJacobianGraph();
  else
#endif
    return disc->getJacobianGraph();
}

Teuchos::RCP<const Epetra_Vector>
FEApp::Application::getInitialSolution() const
{
#if SG_ACTIVE
  if (enable_sg)
    return sg_initial_x;
  else
#endif
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
#if SG_ACTIVE
  if (enable_sg) {
    if (sg_solver_method == "Fully Assembled") {
      return 
        Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
                                          *(sg_disc->getJacobianGraph())));
    }
    else if (sg_solver_method == "Matrix Free Mean Prec") {
      unsigned int sz = sg_basis->size();
      Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_Operator> > > jacs = 
        Teuchos::rcp(new std::vector< Teuchos::RCP<Epetra_Operator> >(sz));
      for (unsigned int i=0; i<sz; i++)
        (*jacs)[i] = 
          Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
					    *(disc->getJacobianGraph())));
      return 
        Teuchos::rcp(new Stokhos::MatrixFreeEpetraOp(disc->getMap(),
						     sg_disc->getMap(),
						     sg_basis,
						     Cijk,
						     jacs));
    }
  }
  else
#endif
    return Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
                                             *(disc->getJacobianGraph())));
}

Teuchos::RCP<Epetra_Operator>
FEApp::Application::createPrec() const
{
#if SG_ACTIVE
  if (enable_sg) {
    if (sg_solver_method == "Fully Assembled") {
      return 
        Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
                                          *(sg_disc->getJacobianGraph())));
    }
    else if (sg_solver_method == "Matrix Free Mean Prec") {
      unsigned int sz = sg_basis->size();
      return Teuchos::rcp(new Stokhos::MeanEpetraOp(disc->getMap(),
						    sg_disc->getMap(),
						    sz,
						    Teuchos::null));
    }
  }
  else
#endif
    return Teuchos::rcp(new Epetra_CrsMatrix(::Copy, 
                                             *(disc->getJacobianGraph())));
}

void
FEApp::Application::computeGlobalResidual(
                                      const Epetra_Vector* xdot,
                                      const Epetra_Vector& x,
                                      const ParamVec* p,
                                      Epetra_Vector& f)
{
  if (enable_sg) {
    computeGlobalSGResidual(xdot, x, p, f);
    return;
  }

  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalResidual");

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  if (p != NULL) {
    for (unsigned int i=0; i<p->size(); ++i) {
      (*p)[i].family->setRealValueForAllTypes((*p)[i].baseValue);
    }
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
				      const ParamVec* p,
				      Epetra_Vector* f,
				      Epetra_Operator& jacOp)
{
  if (enable_sg) {
    computeGlobalSGJacobian(alpha, beta, xdot, x, p, f, jacOp);
    return;
  }

  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalJacobian");

  // Cast jacOp to an Epetra_CrsMatrix
  Epetra_CrsMatrix& jac = dynamic_cast<Epetra_CrsMatrix&>(jacOp);

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  if (p != NULL) {
    for (unsigned int i=0; i<p->size(); ++i) {
      (*p)[i].family->setRealValueForAllTypes((*p)[i].baseValue);
    }
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
				      const ParamVec* p,
				      Epetra_Vector* f,
				      Epetra_Operator& jacOp)
{
  if (enable_sg) {
    computeGlobalSGPreconditioner(alpha, beta, xdot, x, p, f, jacOp);
    return;
  }

  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalPreconditioner");

  // Cast jacOp to an Epetra_CrsMatrix
  Epetra_CrsMatrix& jac = dynamic_cast<Epetra_CrsMatrix&>(jacOp);

  // Scatter x to the overlapped distrbution
  overlapped_x->Import(x, *importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    overlapped_xdot->Import(*xdot, *importer, Insert);

  // Set parameters
  if (p != NULL) {
    for (unsigned int i=0; i<p->size(); ++i) {
      (*p)[i].family->setRealValueForAllTypes((*p)[i].baseValue);
    }
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
			      ParamVec* p,
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
    Teuchos::rcp(p, false);

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
FEApp::Application::computeGlobalSGResidual(
				    const Epetra_Vector* sg_xdot,
				    const Epetra_Vector& sg_x,
				    const ParamVec* p,
				    Epetra_Vector& sg_f)
{
#if SG_ACTIVE
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGResidual");

  // Scatter x to the overlapped distrbution
  sg_overlapped_x->Import(sg_x, *sg_importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    sg_overlapped_xdot->Import(*sg_xdot, *sg_importer, Insert);

  // Set parameters
//   if (p != NULL) {
//     for (unsigned int i=0; i<p->size(); ++i) {
//       (*p)[i].family->setRealValueForAllTypes((*p)[i].baseValue);
//     }
//   }

  // Zero out overlapped residual
  sg_overlapped_f->PutScalar(0.0);
  sg_f.PutScalar(0.0);

  std::string method = params->get("SG Method", "AD");
  if (method == "AD") {

    // Create residual init/post op
    if (sg_res_fill_op == Teuchos::null)
      sg_res_fill_op = 
        Teuchos::rcp(new FEApp::SGResidualOp(disc->getOverlapMap(),
                                             sg_basis,
                                             sg_overlapped_xdot, 
                                             sg_overlapped_x, 
                                             sg_overlapped_f));
    else
      sg_res_fill_op->reset(sg_overlapped_xdot, sg_overlapped_x);
    
    // Get template PDE instantiation
    Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGResidualType> > pde = 
      pdeTM.getAsObject<FEApp::SGResidualType>();

    // Instantiate global fill
    if (sg_res_global_fill == Teuchos::null) 
      sg_res_global_fill = 
        Teuchos::rcp(new FEApp::GlobalFill<FEApp::SGResidualType>(disc->getMesh(), quad, pde, bc, transient));

    // Do global fill
    sg_res_global_fill->computeGlobalFill(*sg_res_fill_op);
  }

  else if (method == "Gauss Quadrature") {
    // Create residual init/post op
    if (sg_res_fill_op == Teuchos::null)
      sg_res_fill_op = 
        Teuchos::rcp(new FEApp::SGResidualOp(disc->getOverlapMap(),
                                             sg_basis,
                                             sg_overlapped_xdot, 
                                             sg_overlapped_x, 
                                             sg_overlapped_f));
    else
      sg_res_fill_op->reset(sg_overlapped_xdot, sg_overlapped_x);

    // Get template PDE instantiation
    Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGResidualType> > pde = 
      pdeTM.getAsObject<FEApp::SGResidualType>();
    Teuchos::RCP< FEApp::AbstractPDE<FEApp::ResidualType> > res_pde = 
      pdeTM.getAsObject<FEApp::ResidualType>();
    
    // Instantiate global fill
    if (sg_res_global_fill == Teuchos::null)
      sg_res_global_fill = 
        Teuchos::rcp(new FEApp::SGGaussQuadResidualGlobalFill(disc->getMesh(), 
                                                              quad, pde, bc, 
                                                              transient,
                                                              sg_basis,
							      sg_quad,
                                                              res_pde,
                                                              p));

    // Do global fill
    sg_res_global_fill->computeGlobalFill(*sg_res_fill_op);
  }

  // Assemble global residual
  sg_f.Export(*sg_overlapped_f, *sg_exporter, Add);
#endif
}

void
FEApp::Application::computeGlobalSGJacobian(
				      double alpha, double beta,
				      const Epetra_Vector* sg_xdot,
				      const Epetra_Vector& sg_x,
				      const ParamVec* p,
				      Epetra_Vector* sg_f,
				      Epetra_Operator& sg_jacOp)
{
#if SG_ACTIVE
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGJacobian");

  // Scatter x to the overlapped distrbution
  sg_overlapped_x->Import(sg_x, *sg_importer, Insert);

  // Scatter xdot to the overlapped distribution
  if (transient)
    sg_overlapped_xdot->Import(*sg_xdot, *sg_importer, Insert);

  // Set parameters
//   if (p != NULL) {
//     for (unsigned int i=0; i<p->size(); ++i) {
//       (*p)[i].family->setRealValueForAllTypes((*p)[i].baseValue);
//     }
//   }

  // Zero out overlapped residual
  Teuchos::RCP<EpetraExt::BlockVector> sg_overlapped_ff;
  if (sg_f != NULL) {
    sg_overlapped_ff = sg_overlapped_f;
    sg_overlapped_ff->PutScalar(0.0);
    sg_f->PutScalar(0.0);
  }

  std::string method = params->get("SG Method", "AD");
  if (method == "AD" || method == "Gauss Quadrature") {

    // Create Jacobian init/post op
    Teuchos::RCP<FEApp::AbstractInitPostOp<FEApp::SGJacobianType> > jac_fill_op;
    if (sg_solver_method == "Fully Assembled") {
      sg_overlapped_jac->PutScalar(0.0);
      if (sg_full_jac_fill_op == Teuchos::null)
        sg_full_jac_fill_op = 
          Teuchos::rcp(new FEApp::SGJacobianOp(disc->getOverlapMap(),
                                               disc->getOverlapJacobianGraph(),
                                               sg_basis,
                                               Cijk,
                                               alpha, beta, 
                                               sg_overlapped_xdot, 
                                               sg_overlapped_x, 
                                               sg_overlapped_ff, 
                                               sg_overlapped_jac));
      else
        sg_full_jac_fill_op->reset(alpha, beta, sg_overlapped_xdot, 
                                   sg_overlapped_x);
      jac_fill_op = sg_full_jac_fill_op;
    }
    else if (sg_solver_method == "Matrix Free Mean Prec") {
      if (sg_mf_jac_fill_op == Teuchos::null)
        sg_mf_jac_fill_op = 
          Teuchos::rcp(new FEApp::SGMatrixFreeJacobianOp(
                                            disc->getOverlapMap(),
                                            disc->getOverlapJacobianGraph(),
                                            sg_basis,
                                            alpha, beta, 
                                            sg_overlapped_xdot, 
                                            sg_overlapped_x, 
                                            sg_overlapped_ff));
      else
        sg_mf_jac_fill_op->reset(alpha, beta, sg_overlapped_xdot, 
                                 sg_overlapped_x);
      jac_fill_op = sg_mf_jac_fill_op;
    }

    // Get template PDE instantiation
    Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGJacobianType> > pde = 
      pdeTM.getAsObject<FEApp::SGJacobianType>();

    // Do global fill
    if (sg_jac_global_fill == Teuchos::null) {
      if (method == "AD") {
	if (sg_jac_global_fill == Teuchos::null)
	  sg_jac_global_fill = 
	    Teuchos::rcp(new FEApp::GlobalFill<FEApp::SGJacobianType>(disc->getMesh(), quad, pde, bc,transient));
      }
      else if (method == "Gauss Quadrature") {
        Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> > jac_pde = 
          pdeTM.getAsObject<FEApp::JacobianType>();
	if (sg_jac_global_fill == Teuchos::null)
	  sg_jac_global_fill = 
	    Teuchos::rcp(new FEApp::SGGaussQuadJacobianGlobalFill(
                                                         disc->getMesh(), 
                                                         quad, pde, bc, 
                                                         transient,
                                                         sg_basis,
							 sg_quad,
                                                         jac_pde,
                                                         p, alpha, beta));
      }
    }
    sg_jac_global_fill->computeGlobalFill(*jac_fill_op);
  }
  
  // Assemble global residual
  if (sg_f != NULL)
    sg_f->Export(*sg_overlapped_f, *sg_exporter, Add);

  
  if (sg_solver_method == "Fully Assembled") {
    // Cast sg_jacOp to an Epetra_CrsMatrix
    Epetra_CrsMatrix& sg_jac = dynamic_cast<Epetra_CrsMatrix&>(sg_jacOp);

    // Assemble global Jacobian
    sg_jac.PutScalar(0.0);
    sg_jac.Export(*sg_overlapped_jac, *sg_exporter, Add);
    sg_jac.FillComplete(true);
  }
  else if (sg_solver_method == "Matrix Free Mean Prec") {
    // Cast sg_jacOp to an Stokhos::MatrixFreeEpetraOp
    Stokhos::MatrixFreeEpetraOp& sg_jac = 
      dynamic_cast<Stokhos::MatrixFreeEpetraOp&>(sg_jacOp);
    
    // Assemble block Jacobians
    std::vector< Teuchos::RCP<Epetra_Operator> > jacs = 
      sg_jac.getOperatorBlocks();
    std::vector< Teuchos::RCP<Epetra_CrsMatrix> > ov_jacs = 
      sg_mf_jac_fill_op->getJacobianBlocks();
    Teuchos::RCP<Epetra_CrsMatrix> jac;
    for (unsigned int i=0; i<jacs.size(); i++) {
      jac = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(jacs[i]);
      jac->PutScalar(0.0);
      jac->Export(*ov_jacs[i], *exporter, Add);
      jac->FillComplete(true);
    }
    mean_jac = jacs[0];
  }
#endif
}

void
FEApp::Application::computeGlobalSGPreconditioner(
				      double alpha, double beta,
				      const Epetra_Vector* sg_xdot,
				      const Epetra_Vector& sg_x,
				      const ParamVec* p,
				      Epetra_Vector* sg_f,
				      Epetra_Operator& sg_jacOp)
{
#if SG_ACTIVE
  TEUCHOS_FUNC_TIME_MONITOR("FEApp::Application::computeGlobalSGPreconditioner");
  
  // We are assuming we have a valid and up-to-date fill
  // from computeGlobalSGJacobian()

  // Assemble global residual
  if (sg_f != NULL)
    sg_f->Export(*sg_overlapped_f, *sg_exporter, Add);
  
  if (sg_solver_method == "Fully Assembled") {
    // Cast sg_jacOp to an Epetra_CrsMatrix
    Epetra_CrsMatrix& sg_jac = dynamic_cast<Epetra_CrsMatrix&>(sg_jacOp);

    // Assemble global Jacobian
    sg_jac.PutScalar(0.0);
    sg_jac.Export(*sg_overlapped_jac, *sg_exporter, Add);
    sg_jac.FillComplete(true);
  }
  else if (sg_solver_method == "Matrix Free Mean Prec") {
    // Cast sg_jacOp to an FEApp::SGMeanPrecOp
    Stokhos::MeanEpetraOp& sg_jac = 
      dynamic_cast<Stokhos::MeanEpetraOp&>(sg_jacOp);

    // Compute preconditioner
    Teuchos::RCP<Epetra_Operator> prec = precFactory->compute(mean_jac);    
    sg_jac.setMeanOperator(prec);
  }
#endif
}


