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

#include "FEApp_SGGaussQuadJacobianGlobalFill.hpp"

#if SG_ACTIVE

FEApp::SGGaussQuadJacobianGlobalFill::
SGGaussQuadJacobianGlobalFill(
      const Teuchos::RCP<const FEApp::Mesh>& elementMesh,
      const Teuchos::RCP<const FEApp::AbstractQuadrature>& quadRule,
      const Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGJacobianType> >& pdeEquations,
      const std::vector< Teuchos::RCP<FEApp::NodeBC> >& nodeBCs,
      bool is_transient,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sgBasis,
      const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sgQuad,
      const Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> >& jacPDEEquations,
      const ParamVec* pvec,
      double alpha_, double beta_):
  GlobalFill<SGJacobianType>(elementMesh, quadRule, pdeEquations, nodeBCs,
                             is_transient),
  sg_basis(sgBasis),
  sg_quad(sgQuad),
  jacPDE(jacPDEEquations),
  p(pvec),
  alpha(alpha_),
  beta(beta_),
  quad_points(sg_quad->getQuadPoints()),
  quad_weights(sg_quad->getQuadWeights()),
  quad_values(sg_quad->getBasisAtQuadPoints()),
  norms(sg_basis->norm_squared()),
  sg_size(sg_basis->size()),
  nqp(quad_points.size()),
  x(ndof),
  xdot(NULL),
  f(ndof),
  xqp(ndof*nqp),
  xdotqp(),
  pqp(p->size()*nqp),
  fqp(ndof*ndof*nqp),
  qv(nqp*sg_size),
  sqv(nqp*sg_size),
  sg_x(ndof*sg_size),
  sg_xdot(),
  sg_p(p->size()*sg_size),
  sg_f(ndof*ndof*sg_size)
{
  if (transient) {
    xdot = new std::vector<FadType>(ndof);
    xdotqp.resize(ndof*nqp);
    sg_xdot.resize(ndof*sg_size);
  }
  for (unsigned int i=0; i<ndof; i++) {
    x[i].resize(ndof);
    x[i].fastAccessDx(i) = beta;
    if (transient) {
      (*xdot)[i].resize(ndof);
      (*xdot)[i].fastAccessDx(i) = alpha;
    }
    elem_f[i].resize(ndof);
    for (unsigned int j=0; j<ndof; j++) {
      elem_f[i].fastAccessDx(j).copyForWrite();
      elem_f[i].fastAccessDx(j).resize(sg_size);
    }
  }

  for (unsigned int qp=0; qp<nqp; qp++)
    for (unsigned int i=0; i<sg_size; i++) {
      qv[qp*sg_size+i] = quad_values[qp][i];
      sqv[qp*sg_size+i] = quad_values[qp][i]/norms[i];
    }
}

FEApp::SGGaussQuadJacobianGlobalFill::
~SGGaussQuadJacobianGlobalFill()
{
  if (transient)
    delete xdot;
}

void
FEApp::SGGaussQuadJacobianGlobalFill::
computeGlobalFill(FEApp::AbstractInitPostOp<FEApp::SGJacobianType>& initPostOp)
{
  // Evaluate parameters at quadrature points
  for (unsigned int i=0; i<p->size(); i++) {
    SGType pv = (*p)[i].family->getValue<FEApp::SGResidualType>();
    for (unsigned int j=0; j<sg_size; j++)
      sg_p[i*sg_size+j] = pv.fastAccessCoeff(j);
  }
  blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, nqp, p->size(), sg_size, 1.0, 
	    &qv[0], sg_size, &sg_p[0], sg_size, 0.0, &pqp[0], nqp);

  // Loop over elements
  Teuchos::RCP<const FEApp::AbstractElement> e;
  for (FEApp::Mesh::const_iterator eit=mesh->begin(); eit!=mesh->end(); ++eit){
    e = *eit;

    // Initialize element solution
    initPostOp.elementInit(*e, neqn, elem_xdot, elem_x);

    // Evaluate x, xdot at all quadrature points
    for (unsigned int j=0; j<ndof; j++)
      for (unsigned int k=0; k<sg_size; k++)
	sg_x[j*sg_size+k] = elem_x[j].val().fastAccessCoeff(k);
    blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, nqp, ndof, sg_size, 1.0, 
	      &qv[0], sg_size, &sg_x[0], sg_size, 0.0, &xqp[0], nqp);
      
    if (transient) {
      for (unsigned int j=0; j<ndof; j++)
	for (unsigned int k=0; k<sg_size; k++)
	  sg_xdot[j*sg_size+k] = (*elem_xdot)[j].val().fastAccessCoeff(k);
      blas.GEMM(Teuchos::TRANS, Teuchos::NO_TRANS, nqp, ndof, sg_size, 1.0, 
		&qv[0], sg_size, &sg_xdot[0], sg_size, 0.0, &xdotqp[0], nqp);
    }

    // Loop over SG Quad points
    for (unsigned int qp=0; qp<nqp; qp++) {

      // Evaluate parameters
      for (unsigned int i=0; i<p->size(); i++)
        (*p)[i].family->setValue<FEApp::JacobianType>(pqp[i*nqp+qp]);

      // Get x, xdot at quadrature points
      for (unsigned int i=0; i<ndof; i++) {
	x[i].val() = xqp[i*nqp+qp];
	if (transient)
	  (*xdot)[i].val() = xdotqp[i*nqp+qp];
      }

      // Zero out residual
      for (unsigned int i=0; i<ndof; i++)
        f[i] = 0.0;

      // Compute element residual
      jacPDE->evaluateElementResidual(*quad, *e, xdot, x, f);

      // Scale residual by quadrature weights
      for (unsigned int i=0; i<ndof; i++)
	for (unsigned int j=0; j<ndof; j++)
	  fqp[(i*ndof+j)*nqp+qp] = f[i].fastAccessDx(j)*quad_weights[qp];

    }

    blas.GEMM(Teuchos::NO_TRANS, Teuchos::NO_TRANS, sg_size, ndof*ndof, nqp, 
	      1.0, &sqv[0], sg_size, &fqp[0], nqp, 0.0, &sg_f[0], sg_size);
    for (unsigned int i=0; i<ndof; i++)
      for (unsigned int j=0; j<ndof; j++)
	for (unsigned int k=0; k<sg_size; k++)
	  elem_f[i].fastAccessDx(j).fastAccessCoeff(k) = 
	    sg_f[(i*ndof+j)*sg_size+k];

    // Post-process element residual
    initPostOp.elementPost(*e, neqn, elem_f);

  }

  // Loop over boundary conditions
  for (std::size_t i=0; i<bc.size(); i++) {

    if (bc[i]->isOwned() || bc[i]->isShared()) {

      // Zero out node residual
      for (unsigned int j=0; j<neqn; j++)
        node_f[j] = 0.0;

      // Initialize node solution
      initPostOp.nodeInit(*bc[i], neqn, node_xdot, node_x);

      // Compute node residual
      bc[i]->getStrategy<SGJacobianType>()->evaluateResidual(node_xdot, 
                                                             node_x, 
                                                             node_f);

      // Post-process node residual
      initPostOp.nodePost(*bc[i], neqn, node_f);

    }
    
  }

  // Finalize fill
  initPostOp.finalizeFill();

}

#endif
