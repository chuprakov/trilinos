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

template <typename EvalT>
FEApp::GlobalFill<EvalT>::
GlobalFill(
      const Teuchos::RCP<const FEApp::Mesh>& elementMesh,
      const Teuchos::RCP<const FEApp::AbstractQuadrature>& quadRule,
      const Teuchos::RCP< FEApp::AbstractPDE<EvalT> >& pdeEquations,
      const std::vector< Teuchos::RCP<FEApp::NodeBC> >& nodeBCs,
      bool is_transient):
  mesh(elementMesh),
  quad(quadRule),
  pde(pdeEquations),
  bc(nodeBCs),
  transient(is_transient),
  nnode(0),
  neqn(pde->numEquations()),
  ndof(0),
  elem_x(),
  elem_xdot(NULL),
  elem_f(),
  node_x(),
  node_xdot(NULL),
  node_f()
{
  Teuchos::RCP<const FEApp::AbstractElement> e0 = *(mesh->begin());
  nnode = e0->numNodes();
  ndof = nnode*neqn;
  elem_x.resize(ndof);
  if (transient)
    elem_xdot = new std::vector<ScalarT>(ndof);
  elem_f.resize(ndof);

  node_x.resize(neqn);
  if (transient)
    node_xdot = new std::vector<ScalarT>(neqn);
  node_f.resize(neqn);
}

template <typename EvalT>
FEApp::GlobalFill<EvalT>::
~GlobalFill()
{
  if (transient) {
    delete elem_xdot;
    delete node_xdot;
  }
}

template <typename EvalT>
void
FEApp::GlobalFill<EvalT>::
computeGlobalFill(FEApp::AbstractInitPostOp<EvalT>& initPostOp)
{
  // Loop over elements
  Teuchos::RCP<const FEApp::AbstractElement> e;
  for (FEApp::Mesh::const_iterator eit=mesh->begin(); eit!=mesh->end(); ++eit){
    e = *eit;

    // Zero out element residual
    for (unsigned int i=0; i<ndof; i++)
      elem_f[i] = 0.0;

    // Initialize element solution
    initPostOp.elementInit(*e, neqn, elem_xdot, elem_x);

    // Compute element residual
    pde->evaluateElementResidual(*quad, *e, elem_xdot, elem_x, elem_f);

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
      bc[i]->getStrategy<EvalT>()->evaluateResidual(node_xdot, node_x, 
                                                    node_f);

      // Post-process node residual
      initPostOp.nodePost(*bc[i], neqn, node_f);

    }
    
  }

  // Finalize fill
  initPostOp.finalizeFill();

}
