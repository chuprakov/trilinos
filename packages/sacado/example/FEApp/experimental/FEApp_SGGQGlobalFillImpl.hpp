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
FEApp::SGGQGlobalFill<EvalT>::
SGGQGlobalFill(
      const Teuchos::RCP<const FEApp::Mesh>& elementMesh,
      const Teuchos::RCP<const FEApp::AbstractQuadrature>& quadRule,
      const Teuchos::RCP< FEApp::AbstractPDE<EvalT> >& pdeEquations,
      const std::vector< Teuchos::RCP<FEApp::NodeBC> >& nodeBCs,
      bool is_transient,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sgBasis,
      const Teuchos::RCP<const Stokhos::Quadrature<int,double> >& sgQuad,
      const Teuchos::RCP<const ParamVec>& pvec):
  GlobalFill<EvalT>(elementMesh, quadRule, pdeEquations, nodeBCs, 
                    is_transient),
  sg_basis(sgBasis),
  sg_quad(sgQuad),
  p(pvec)
{
}

template <typename EvalT>
FEApp::SGGQGlobalFill<EvalT>::
~SGGQGlobalFill()
{
}

template <typename EvalT>
void
FEApp::SGGQGlobalFill<EvalT>::
computeGlobalFill(FEApp::AbstractInitPostOp<EvalT>& initPostOp)
{
  const std::vector< std::vector<double> >& quad_points = 
    sg_quad->getQuadPoints();
  unsigned int nqp = quad_points.size();

  // Loop over SG Quad points
  for (unsigned int qp=0; qp<nqp; qp++) {

    // Evaluate parameters
    for (unsigned int i=0; i<p->size(); i++)
      (*p)[i].family->template setValue<EvalT>(quad_points[qp][i]);

    // Set quad point index
    initPostOp.setQuadPointIndex(qp);

    // Loop over elements
    Teuchos::RCP<const FEApp::AbstractElement> e;
    for (FEApp::Mesh::const_iterator eit=this->mesh->begin(); 
         eit!=this->mesh->end(); ++eit) {
      e = *eit;

      // Zero out residual
      for (unsigned int i=0; i<this->ndof; i++)
        this->elem_f[i] = 0.0;
      
      // Initialize element solution
      initPostOp.elementInit(*e, this->neqn, this->elem_xdot, this->elem_x);
      
      // Compute element residual
      this->pde->evaluateElementResidual(*(this->quad), *e, this->elem_xdot, 
                                         this->elem_x, this->elem_f);
      
      // Post-process element residual
      initPostOp.elementPost(*e, this->neqn, this->elem_f);
      
    }

  }


//   // Loop over SG Quad points
//   for (unsigned int qp=0; qp<nqp; qp++) {

//     // Evaluate parameters
//     for (unsigned int i=0; i<p->size(); i++)
//       (*p)[i].family->setValue<EvalT>(quad_points[qp][i]);

    // Set quad point index
//     initPostOp.setQuadPointIndex(qp);
    initPostOp.setQuadPointIndex(0);

    // Loop over boundary conditions
    for (std::size_t i=0; i<this->bc.size(); i++) {

      if (this->bc[i]->isOwned() || this->bc[i]->isShared()) {
        
        // Zero out node residual
        for (unsigned int j=0; j<this->neqn; j++)
          this->node_f[j] = 0.0;

        // Initialize node solution
        initPostOp.nodeInit(*(this->bc[i]), this->neqn, this->node_xdot, 
                            this->node_x);

        // Compute node residual
        this->bc[i]->template getStrategy<EvalT>()->evaluateResidual(
                                                            this->node_xdot, 
                                                            this->node_x, 
                                                            this->node_f);

        // Post-process node residual
        initPostOp.nodePost(*this->bc[i], this->neqn, this->node_f);

      }

    }
    
//   }

  // Finalize fill
  initPostOp.finalizeFill();

}
