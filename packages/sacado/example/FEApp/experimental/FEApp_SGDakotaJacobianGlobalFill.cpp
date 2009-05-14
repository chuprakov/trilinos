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

#include "FEApp_SGDakotaJacobianGlobalFill.hpp"

#if SG_ACTIVE

#ifdef HAVE_DAKOTA

#include "DakotaInterface.H"
#include "DakotaIterator.H"
#include "NonDExpansion.H"
#include "DataFitSurrModel.H"

FEApp::SGDakotaJacobianGlobalFill::
SGDakotaJacobianGlobalFill(
      const Teuchos::RCP<const FEApp::Mesh>& elementMesh,
      const Teuchos::RCP<const FEApp::AbstractQuadrature>& quadRule,
      const Teuchos::RCP< FEApp::AbstractPDE<FEApp::SGJacobianType> >& pdeEquations,
      const std::vector< Teuchos::RCP<FEApp::NodeBC> >& nodeBCs,
      bool is_transient,
      const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sgBasis,
      const Teuchos::RCP< FEApp::AbstractPDE<FEApp::JacobianType> >& jacPDEEquations,
      const ParamVec* pvec,
      double alpha,
      double beta):
  GlobalFill<SGJacobianType>(elementMesh, quadRule, pdeEquations, nodeBCs,
                             is_transient),
  sg_basis(sgBasis),
  jacPDE(jacPDEEquations),
  p(pvec),
  sg_size(sg_basis->size()),
  parallel_lib(),
  problem_db(parallel_lib)
{
  const char* dakota_in = "dakota_jac.in";
  const char* dakota_out = "dakota_jac.out";
  const char* dakota_err = "dakota_jac.err";
  const char* dakota_restart_out = "dakota_jac.restart";
  parallel_lib.specify_outputs_restart(dakota_out, dakota_err, NULL,
                                       dakota_restart_out, 0);
  problem_db.manage_inputs(dakota_in);
  problem_db.check_input();
  selected_strategy = Dakota::Strategy(problem_db);
  Dakota::Model& first_model = *(problem_db.model_list().begin());
  Dakota::Interface& interface = first_model.interface();
  elemInterface = 
    new FEApp::DakotaElementJacobianInterface(problem_db, jacPDE, *p, 
                                              quad, ndof,alpha,beta);
  interface.assign_rep(elemInterface, false);
}

FEApp::SGDakotaJacobianGlobalFill::
~SGDakotaJacobianGlobalFill()
{
}

void
FEApp::SGDakotaJacobianGlobalFill::
computeGlobalFill(FEApp::AbstractInitPostOp<FEApp::SGJacobianType>& initPostOp)
{
  // Loop over elements
  Teuchos::RCP<const FEApp::AbstractElement> e;
  for (FEApp::Mesh::const_iterator eit=mesh->begin(); eit!=mesh->end(); ++eit){
    e = *eit;

    // Initialize element solution
    initPostOp.elementInit(*e, neqn, elem_xdot, elem_x);

    // Reset Dakota element interface
  elemInterface->reset(e, Teuchos::rcp(elem_xdot,false), Teuchos::rcp(&elem_x,false));

    // Run Dakota on this element
    selected_strategy.run_strategy();

    // Extract SG components from Dakota
    // We are assuming Dakota and Sacado order them the same way
    Dakota::Iterator& it = *(problem_db.iterator_list().begin());
    Dakota::NonDExpansion& nondexp = 
      dynamic_cast<Dakota::NonDExpansion&>(*(it.iterator_rep()));
    Dakota::Model& sg_model = nondexp.get_uSpaceModel();
    const Dakota::RealVectorArray& f_coeffs = 
      sg_model.approximation_coefficients();
    for (std::size_t i=0; i<ndof; i++) {
      elem_f[i].resize(ndof);
      for (std::size_t j=0; j<ndof; j++) {
        elem_f[i].fastAccessDx(j).copyForWrite();
        elem_f[i].fastAccessDx(j).resize(sg_size);
        for (std::size_t k=0; k<sg_size; k++)
          elem_f[i].fastAccessDx(j).fastAccessCoeff(k) = 
            f_coeffs[j+i*ndof][k];
      }
    }

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

#endif
