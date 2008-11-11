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

#include "FEApp_InitPostOps.hpp"
#include "Teuchos_TestForException.hpp"
//#include "Teuchos_Exceptions.hpp"
#include "Epetra_Map.h"
#include "EpetraExt_MatrixMatrix.h"

FEApp::ResidualOp::ResidualOp(
        const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
	      const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
	      const Teuchos::RCP<Epetra_Vector>& overlapped_f) :
  xdot(overlapped_xdot),
  x(overlapped_x),
  f(overlapped_f)
{
}

FEApp::ResidualOp::~ResidualOp()
{
}

void
FEApp::ResidualOp::elementInit(const FEApp::AbstractElement& e,
                               unsigned int neqn,
                               std::vector<double>* elem_xdot,
                               std::vector<double>& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
    firstDOF = x->Map().LID(node_GID*neqn);
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = (*x)[firstDOF+j];
      if (elem_xdot != NULL)
        (*elem_xdot)[neqn*i+j] = (*xdot)[firstDOF+j];
    }
  }
}

void
FEApp::ResidualOp::elementPost(const FEApp::AbstractElement& e,
                               unsigned int neqn,
                               std::vector<double>& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual into global residual
  int row;
  unsigned int lrow;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
      f->SumIntoGlobalValue(row, 0, elem_f[lrow]);
    }
  }
}

void
FEApp::ResidualOp::nodeInit(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector<double>* node_xdot,
                            std::vector<double>& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
  unsigned int firstDOF = x->Map().LID(node_GID*neqn);

  // Copy node solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = (*x)[firstDOF+j];
    if (node_xdot != NULL)
      (*node_xdot)[j] = (*xdot)[firstDOF+j];
  }
}

void
FEApp::ResidualOp::nodePost(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector<double>& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  // Replace global residual
  for (unsigned int j=0; j<offsets.size(); j++) {
    int row = static_cast<int>(firstDOF + offsets[j]);
    if (bc.isOwned())
      f->ReplaceGlobalValue(row, 0, node_f[offsets[j]]);
    else if (bc.isShared())
      f->ReplaceGlobalValue(row, 0, 0.0);
  }
}

FEApp::JacobianOp::JacobianOp(
             double alpha, double beta,
             const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
             const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
             const Teuchos::RCP<Epetra_Vector>& overlapped_f,
             const Teuchos::RCP<Epetra_CrsMatrix>& overlapped_jac) :
  m_coeff(alpha),
  j_coeff(beta),
  xdot(overlapped_xdot),
  x(overlapped_x),
  f(overlapped_f),
  jac(overlapped_jac)
{
}

FEApp::JacobianOp::~JacobianOp()
{
}

void
FEApp::JacobianOp::elementInit(const FEApp::AbstractElement& e,
                               unsigned int neqn,
                               std::vector< FadType >* elem_xdot,
                               std::vector< FadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of dof
  unsigned int ndof = nnode*neqn;

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
    firstDOF = x->Map().LID(node_GID*neqn);
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = 
        FadType(ndof, (*x)[firstDOF+j]);
      elem_x[neqn*i+j].fastAccessDx(neqn*i+j) = j_coeff;
      if (elem_xdot != NULL) {
        (*elem_xdot)[neqn*i+j] = 
          FadType(ndof, (*xdot)[firstDOF+j]);
        (*elem_xdot)[neqn*i+j].fastAccessDx(neqn*i+j) = m_coeff;
      }
    }
  }
}

void
FEApp::JacobianOp::elementPost(
			    const FEApp::AbstractElement& e,
			    unsigned int neqn,
			    std::vector< FadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual and Jacobian into global residual, Jacobian
  // Loop over nodes in element
  int row, col;
  unsigned int lrow, lcol;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
        f->SumIntoGlobalValue(row, 0, elem_f[lrow].val());
	
      // Check derivative array is nonzero
      if (elem_f[lrow].hasFastAccess()) {

        // Loop over nodes in element
        for (unsigned int node_col=0; node_col<nnode; node_col++){
          
          // Loop over equations per node
          for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
            lcol = neqn*node_col+eq_col;
            
            // Global column
            col = static_cast<int>(e.nodeGID(node_col)*neqn + eq_col);
            
            // Sum Jacobian
            jac->SumIntoGlobalValues(row, 1, 
                                     &(elem_f[lrow].fastAccessDx(lcol)),
                                     &col);
            
          } // column equations
          
        } // column nodes
        
      } // has fast access
      
    } // row equations
    
  } // row node
  
}

void
FEApp::JacobianOp::nodeInit(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector< FadType >* node_xdot,
                            std::vector< FadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
  unsigned int firstDOF = x->Map().LID(node_GID*neqn);

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = FadType(neqn, (*x)[firstDOF+j]);
    node_x[j].fastAccessDx(j) = j_coeff;
    if (node_xdot != NULL) {
      (*node_xdot)[j] = FadType(neqn, (*xdot)[firstDOF+j]);
      (*node_xdot)[j].fastAccessDx(j) = m_coeff;
    }
  }
}

void
FEApp::JacobianOp::nodePost(const FEApp::NodeBC& bc,
                            unsigned int neqn,
                            std::vector< FadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int num_entries;
  double* row_view = 0;
  int row, col;

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
        f->ReplaceGlobalValue(row, 0, node_f[offsets[eq_row]].val());
      else if (bc.isShared())
        f->ReplaceGlobalValue(row, 0, 0.0);
    }

    // Always zero out row (This takes care of the not-owned case)
    if (bc.isOwned() || bc.isShared()) {
      jac->ExtractGlobalRowView(row, num_entries, row_view);
      for (int k=0; k<num_entries; k++)
        row_view[k] = 0.0;
    }
	
    // Check derivative array is nonzero
    if (node_f[offsets[eq_row]].hasFastAccess()) {
	    
      // Loop over equations per node
      for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	      
        // Global column
        col = static_cast<int>(firstDOF + eq_col);

        // Replace Jacobian
        if (bc.isOwned())
          jac->ReplaceGlobalValues(row, 1, 
                                   &(node_f[eq_row].fastAccessDx(eq_col)),
                                   &col);

      } // column equations
	
    } // has fast access
    
  } // row equations

}

FEApp::TangentOp::TangentOp(
        double alpha, double beta, bool sum_derivs_,
        const Teuchos::RCP<const Epetra_Vector>& overlapped_xdot,
        const Teuchos::RCP<const Epetra_Vector>& overlapped_x,
        const Teuchos::RCP<ParamVec>& p,
        const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vx,
        const Teuchos::RCP<const Epetra_MultiVector>& overlapped_Vxdot,
        const Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,double> >& Vp_,
        const Teuchos::RCP<Epetra_Vector>& overlapped_f,
        const Teuchos::RCP<Epetra_MultiVector>& overlapped_JV,
        const Teuchos::RCP<Epetra_MultiVector>& overlapped_fp) :
  m_coeff(alpha),
  j_coeff(beta),
  sum_derivs(sum_derivs_),
  xdot(overlapped_xdot),
  x(overlapped_x),
  params(p),
  Vx(overlapped_Vx),
  Vxdot(overlapped_Vxdot),
  Vp(Vp_),
  f(overlapped_f),
  JV(overlapped_JV),
  fp(overlapped_fp),
  num_cols_x(0),
  num_cols_p(0),
  num_cols_tot(0),
  param_offset(0)
{
  if (Vx != Teuchos::null)
    num_cols_x = Vx->NumVectors();
  else if (Vxdot != Teuchos::null)
    num_cols_x = Vxdot->NumVectors();
  if (params != Teuchos::null) {
    if (Vp != Teuchos::null)
      num_cols_p = Vp->numCols();
    else
      num_cols_p = params->size();
  }
  if (sum_derivs) {
    num_cols_tot = num_cols_x;
    param_offset = 0;
  }
  else {
    num_cols_tot = num_cols_x + num_cols_p;
    param_offset = num_cols_x;
  }

  TEST_FOR_EXCEPTION(sum_derivs && (num_cols_x != 0) && (num_cols_p != 0) && 
                     (num_cols_x != num_cols_p),
                     std::logic_error,
                     "Seed matrices Vx and Vp must have the same number " << 
                     " of columns when sum_derivs is true and both are "
                     << "non-null!" << std::endl);
}

FEApp::TangentOp::~TangentOp()
{
}

void
FEApp::TangentOp::elementInit(
			 const FEApp::AbstractElement& e,
			 unsigned int neqn,
			 std::vector< FadType >* elem_xdot,
			 std::vector< FadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
    firstDOF = x->Map().LID(node_GID*neqn);

    for (unsigned int j=0; j<neqn; j++) {
      if (Vx != Teuchos::null && j_coeff != 0.0) {
        elem_x[neqn*i+j] = 
          FadType(num_cols_tot, (*x)[firstDOF+j]);
        for (int k=0; k<num_cols_x; k++)
          elem_x[neqn*i+j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
      }
      else
        elem_x[neqn*i+j] = 
          FadType((*x)[firstDOF+j]);
      
      if (elem_xdot != NULL) {
        if (Vxdot != Teuchos::null && m_coeff != 0.0) {
          (*elem_xdot)[neqn*i+j] = 
            FadType(num_cols_tot, (*xdot)[firstDOF+j]);
          for (int k=0; k<num_cols_x; k++)
            (*elem_xdot)[neqn*i+j].fastAccessDx(k) = 
              m_coeff*(*Vxdot)[k][firstDOF+j];
        }
        else
          (*elem_xdot)[neqn*i+j] = 
            FadType((*xdot)[firstDOF+j]);
      }   
    }
  }

  if (params != Teuchos::null) {
    FadType p;
    for (unsigned int i=0; i<params->size(); i++) {
      p = FadType(num_cols_tot, (*params)[i].baseValue);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::TangentType>(p);
    }
  }
}

void
FEApp::TangentOp::elementPost(const FEApp::AbstractElement& e,
                              unsigned int neqn,
                              std::vector< FadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  int row;
  unsigned int lrow;
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (f != Teuchos::null)
        f->SumIntoGlobalValue(row, 0, elem_f[lrow].val());
	
      // Extract derivative components for each column
      if (JV != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
          JV->SumIntoGlobalValue(row, col, elem_f[lrow].dx(col));
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
          fp->SumIntoGlobalValue(row, col, elem_f[lrow].dx(col+param_offset));
      
    } // row equations
    
  } // row node

}

void
FEApp::TangentOp::nodeInit(const FEApp::NodeBC& bc,
                           unsigned int neqn,
                           std::vector< FadType >* node_xdot,
                           std::vector< FadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
  unsigned int firstDOF = x->Map().LID(node_GID*neqn);

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    if (Vx != Teuchos::null && j_coeff != 0.0) {
      node_x[j] = FadType(num_cols_tot, (*x)[firstDOF+j]);
      for (int k=0; k < num_cols_x; k++)
        node_x[j].fastAccessDx(k) = j_coeff*(*Vx)[k][firstDOF+j];
    }
    else
      node_x[j] = FadType((*x)[firstDOF+j]);

    if (node_xdot != NULL) {
      if (Vxdot != Teuchos::null && m_coeff != 0.0) {
        (*node_xdot)[j] = 
          FadType(num_cols_tot, (*xdot)[firstDOF+j]);
        for (int k=0; k < num_cols_x; k++)
          (*node_xdot)[j].fastAccessDx(k) = m_coeff*(*Vxdot)[k][firstDOF+j];
      }
      else
        (*node_xdot)[j] = 
          FadType((*xdot)[firstDOF+j]);
    }
  }
  
  if (params != Teuchos::null) {
    FadType p;
    for (unsigned int i=0; i<params->size(); i++) {
      p = FadType(num_cols_tot, (*params)[i].baseValue);
      if (Vp != Teuchos::null) 
        for (int k=0; k<num_cols_p; k++)
          p.fastAccessDx(param_offset+k) = (*Vp)(i,k);
      else
        p.fastAccessDx(param_offset+i) = 1.0;
      (*params)[i].family->setValue<FEApp::TangentType>(p);
    }
  }
}

void
FEApp::TangentOp::nodePost(const FEApp::NodeBC& bc,
			    unsigned int neqn,
			    std::vector< FadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int row;

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (f != Teuchos::null) {
      if (bc.isOwned())
        f->ReplaceGlobalValue(row, 0, node_f[offsets[eq_row]].val());
      else if (bc.isShared())
        f->ReplaceGlobalValue(row, 0, 0.0);
    }
    
    if (bc.isOwned()) {
      
      // Extract derivative components for each column
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
          JV->ReplaceGlobalValue(row, col, node_f[offsets[eq_row]].dx(col));
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
          fp->ReplaceGlobalValue(row, col, 
                                 node_f[offsets[eq_row]].dx(col+param_offset));
    }
    else if (bc.isShared()) {
      if (JV != Teuchos::null && Vx != Teuchos::null)
        for (int col=0; col<num_cols_x; col++)
          JV->ReplaceGlobalValue(row, col, 0.0);
      
      if (fp != Teuchos::null)
        for (int col=0; col<num_cols_p; col++)
          fp->ReplaceGlobalValue(row, col, 0.0);
    }
    
  } // row equations

}

#if SG_ACTIVE

FEApp::SGResidualOp::SGResidualOp(
	      const Teuchos::RCP<const Epetra_Map>& base_map,
        const Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> >& sg_basis_,
        const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
        const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x,
        const Teuchos::RCP<EpetraExt::BlockVector>& sg_overlapped_f) :
  nblock(sg_basis_->size()),
  map(base_map),
  sg_basis(sg_basis_),
  sg_xdot(sg_overlapped_xdot),
  sg_x(sg_overlapped_x),
  sg_f(sg_overlapped_f),
  xdot(sg_xdot != Teuchos::null ? nblock : 0),
  x(nblock)
{
  for (unsigned int block=0; block<nblock; block++) {
    // Allocate vectors for blocks
    if (sg_xdot != Teuchos::null)
      xdot[block] = Teuchos::rcp(new Epetra_Vector(*map));
    x[block] = Teuchos::rcp(new Epetra_Vector(*map));

    // Extract blocks
    sg_x->ExtractBlockValues(*x[block], block);
    if (sg_xdot != Teuchos::null)
      sg_xdot->ExtractBlockValues(*xdot[block], block);
  }
}

FEApp::SGResidualOp::~SGResidualOp()
{
}

void
FEApp::SGResidualOp::reset(
	  const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
	  const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x)
{
  sg_xdot = sg_overlapped_xdot;
  sg_x = sg_overlapped_x;

  // Extract blocks
  for (unsigned int block=0; block<nblock; block++) {
    sg_x->ExtractBlockValues(*x[block], block);
    if (sg_xdot != Teuchos::null)
      sg_xdot->ExtractBlockValues(*xdot[block], block);
  }
}

void
FEApp::SGResidualOp::elementInit(const FEApp::AbstractElement& e,
			       unsigned int neqn,
			       std::vector<SGType>* elem_xdot,
			       std::vector<SGType>& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Allocate appropriate sizes for coefficients
  for (unsigned int i=0; i<nnode; i++) {
    for (unsigned int j=0; j<neqn; j++) {
      if (elem_x[neqn*i+j].size() != nblock)
        elem_x[neqn*i+j].resize(nblock);
      if (elem_xdot != NULL && (*elem_xdot)[neqn*i+j].size() != nblock)
        (*elem_xdot)[neqn*i+j].resize(nblock);
      elem_x[neqn*i+j].copyForWrite();
      if (elem_xdot != NULL)
        (*elem_xdot)[neqn*i+j].copyForWrite();
    }
  }

  for (unsigned int block=0; block<nblock; block++) {

    // Copy element solution
    for (unsigned int i=0; i<nnode; i++) {
      node_GID = e.nodeGID(i);
      firstDOF = map->LID(node_GID*neqn);
      for (unsigned int j=0; j<neqn; j++) {
        elem_x[neqn*i+j].fastAccessCoeff(block) = (*x[block])[firstDOF+j];
        if (elem_xdot != NULL)
          (*elem_xdot)[neqn*i+j].fastAccessCoeff(block) = 
            (*xdot[block])[firstDOF+j];
      }
    }
    
  }
}

void
FEApp::SGResidualOp::elementPost(const FEApp::AbstractElement& e,
                                 unsigned int neqn,
                                 std::vector<SGType>& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of SG blocks
  unsigned int nblock = sg_basis->size();

  int row;
  unsigned int lrow;
  double c;

  for (unsigned int block=0; block<nblock; block++) {

    // Sum element residual into global residual
    for (unsigned int node_row=0; node_row<nnode; node_row++) {
      for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
        lrow = neqn*node_row+eq_row;
        row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);
        c =  elem_f[lrow].coeff(block);
        sg_f->BlockSumIntoGlobalValues(1, &c, &row, block);
      }
    }
    
  }
}

void
FEApp::SGResidualOp::nodeInit(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector<SGType>* node_xdot,
                              std::vector<SGType>& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
  unsigned int firstDOF = map->LID(node_GID*neqn);

  // Allocate appropriate sizes for coefficients
  for (unsigned int j=0; j<neqn; j++) {
    if (node_x[j].size() != nblock)
      node_x[j].resize(nblock);
    if (node_xdot != NULL && node_x[j].size() != nblock)
      (*node_xdot)[j].resize(nblock);
    node_x[j].copyForWrite();
    if (node_xdot != NULL)
      (*node_xdot)[j].copyForWrite();
  }

  for (unsigned int block=0; block<nblock; block++) {

    // Copy node solution
    for (unsigned int j=0; j<neqn; j++) {
      node_x[j].fastAccessCoeff(block) = (*x[block])[firstDOF+j];
      if (node_xdot != NULL)
        (*node_xdot)[j].fastAccessCoeff(block) = (*xdot[block])[firstDOF+j];
    }

  }
}

void
FEApp::SGResidualOp::nodePost(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector<SGType>& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  double c;
  double zero = 0.0;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  for (unsigned int block=0; block<nblock; block++) {

    // Replace global residual
    for (unsigned int j=0; j<offsets.size(); j++) {
      int row = static_cast<int>(firstDOF + offsets[j]);
      if (bc.isOwned()) {
        c = node_f[offsets[j]].coeff(block);
        sg_f->BlockReplaceGlobalValues(1, &c, &row, block);
      }
      else if (bc.isShared())
        sg_f->BlockReplaceGlobalValues(1, &zero, &row, block);
    }

  }
}

void
FEApp::SGResidualOp::finalizeFill()
{
//   // Load Epetra vectors into block vector
//   for (unsigned int block=0; block<nblock; block++) 
//     sg_f->LoadBlockValues(*f[block], block);
}

FEApp::SGJacobianOp::SGJacobianOp(
       const Teuchos::RCP<const Epetra_Map>& base_map,
       const Teuchos::RCP<const Epetra_CrsGraph>& base_graph,
       const Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> >& sg_basis_,
       const Teuchos::RCP<const FEApp::SGJacobianOp::tp_type>& Cijk_,
       double alpha, double beta,
       const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
       const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x,
       const Teuchos::RCP<EpetraExt::BlockVector>& sg_overlapped_f,
       const Teuchos::RCP<EpetraExt::BlockCrsMatrix>& sg_overlapped_jac) :
  nblock(sg_basis_->size()),
  map(base_map),
  graph(base_graph),
  sg_basis(sg_basis_),
  Cijk(Cijk_),
  m_coeff(alpha),
  j_coeff(beta),
  sg_xdot(sg_overlapped_xdot),
  sg_x(sg_overlapped_x),
  sg_f(sg_overlapped_f),
  sg_jac(sg_overlapped_jac),
  xdot(sg_xdot != Teuchos::null ? nblock : 0),
  x(nblock)
{
  for (unsigned int block=0; block<nblock; block++) {
    // Allocate vectors for blocks
    if (sg_xdot != Teuchos::null)
      xdot[block] = Teuchos::rcp(new Epetra_Vector(*map));
    x[block] = Teuchos::rcp(new Epetra_Vector(*map));

    // Extract blocks
    sg_x->ExtractBlockValues(*x[block], block);
    if (sg_xdot != Teuchos::null)
      sg_xdot->ExtractBlockValues(*xdot[block], block);
  }
}

FEApp::SGJacobianOp::~SGJacobianOp()
{
}

void
FEApp::SGJacobianOp::reset(
	  double alpha, double beta,
	  const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
	  const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x)
{
  m_coeff = alpha;
  j_coeff = beta;
  sg_xdot = sg_overlapped_xdot;
  sg_x = sg_overlapped_x;

  // Extract blocks
  for (unsigned int block=0; block<nblock; block++) {
    sg_x->ExtractBlockValues(*x[block], block);
    if (sg_xdot != Teuchos::null)
      sg_xdot->ExtractBlockValues(*xdot[block], block);
  }
}

void
FEApp::SGJacobianOp::elementInit(const FEApp::AbstractElement& e,
				 unsigned int neqn,
				 std::vector< SGFadType >* elem_xdot,
				 std::vector< SGFadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of dof
  unsigned int ndof = nnode*neqn;

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
    firstDOF = map->LID(node_GID*neqn);
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = SGFadType(ndof, 0.0);
      elem_x[neqn*i+j].fastAccessDx(neqn*i+j) = j_coeff;
      elem_x[neqn*i+j].val().resize(nblock);
      elem_x[neqn*i+j].val().copyForWrite();
      for (unsigned int block=0; block<nblock; block++)
        elem_x[neqn*i+j].val().fastAccessCoeff(block) = 
          (*x[block])[firstDOF+j];
      if (elem_xdot != NULL) {
        (*elem_xdot)[neqn*i+j] = SGFadType(ndof, 0.0);
        (*elem_xdot)[neqn*i+j].fastAccessDx(neqn*i+j) = m_coeff;
        (*elem_xdot)[neqn*i+j].val().resize(nblock);
        (*elem_xdot)[neqn*i+j].val().copyForWrite();
        for (unsigned int block=0; block<nblock; block++)
          (*elem_xdot)[neqn*i+j].val().fastAccessCoeff(block) = 
            (*xdot[block])[firstDOF+j];
      }
      
    }
  }

}

void
FEApp::SGJacobianOp::elementPost(const FEApp::AbstractElement& e,
                                 unsigned int neqn,
                                 std::vector< SGFadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual and Jacobian into global residual, Jacobian
  // Loop over nodes in element
  int row, col;
  unsigned int lrow, lcol, nl, i, j;
  double v1, v2, cijk;
  
  for (unsigned int node_row=0; node_row<nnode; node_row++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;

      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (sg_f != Teuchos::null)
        for (unsigned int k=0; k<nblock; k++) {
          v1 = elem_f[lrow].val().coeff(k);
          sg_f->BlockSumIntoGlobalValues(1, &v1, &row, k);
        }
      
      // Check derivative array is nonzero
      if (elem_f[lrow].hasFastAccess()) {
        
        // Loop over nodes in element
        for (unsigned int node_col=0; node_col<nnode; node_col++){
          
          // Loop over equations per node
          for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
            lcol = neqn*node_col+eq_col;
            
            // Global column
            col = static_cast<int>(e.nodeGID(node_col)*neqn + eq_col);
            
            // Sum Jacobian
            for (unsigned int k=0; k<nblock; k++) {
              v1 = elem_f[lrow].fastAccessDx(lcol).coeff(k);
              nl = Cijk->num_values(k);
              for (unsigned int l=0; l<nl; l++) {
                Cijk->triple_value(k,l,i,j,cijk); 
                v2 = v1 * cijk / Cijk->norm_squared(i);
                sg_jac->BlockSumIntoGlobalValues(row, 1, &v2, &col, i, j);
              }
            } // SG blocks
            
          } // column equations
          
        } // column nodes
        
      } // has fast access
      
    } // row equations
      
  } // row node

}

void
FEApp::SGJacobianOp::nodeInit(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector< SGFadType >* node_xdot,
                              std::vector< SGFadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
  unsigned int firstDOF = map->LID(node_GID*neqn);

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = SGFadType(neqn, 0.0);
    node_x[j].fastAccessDx(j) = j_coeff;
    node_x[j].val().resize(nblock);
    node_x[j].val().copyForWrite();
    for (unsigned int block=0; block<nblock; block++)
      node_x[j].val().fastAccessCoeff(block) = (*x[block])[firstDOF+j];
    if (node_xdot != NULL) {
      (*node_xdot)[j] = SGFadType(neqn, 0.0);
      (*node_xdot)[j].fastAccessDx(j) = m_coeff;
      (*node_xdot)[j].val().resize(nblock);
      (*node_xdot)[j].val().copyForWrite();
      for (unsigned int block=0; block<nblock; block++)
        (*node_xdot)[j].val().fastAccessCoeff(block) = 
          (*xdot[block])[firstDOF+j];
    }
  }
}

void
FEApp::SGJacobianOp::nodePost(const FEApp::NodeBC& bc,
                              unsigned int neqn,
                              std::vector< SGFadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int num_entries;
  double* row_view = 0;
  int row, col;
  unsigned int nl, i, j;
  double v1, v2, cijk, zero;
  zero = 0.0;

  // Loop over equations per node
  for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

    // Global row
    row = static_cast<int>(firstDOF + offsets[eq_row]);
    
    // Replace residual
    if (sg_f != Teuchos::null) {
      if (bc.isOwned()) {
        for (unsigned int k=0; k<nblock; k++) {
          v1 = node_f[offsets[eq_row]].val().coeff(k);
          sg_f->BlockReplaceGlobalValues(1, &v1, &row, k);
        }
      }
      else if (bc.isShared()) {
        for (unsigned int k=0; k<nblock; k++)
          sg_f->BlockReplaceGlobalValues(1, &zero, &row, k);
      }
    }
    
    // Always zero out row (This takes care of the not-owned case)
    // Do this for each SG block
    if (bc.isOwned() || bc.isShared()) {
      for (unsigned int i=0; i<nblock; i++)
        for (unsigned int j=0; j<nblock; j++) {
          sg_jac->BlockExtractGlobalRowView(row, num_entries, row_view, i, j);
          for (int k=0; k<num_entries; k++)
            row_view[k] = 0.0;
        }
    }
    
    // Check derivative array is nonzero
    if (node_f[offsets[eq_row]].hasFastAccess()) {
	    
      // Loop over equations per node
      for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
	      
        // Global column
        col = static_cast<int>(firstDOF + eq_col);
        
        // Replace Jacobian
        if (bc.isOwned()) {
          for (unsigned int k=0; k<nblock; k++) {
            v1 = node_f[eq_row].fastAccessDx(eq_col).coeff(k);
            nl = Cijk->num_values(k);
            for (unsigned int l=0; l<nl; l++) {
              Cijk->triple_value(k,l,i,j,cijk); 
              v2 = v1 * cijk / Cijk->norm_squared(i);
              sg_jac->BlockSumIntoGlobalValues(row, 1, &v2, &col, i, j);
            }
          }
        } // SG blocks
        
      } // column equations
      
    } // has fast access
    
  } // row equations

}

void
FEApp::SGJacobianOp::finalizeFill()
{
}

FEApp::SGMatrixFreeJacobianOp::SGMatrixFreeJacobianOp(
       const Teuchos::RCP<const Epetra_Map>& base_map,
       const Teuchos::RCP<const Epetra_CrsGraph>& base_graph,
       const Teuchos::RCP<const Stokhos::OrthogPolyBasis<double> >& sg_basis_,
       double alpha, double beta,
       const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
       const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x,
       const Teuchos::RCP<EpetraExt::BlockVector>& sg_overlapped_f) :
  nblock(sg_basis_->size()),
  map(base_map),
  graph(base_graph),
  sg_basis(sg_basis_),
  m_coeff(alpha),
  j_coeff(beta),
  sg_xdot(sg_overlapped_xdot),
  sg_x(sg_overlapped_x),
  sg_f(sg_overlapped_f),
  xdot(sg_xdot != Teuchos::null ? nblock : 0),
  x(nblock),
  jac(nblock)
{
  for (unsigned int block=0; block<nblock; block++) {
    // Allocate vectors for blocks
    if (sg_xdot != Teuchos::null)
      xdot[block] = Teuchos::rcp(new Epetra_Vector(*map));
    x[block] = Teuchos::rcp(new Epetra_Vector(*map));
    jac[block] = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));

    // Extract blocks
    sg_x->ExtractBlockValues(*x[block], block);
    if (sg_xdot != Teuchos::null)
      sg_xdot->ExtractBlockValues(*xdot[block], block);
  }
}

FEApp::SGMatrixFreeJacobianOp::~SGMatrixFreeJacobianOp()
{
}

void
FEApp::SGMatrixFreeJacobianOp::reset(
	  double alpha, double beta,
	  const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_xdot,
	  const Teuchos::RCP<const EpetraExt::BlockVector>& sg_overlapped_x)
{
  m_coeff = alpha;
  j_coeff = beta;
  sg_xdot = sg_overlapped_xdot;
  sg_x = sg_overlapped_x;

  // Extract blocks
  for (unsigned int block=0; block<nblock; block++) {
    sg_x->ExtractBlockValues(*x[block], block);
    if (sg_xdot != Teuchos::null)
      sg_xdot->ExtractBlockValues(*xdot[block], block);
    jac[block]->PutScalar(0.0);
  }
}

std::vector< Teuchos::RCP<Epetra_CrsMatrix> >&
FEApp::SGMatrixFreeJacobianOp::getJacobianBlocks()
{
  return jac;
}


void
FEApp::SGMatrixFreeJacobianOp::elementInit(const FEApp::AbstractElement& e,
                                           unsigned int neqn,
                                           std::vector< SGFadType >* elem_xdot,
                                           std::vector< SGFadType >& elem_x)
{
  // Global node ID
  unsigned int node_GID;

  // Local ID of first DOF
  unsigned int firstDOF;

  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Number of dof
  unsigned int ndof = nnode*neqn;

  // Copy element solution
  for (unsigned int i=0; i<nnode; i++) {
    node_GID = e.nodeGID(i);
    firstDOF = map->LID(node_GID*neqn);
    for (unsigned int j=0; j<neqn; j++) {
      elem_x[neqn*i+j] = SGFadType(ndof, 0.0);
      elem_x[neqn*i+j].fastAccessDx(neqn*i+j) = j_coeff;
      elem_x[neqn*i+j].val().resize(nblock);
      elem_x[neqn*i+j].val().copyForWrite();
      for (unsigned int block=0; block<nblock; block++)
        elem_x[neqn*i+j].val().fastAccessCoeff(block) = 
          (*x[block])[firstDOF+j];
      if (elem_xdot != NULL) {
        (*elem_xdot)[neqn*i+j] = SGFadType(ndof, 0.0);
        (*elem_xdot)[neqn*i+j].fastAccessDx(neqn*i+j) = m_coeff;
        (*elem_xdot)[neqn*i+j].val().resize(nblock);
        (*elem_xdot)[neqn*i+j].val().copyForWrite();
        for (unsigned int block=0; block<nblock; block++)
          (*elem_xdot)[neqn*i+j].val().fastAccessCoeff(block) = 
            (*xdot[block])[firstDOF+j];
      }
      
    }
  }

}

void
FEApp::SGMatrixFreeJacobianOp::elementPost(const FEApp::AbstractElement& e,
				 unsigned int neqn,
				 std::vector< SGFadType >& elem_f)
{
  // Number of nodes
  unsigned int nnode = e.numNodes();

  // Sum element residual and Jacobian into global residual, Jacobian
  // Loop over nodes in element
  int row, col;
  unsigned int lrow, lcol;
  double c;
  
  for (unsigned int node_row=0; node_row<nnode; node_row++) {
    
    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<neqn; eq_row++) {
      lrow = neqn*node_row+eq_row;
      
      // Global row
      row = static_cast<int>(e.nodeGID(node_row)*neqn + eq_row);

      // Sum residual
      if (sg_f != Teuchos::null) 
        for (unsigned int block=0; block<nblock; block++) {
          c = elem_f[lrow].val().coeff(block);
          sg_f->BlockSumIntoGlobalValues(1, &c, &row, block);
        }
      
      // Check derivative array is nonzero
      if (elem_f[lrow].hasFastAccess()) {
        
        // Loop over nodes in element
        for (unsigned int node_col=0; node_col<nnode; node_col++){
          
          // Loop over equations per node
          for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
            lcol = neqn*node_col+eq_col;
            
            // Global column
            col = static_cast<int>(e.nodeGID(node_col)*neqn + eq_col);
            
            // Sum Jacobian
            for (unsigned int block=0; block<nblock; block++) {
              c = elem_f[lrow].fastAccessDx(lcol).coeff(block);
              jac[block]->SumIntoGlobalValues(row, 1, &c, &col);
            }
            
          } // column equations
          
        } // column nodes
        
      } // has fast access
      
    } // row equations
      
  } // row node

}

void
FEApp::SGMatrixFreeJacobianOp::nodeInit(const FEApp::NodeBC& bc,
                                        unsigned int neqn,
                                        std::vector< SGFadType >* node_xdot,
                                        std::vector< SGFadType >& node_x)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Local ID of first DOF
  unsigned int firstDOF = map->LID(node_GID*neqn);

  // Copy element solution
  for (unsigned int j=0; j<neqn; j++) {
    node_x[j] = SGFadType(neqn, 0.0);
    node_x[j].fastAccessDx(j) = j_coeff;
    node_x[j].val().resize(nblock);
    node_x[j].val().copyForWrite();
    for (unsigned int block=0; block<nblock; block++)
      node_x[j].val().fastAccessCoeff(block) = (*x[block])[firstDOF+j];
    if (node_xdot != NULL) {
      (*node_xdot)[j] = SGFadType(neqn, 0.0);
      (*node_xdot)[j].fastAccessDx(j) = m_coeff;
      (*node_xdot)[j].val().resize(nblock);
      (*node_xdot)[j].val().copyForWrite();
      for (unsigned int block=0; block<nblock; block++)
        (*node_xdot)[j].val().fastAccessCoeff(block) = 
          (*xdot[block])[firstDOF+j];
    }
  }
}

void
FEApp::SGMatrixFreeJacobianOp::nodePost(const FEApp::NodeBC& bc,
                                        unsigned int neqn,
                                        std::vector< SGFadType >& node_f)
{
  // Global node ID
  unsigned int node_GID = bc.getNodeGID();

  // Global ID of first DOF
  unsigned int firstDOF = node_GID*neqn;

  // Residual offsets
  const std::vector<unsigned int>& offsets = bc.getOffsets();

  int num_entries;
  double* row_view = 0;
  int row, col;
  double c;
  double zero = 0.0;

  for (unsigned int block=0; block<nblock; block++) {

    // Loop over equations per node
    for (unsigned int eq_row=0; eq_row<offsets.size(); eq_row++) {

      // Global row
      row = static_cast<int>(firstDOF + offsets[eq_row]);
    
      // Replace residual
      if (sg_f != Teuchos::null) {
        if (bc.isOwned()) {
          c = node_f[offsets[eq_row]].val().coeff(block);
          sg_f->BlockReplaceGlobalValues(1, &c, &row, block);
        }
        else if (bc.isShared())
          sg_f->BlockReplaceGlobalValues(1, &zero, &row, block);
      }
      
      // Always zero out row (This takes care of the not-owned case)
      if (bc.isOwned() || bc.isShared()) {
        jac[block]->ExtractGlobalRowView(row, num_entries, row_view);
        for (int k=0; k<num_entries; k++)
          row_view[k] = 0.0;
      }
      
      // Check derivative array is nonzero
      if (node_f[offsets[eq_row]].hasFastAccess()) {
        
        // Loop over equations per node
        for (unsigned int eq_col=0; eq_col<neqn; eq_col++) {
          
          // Global column
          col = static_cast<int>(firstDOF + eq_col);
          
          // Replace Jacobian
          if (bc.isOwned()) {
            c = node_f[eq_row].fastAccessDx(eq_col).coeff(block);
            jac[block]->ReplaceGlobalValues(row, 1, &c, &col);
          }
          
        } // column equations
        
      } // has fast access
      
    } // row equations

  } // SG blocks
  
}

void
FEApp::SGMatrixFreeJacobianOp::finalizeFill()
{
}

#endif // SG_ACTIVE
