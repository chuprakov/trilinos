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

#include "Epetra_config.h"
#include "Epetra_MultiVector.h"

#include "FEApp_SGMatrixFreeOp.hpp"

#if SG_ACTIVE

FEApp::SGMatrixFreeOp::SGMatrixFreeOp(
   const Teuchos::RCP<const Epetra_Map>& base_map_,
   const Teuchos::RCP<const Epetra_Map>& sg_map_,
   const Teuchos::RCP<const tp_type >& Cijk_,
   const Teuchos::RCP<std::vector< Teuchos::RCP<Epetra_CrsMatrix> > >& jacs_) :
  label("FEApp::SGMatrixFreeOp"),
  base_map(base_map_),
  sg_map(sg_map_),
  Cijk(Cijk_),
  jacs(jacs_),
  useTranspose(false),
  num_blocks(jacs->size()),
  sg_input(),
  sg_result(),
  input_block(num_blocks),
  result_block(num_blocks),
  tmp()
{
}

FEApp::SGMatrixFreeOp::~SGMatrixFreeOp()
{
}

std::vector< Teuchos::RCP<Epetra_CrsMatrix> >&
FEApp::SGMatrixFreeOp::getJacobianBlocks()
{
  return *jacs;
}

int 
FEApp::SGMatrixFreeOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  for (unsigned int i=0; i<num_blocks; i++)
    (*jacs)[i]->SetUseTranspose(useTranspose);

  return 0;
}

int 
FEApp::SGMatrixFreeOp::Apply(const Epetra_MultiVector& Input, 
                             Epetra_MultiVector& Result) const
{
  int m = Input.NumVectors();
  if (sg_input == Teuchos::null || sg_input->NumVectors() != m) {
    sg_input = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
    sg_result = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
    for (unsigned int i=0; i<num_blocks; i++) {
      input_block[i] = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
      result_block[i] = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
    }
    tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
  }

  // Fill input blocks
  sg_input->Scale(1.0, Input);
  for (unsigned int i=0; i<num_blocks; i++) {
    sg_input->ExtractBlockValues(*input_block[i], i);
    result_block[i]->PutScalar(0.0);
  }

  // Apply block SG Jacobian via
  // w_i = 
  //    \sum_{j=0}^P \sum_{k=0}^P J_k v_j < \psi_i \psi_j \psi_k > / <\psi_i^2>
  // for i=0,...,P where P = num_blocks w_j is the jth input block, w_i
  // is the ith result block, and J_k is the kth block Jacobian
  double cijk;
  unsigned int i, j;
  for (unsigned int k=0; k<num_blocks; k++) {
    unsigned int nl = Cijk->num_values(k);
    for (unsigned int l=0; l<nl; l++) {
      Cijk->triple_value(k, l, i, j, cijk);
      cijk /= Cijk->norm_squared(i);
      (*jacs)[k]->Apply(*input_block[j], *tmp);
      result_block[i]->Update(cijk, *tmp, 1.0);
    }
  }

  // Get result from blocks
  for (unsigned int i=0; i<num_blocks; i++)
    sg_result->LoadBlockValues(*result_block[i], i);
  Result.Scale(1.0, *sg_result);

  return 0;
}

int 
FEApp::SGMatrixFreeOp::ApplyInverse(const Epetra_MultiVector& Input, 
				    Epetra_MultiVector& Result) const
{
  throw "SGMatrixFreeOp::ApplyInverse not defined!";
  return -1;
}

double 
FEApp::SGMatrixFreeOp::NormInf() const
{
  return 1.0;
}


const char* 
FEApp::SGMatrixFreeOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
FEApp::SGMatrixFreeOp::UseTranspose() const
{
  return useTranspose;
}

bool 
FEApp::SGMatrixFreeOp::HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
FEApp::SGMatrixFreeOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
FEApp::SGMatrixFreeOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
FEApp::SGMatrixFreeOp::OperatorRangeMap() const
{
  return *sg_map;
}

#endif
