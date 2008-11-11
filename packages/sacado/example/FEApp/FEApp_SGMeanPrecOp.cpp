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

#include "FEApp_SGMeanPrecOp.hpp"

#if SG_ACTIVE

FEApp::SGMeanPrecOp::SGMeanPrecOp(
   const Teuchos::RCP<const Epetra_Map>& base_map_,
   const Teuchos::RCP<const Epetra_Map>& sg_map_,
   unsigned int num_blocks_,
   const Teuchos::RCP<Epetra_CrsMatrix>& mean_jac_,
   const Teuchos::RCP<Teuchos::ParameterList>& precParams_) :
  label("FEApp::SGMeanPrecOp"),
  base_map(base_map_),
  sg_map(sg_map_),
  mean_jac(mean_jac_),
  precParams(precParams_),
  useTranspose(false),
  num_blocks(num_blocks_),
  sg_input(),
  sg_result(),
  input_block(num_blocks),
  result_block(num_blocks)
{
  reset();
}

FEApp::SGMeanPrecOp::~SGMeanPrecOp()
{
}

int
FEApp::SGMeanPrecOp::reset()
{
  // Compute Ifpack preconditioner
  Ifpack Factory;
  std::string prec = precParams->get("Ifpack Preconditioner", "ILU");
  int overlap = precParams->get("Overlap", 0);
  ifpackPrec = Teuchos::rcp(Factory.Create(prec, mean_jac.get(), overlap));
  ifpackPrec->SetParameters(*precParams);
  int err = ifpackPrec->Initialize();   
  err = ifpackPrec->Compute();

  return err;
}

Teuchos::RCP<Epetra_CrsMatrix>
FEApp::SGMeanPrecOp::getMeanJacobian()
{
  return mean_jac;
}

int 
FEApp::SGMeanPrecOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  mean_jac->SetUseTranspose(useTranspose);

  return 0;
}

int 
FEApp::SGMeanPrecOp::Apply(const Epetra_MultiVector& Input, 
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
  }

  // Fill input blocks
  sg_input->Scale(1.0, Input);
  for (unsigned int i=0; i<num_blocks; i++) {
    sg_input->ExtractBlockValues(*input_block[i], i);
    result_block[i]->PutScalar(0.0);
  }

  // Apply mean Jacobian block
  for (unsigned int i=0; i<num_blocks; i++) {
    ifpackPrec->Apply(*input_block[i], *result_block[i]);
  }

  // Get result from blocks
  for (unsigned int i=0; i<num_blocks; i++)
    sg_result->LoadBlockValues(*result_block[i], i);
  Result.Scale(1.0, *sg_result);

  return 0;
}

int 
FEApp::SGMeanPrecOp::ApplyInverse(const Epetra_MultiVector& Input, 
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
  }

  // Fill input blocks
  sg_input->Scale(1.0, Input);
  for (unsigned int i=0; i<num_blocks; i++) {
    sg_input->ExtractBlockValues(*input_block[i], i);
    result_block[i]->PutScalar(0.0);
  }

  // Apply mean Jacobian block
  for (unsigned int i=0; i<num_blocks; i++) {
    ifpackPrec->ApplyInverse(*input_block[i], *result_block[i]);
  }

  // Get result from blocks
  for (unsigned int i=0; i<num_blocks; i++)
    sg_result->LoadBlockValues(*result_block[i], i);
  Result.Scale(1.0, *sg_result);

  return 0;
}

double 
FEApp::SGMeanPrecOp::NormInf() const
{
  return mean_jac->NormInf();
}


const char* 
FEApp::SGMeanPrecOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
FEApp::SGMeanPrecOp::UseTranspose() const
{
  return useTranspose;
}

bool 
FEApp::SGMeanPrecOp::HasNormInf() const
{
  return true;
}

const Epetra_Comm & 
FEApp::SGMeanPrecOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
FEApp::SGMeanPrecOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
FEApp::SGMeanPrecOp::OperatorRangeMap() const
{
  return *sg_map;
}

#endif
