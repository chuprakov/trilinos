// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_ApproxGaussSeidelPreconditioner.hpp"
#include "Epetra_config.h"
#include "Teuchos_TimeMonitor.hpp"

Stokhos::ApproxGaussSeidelPreconditioner::
ApproxGaussSeidelPreconditioner(
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  const Teuchos::RCP<Stokhos::PreconditionerFactory>& prec_factory_,
  const Teuchos::RCP<Teuchos::ParameterList>& params_) :
  label("Stokhos Approximate Gauss-Seidel Preconditioner"),
  base_map(base_map_),
  sg_map(sg_map_),
  prec_factory(prec_factory_),
  mean_prec(),
  useTranspose(false),
  sg_op(),
  sg_poly(),
  Cijk(),
  symmetric(false),
  only_use_linear(true),
  mat_vec_tmp(),
  rhs_block()
{
  symmetric = params_->get("Symmetric Gauss-Seidel", false);
  only_use_linear = params_->get("Only Use Linear Terms", true);
}

Stokhos::ApproxGaussSeidelPreconditioner::
~ApproxGaussSeidelPreconditioner()
{
}

void
Stokhos::ApproxGaussSeidelPreconditioner::
setupPreconditioner(const Teuchos::RCP<Stokhos::SGOperator>& sg_op_, 
		    const Epetra_Vector& x)
{
  sg_op = sg_op_;
  sg_poly = sg_op->getSGPolynomial();
  mean_prec = prec_factory->compute(sg_poly->getCoeffPtr(0));
  label = std::string("Stokhos Approximate Gauss-Seidel Preconditioner:\n") + 
    std::string("		***** ") + 
    std::string(mean_prec->Label());
  Cijk = sg_op->getTripleProduct();
}

int 
Stokhos::ApproxGaussSeidelPreconditioner::
SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;

  return 0;
}

int 
Stokhos::ApproxGaussSeidelPreconditioner::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  return sg_op->Apply(Input, Result);
}

int 
Stokhos::ApproxGaussSeidelPreconditioner::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Total Approximate Gauss-Seidel Time");

  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  } 

  int m = input->NumVectors();
  if (mat_vec_tmp == Teuchos::null || mat_vec_tmp->NumVectors() != m)
    mat_vec_tmp = Teuchos::rcp(new Epetra_MultiVector(*base_map, m));
  if (rhs_block == Teuchos::null || rhs_block->NumVectors() != m)
    rhs_block = 
      Teuchos::rcp(new EpetraExt::BlockMultiVector(*base_map, *sg_map, m));
  
  // Extract blocks
  EpetraExt::BlockMultiVector input_block(View, *base_map, *input);
  EpetraExt::BlockMultiVector result_block(View, *base_map, Result);

  result_block.PutScalar(0.0);

  int sz = sg_poly->basis()->size();
  int k_limit = sg_poly->size();
  if (only_use_linear)
    k_limit = sg_poly->basis()->dimension() + 1;
  const Teuchos::Array<double>& norms = sg_poly->basis()->norm_squared();

  rhs_block->Update(1.0, input_block, 0.0);

  for (int i=0; i<sz; i++) {

    Teuchos::RCP<Epetra_MultiVector> res_i = result_block.GetBlock(i);
    {
      // Apply deterministic preconditioner
      TEUCHOS_FUNC_TIME_MONITOR("Total AGS Deterministic Preconditioner Time");
      mean_prec->ApplyInverse(*(rhs_block->GetBlock(i)), *res_i);
    }

    for (Cijk_type::ik_iterator k_it = Cijk->k_begin(i);
	 k_it != Cijk->k_end(i); ++k_it) {
      int k = index(k_it);
      if (k!=0 && k<k_limit) {
	bool do_mat_vec = false;
	for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	     j_it != Cijk->j_end(k_it); ++j_it) {
	  int j = index(j_it);
	  if (j > i) {
	    do_mat_vec = true;
	    break;
	  }
	}
	if (do_mat_vec) {
	  (*sg_poly)[k].Apply(*res_i, *mat_vec_tmp);
	  for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    int j = index(j_it);
	    double c = value(j_it);
	    if (j > i) {
	      rhs_block->GetBlock(j)->Update(-1.0*c/norms[j], *mat_vec_tmp, 1.0);
	    }
	  }
	}
      }
    }
  }

  // For symmetric Gauss-Seidel
  if (symmetric) {
    
    rhs_block->Update(1.0, input_block, 0.0);

    for (int i=sz-1; i>=0; --i) {
       
      Teuchos::RCP<Epetra_MultiVector> res_i = result_block.GetBlock(i);
      {
	// Apply deterministic preconditioner
	TEUCHOS_FUNC_TIME_MONITOR("Total AGS Deterministic Preconditioner Time");
	mean_prec->ApplyInverse(*(rhs_block->GetBlock(i)), *res_i);
      }

      for (Cijk_type::ik_iterator k_it = Cijk->k_begin(i);
	   k_it != Cijk->k_end(i); ++k_it) {
	int k = index(k_it);
	if (k!=0 && k<k_limit) {
	  bool do_mat_vec = false;
	  for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
	       j_it != Cijk->j_end(k_it); ++j_it) {
	    int j = index(j_it);
	    if (j > i) {
	      do_mat_vec = true;
	      break;
	    }
	  }
	  if (do_mat_vec) {
	    (*sg_poly)[k].Apply(*res_i, *mat_vec_tmp);
	    for (Cijk_type::ikj_iterator j_it = Cijk->j_begin(k_it);
		 j_it != Cijk->j_end(k_it); ++j_it) {
	      int j = index(j_it);
	      double c = value(j_it);
	      if (j > i) {
		rhs_block->GetBlock(j)->Update(-1.0*c/norms[j], *mat_vec_tmp, 1.0);
	      }
	    }
	  }
	}
      }
    }
  }

  if (made_copy)
    delete input;

  return 0; 
}

double 
Stokhos::ApproxGaussSeidelPreconditioner::
NormInf() const
{
  return sg_op->NormInf();
}


const char* 
Stokhos::ApproxGaussSeidelPreconditioner::
Label() const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::ApproxGaussSeidelPreconditioner::
UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::ApproxGaussSeidelPreconditioner::
HasNormInf() const
{
  return sg_op->HasNormInf();
}

const Epetra_Comm & 
Stokhos::ApproxGaussSeidelPreconditioner::
Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::ApproxGaussSeidelPreconditioner::
OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::ApproxGaussSeidelPreconditioner::
OperatorRangeMap() const
{
  return *sg_map;
}
