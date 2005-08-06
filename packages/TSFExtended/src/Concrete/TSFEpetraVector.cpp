/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/

#include "TSFEpetraVector.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "TSFVectorImpl.hpp"
#include "RTOp_parallel_helpers.h"
#include "RTOpPack_MPI_apply_op.hpp"
#include "Teuchos_Workspace.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_dyn_cast.hpp"
#include "Thyra_MPIVectorStd.hpp"



using namespace Teuchos;
using namespace TSFExtended;
using namespace Thyra;
using Thyra::Index;


EpetraVector::EpetraVector(const RefCountPtr<const VectorSpaceBase<double> >& vs)
  : MPIVectorStd<double>(), epetraVec_(), mpiVecSpace_(), epetraMap_()
{
  const EpetraVectorSpace* epvs 
    = dynamic_cast<const EpetraVectorSpace*>(vs.get());
  TEST_FOR_EXCEPTION(epvs==0, runtime_error,
                     "could not cast vector space to EpetraVectorSpace in "
                     "EpetraVector ctor");

  mpiVecSpace_ = rcp_dynamic_cast<const MPIVectorSpaceBase<double> >(vs);

  epetraMap_ = epvs->epetraMap();
  epetraVec_ = rcp(new Epetra_Vector(*epetraMap_, true));

  RefCountPtr<double> data = rcp(&(epetraVec_->operator[](0)), false);
  initialize(mpiVecSpace_, data, 1);
}


double& EpetraVector::operator[](Index globalIndex) 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(globalIndex)];
}

void EpetraVector::setElement(Index index, const double& value)
{
  epetraVec()->ReplaceGlobalValues(1, const_cast<double*>(&value), 
                                    const_cast<int*>(&index));
}

void EpetraVector::addToElement(Index index, const double& value)
{
  epetraVec()->SumIntoGlobalValues(1, const_cast<double*>(&value), 
                                    const_cast<int*>(&index));
}

const double& EpetraVector::getElement(Index index) const 
{
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  return (*epetraVec())[myMap.LID(index)];
}

void EpetraVector::getElements(const Index* globalIndices, int numElems,
                               vector<double>& elems) const
{
  elems.resize(numElems);
  const Epetra_BlockMap& myMap = epetraVec()->Map();
  RefCountPtr<const Epetra_Vector> epv = epetraVec();

  for (int i=0; i<numElems; i++)
    {
      elems[i] = (*epv)[myMap.LID(globalIndices[i])];
    }
}

void EpetraVector::setElements(size_t numElems, const Index* globalIndices,
                               const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->ReplaceGlobalValues(numElems, globalIndices, values);
}

void EpetraVector::addToElements(size_t numElems, const Index* globalIndices,
                                 const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  int ierr = vec->SumIntoGlobalValues(numElems, globalIndices, values);
}

void EpetraVector::finalizeAssembly()
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetraVec().get());
  vec->GlobalAssemble();
}

const Epetra_Vector& EpetraVector::getConcrete(const TSFExtended::Vector<double>& tsfVec)
{
  const EpetraVector* epv 
    = dynamic_cast<const EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
                     "EpetraVector::getConcrete called on a vector that "
                     "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}

Epetra_Vector& EpetraVector::getConcrete(TSFExtended::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
                     "EpetraVector::getConcrete called on a vector that "
                     "could not be cast to an EpetraVector");
  return *(epv->epetraVec());
}


Epetra_Vector* EpetraVector::getConcretePtr(TSFExtended::Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
                     "EpetraVector::getConcrete called on a vector that "
                     "could not be cast to an EpetraVector");
  return epv->epetraVec().get();
}


