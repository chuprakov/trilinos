#include "TSFEpetraVector.hpp"
#include "TSFVector.hpp"

using namespace Teuchos;
using namespace TSFExtended;
using TSFCore::Index;


EpetraVector::EpetraVector(const RefCountPtr<Epetra_Vector>& vec,
                           const RefCountPtr<const TSFCore::EpetraVectorSpace>& vs)
   : TSFCore::EpetraVector(vec, vs) 
{;}

void EpetraVector::setElement(Index index, const double& value)
{
  (*epetra_vec())[index] = value;
}

void EpetraVector::addToElement(Index index, const double& value)
{
  (*epetra_vec())[index] += value;
}

const double& EpetraVector::getElement(Index index) const 
{
  return (*epetra_vec())[index];
}

void EpetraVector::setElements(size_t numElems, const Index* globalIndices,
                               const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetra_vec().get());
  int ierr = vec->ReplaceGlobalValues(numElems, globalIndices, values);
}

void EpetraVector::addToElements(size_t numElems, const Index* globalIndices,
                                 const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetra_vec().get());
  int ierr = vec->SumIntoGlobalValues(numElems, globalIndices, values);
}

void EpetraVector::finalizeAssembly()
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetra_vec().get());
  vec->GlobalAssemble();
}

const Epetra_Vector& EpetraVector::getConcrete(const Vector<double>& tsfVec)
{
  const EpetraVector* epv 
    = dynamic_cast<const EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
                     "EpetraVector::getConcrete called on a vector that "
                     "could not be cast to an EpetraVector");
  return *(epv->epetra_vec());
}

Epetra_Vector& EpetraVector::getConcrete(Vector<double>& tsfVec)
{
  EpetraVector* epv 
    = dynamic_cast<EpetraVector*>(tsfVec.ptr().get());
  TEST_FOR_EXCEPTION(epv==0, std::runtime_error,
                     "EpetraVector::getConcrete called on a vector that "
                     "could not be cast to an EpetraVector");
  return *(epv->epetra_vec());
}

