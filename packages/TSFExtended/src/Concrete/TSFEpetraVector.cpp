#include "TSFEpetraVector.hpp"

using namespace Teuchos;
using namespace TSF;
using namespace TSF::Internal;

EpetraVector::EpetraVector(const RefCountPtr<Epetra_FEVector>& vec)
  : TSFCore::EpetraVector(vec), 
    PrintableVector()
{;}

void EpetraVector::setElement(int index, const double& value)
{
  (*epetra_vec())[index] = value;
}

double EpetraVector::getElement(int index) const 
{
  return (*epetra_vec())[index];
}

void EpetraVector::setElements(int numElems, const int* globalIndices,
                               const double* values)
{
  Epetra_FEVector* vec = dynamic_cast<Epetra_FEVector*>(epetra_vec().get());
  int ierr = vec->ReplaceGlobalValues(numElems, globalIndices, values);
}

void EpetraVector::addToElements(int numElems, const int* globalIndices,
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
