// /////////////////////////////////////////////////////////////////
// Epetra_NonlinearProblemFirstOrder.cpp

#include "Epetra_NonlinearProblemFirstOrder.hpp"

namespace Epetra {

// Adjoints supported?

bool NonlinearProblemFirstOrder::adjointSupported() const
{
  return true;
}

// Factories for linear operators

Teuchos::RefCountPtr<Epetra_Operator>
NonlinearProblemFirstOrder::create_DcDu(int l) const
{
  assert(0);
  return Teuchos::null;
}

// Transpose arguments

ETransp NonlinearProblemFirstOrder::opDcDu(int l) const
{
  assert(0);
  return NOTRANS;
}

// Calculation methods

void NonlinearProblemFirstOrder::calc_Dg(
  const Epetra_Vector           &y
  ,const Epetra_Vector*         u[]
  ,Epetra_Vector                *g
  ,Epetra_MultiVector           *DgDy
  ,Epetra_MultiVector*          DgDu[]
  ) const
{
  assert(0);
}

} // namespace Epetra
