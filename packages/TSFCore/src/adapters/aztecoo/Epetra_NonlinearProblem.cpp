// /////////////////////////////////////////////////////////////////
// Epetra_NonlinearProblem.cpp

#include "Epetra_NonlinearProblem.hpp"

namespace Epetra {

NonlinearProblem::Scalar
NonlinearProblem::infiniteBound()
{
	return 1e+50;
}
	
// Basic information

int NonlinearProblem::Nu() const
{
  assert(0);
  return 0;
}
	
int NonlinearProblem::numResponseFunctions() const
{
  assert(0);
  return 0;
}

// VectorSpaces

Teuchos::RefCountPtr<const Epetra_Map>
NonlinearProblem::map_u(int l) const
{
  assert(0);
  return Teuchos::null;
}

Teuchos::RefCountPtr<const Epetra_Map>
NonlinearProblem::map_g() const
{
  assert(0);
  return Teuchos::null;
}

// Bounds

const Epetra_Vector&
NonlinearProblem::uL(int l) const
{
  assert(0);
  return uL(l);
}

const Epetra_Vector&
NonlinearProblem::uU(int l) const
{
  assert(0);
  return uU(l);
}

const Epetra_Vector&
NonlinearProblem::gL() const
{
  assert(0);
  return gL();
}

const Epetra_Vector&
NonlinearProblem::gU() const
{
  assert(0);
  return gU();
}

// Initial values (guesses) for state and auxiliary variables

const Epetra_Vector&
NonlinearProblem::u0(int l) const
{
  assert(0);
  return u0(l);
}

// Calculation methods

void NonlinearProblem::calc_g(
  const Epetra_Vector     &y
  ,const Epetra_Vector*   u[]
  ,Epetra_Vector          *g
  ) const
{
  assert(0);
}

} // namespace Epetra
