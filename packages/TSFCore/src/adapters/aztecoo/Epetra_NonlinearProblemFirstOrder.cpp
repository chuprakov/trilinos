// /////////////////////////////////////////////////////////////////
// Epetra_NonlinearProblemFirstOrder.cpp

#include "Epetra_NonlinearProblemFirstOrder.hpp"

#define EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED( FNC_NAME ) \
  TEST_FOR_EXCEPTION( \
    true, std::logic_error \
    ,"Error, the function " << (FNC_NAME) << " given a default implementation in " \
    "Epetra::NonlinearProblemFirstOrder is not overridden in the class " \
    << typeid(*this).name() << " or should not be called by the client!" \
  )

namespace Epetra {

// Adjoints supported?

bool NonlinearProblemFirstOrder::adjointSupported() const
{
  return true;
}

// Factories for linear operators

bool NonlinearProblemFirstOrder::use_EO_DcDu(int l) const
{
  EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED("use_EO_DcDu(l)");
  return false;
}

Teuchos::RefCountPtr<Epetra_Operator>
NonlinearProblemFirstOrder::create_DcDu_op(int l) const
{
  EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED("create_DcDu_op(l)");
  return Teuchos::null;
}

bool NonlinearProblemFirstOrder::specialized_DcDy_prec() const
{
  return false;
}

Teuchos::RefCountPtr<Epetra_Operator>
NonlinearProblemFirstOrder::create_DcDy_prec() const
{
  return Teuchos::null;
}

Teuchos::RefCountPtr<Epetra_MultiVector>
NonlinearProblemFirstOrder::create_DcDu_mv(int l) const
{
  if(this->use_EO_DcDu(l)) {
    return Teuchos::null;
  }
  else {
    return Teuchos::rcp(
      new Epetra_MultiVector(
        *this->map_c()
        ,map_u(l)->NumGlobalElements()
        )
      );
  }
}

// Transpose arguments

ETransp NonlinearProblemFirstOrder::opDcDu(int l) const
{
  TEST_FOR_EXCEPTION(
    this->use_EO_DcDu(l), std::logic_error
    ,"NonlinearProblemFirstOrder::opDcDu(l): Error, the default implementation for "
    "this function is for the creation of Epetra_MultiVector objects for DcDu.  For "
    "different behavior, this function must be overridden by the subclass "
    << typeid(*this).name() << " or one of its base classes!" \
    );
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
  EPETRA_NONLINEARPROBLEMFIRSTORDER_NOT_DEFINED("calc_Dg(...)");
}

} // namespace Epetra
