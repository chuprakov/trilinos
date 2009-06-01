#include "Teuchos_Exceptions.hpp"


#ifndef CTRILINOS_EXCEPTIONS_HPP
#define CTRILINOS_EXCEPTIONS_HPP


namespace CTrilinos {


class CTrilinosTypeMismatchError : public Teuchos::ExceptionBase
{
  public:
    CTrilinosTypeMismatchError(const std::string& what_arg) : Teuchos::ExceptionBase(what_arg) {}
};


} // namespace CTrilinos


#endif

