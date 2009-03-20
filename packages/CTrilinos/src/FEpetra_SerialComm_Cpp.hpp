#include "FEpetra_SerialComm.h"
#include "Teuchos_RCP.hpp"

#ifndef FEPETRA_SERIALCOMM_CPP_HPP
#define FEPETRA_SERIALCOMM_CPP_HPP

class Epetra_SerialComm;

namespace FEpetra{
  using Teuchos::RCP;

  const RCP<const Epetra_SerialComm> getSerialComm( SerialCommID id);
} // namespace FEpetra

#endif // FEPETRA_SERIALCOMM_CPP_HPP
~                                         
