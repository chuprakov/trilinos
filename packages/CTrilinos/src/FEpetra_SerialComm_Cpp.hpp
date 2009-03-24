#ifndef FEPETRA_SERIALCOMM_CPP_HPP
#define FEPETRA_SERIALCOMM_CPP_HPP

#include "FEpetra_SerialComm.h"
#include "Teuchos_RCP.hpp"

class Epetra_SerialComm;

namespace FEpetra{
  using Teuchos::RCP;

  const RCP<const Epetra_SerialComm> getSerialComm( FEpetra_SerialCommID id);
} // namespace FEpetra

#endif // FEPETRA_SERIALCOMM_CPP_HPP
