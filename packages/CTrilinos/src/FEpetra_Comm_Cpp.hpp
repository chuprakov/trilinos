#ifndef FEPETRA_COMM_CPP_HPP
#define FEPETRA_COMM_CPP_HPP

#include "FEpetra_Comm.h"
#include "Teuchos_RCP.hpp"

class Epetra_Comm;

namespace FEpetra{
  // using Teuchos::RCP;

  const Teuchos::RCP<const Epetra_Comm> getComm( CommID id );
} // namespace FEpetra

#endif // FEPETRA_COMM_CPP_HPP
~                                         
