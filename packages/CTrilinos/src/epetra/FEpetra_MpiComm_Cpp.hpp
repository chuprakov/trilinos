#ifdef HAVE_MPI

#ifndef FEPETRA_MPICOMM_CPP_HPP
#define FEPETRA_MPICOMM_CPP_HPP

#include "FEpetra_MpiComm.h"
#include "Teuchos_RCP.hpp"

class Epetra_MpiComm;

namespace FEpetra{
  // using Teuchos::RCP;

  const Teuchos::RCP<const Epetra_MpiComm> getMpiComm( FEpetra_MpiCommID id );
} // namespace FEpetra

#endif // FEPETRA_MPICOMM_CPP_HPP        

#endif /* HAVE_MPI */
