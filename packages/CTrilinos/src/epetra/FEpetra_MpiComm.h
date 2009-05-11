#ifdef HAVE_MPI

/* True C header file! */
#ifndef FEPETRA_MPICOMM_H
#define FEPETRA_MPICOMM_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif
  typedef int FEpetra_MpiCommID;

  FEpetra_MpiCommID Epetra_MpiComm_Create(MPI_Comm comm);

  void Epetra_MpiComm_Destroy( FEpetra_MpiCommID id );
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FEPETRA_MPICOMM_H */

#endif /* HAVE_MPI */
