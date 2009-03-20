/* True C header file! */
#ifndef FEPETRA_COMM_H
#define FEPETRA_COMM_H
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif
  typedef int CommID;

  int Epetra_Comm_SumAll(double *PartialSums, double *GlobalSums,int Count,CommID comm);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FEPETRA_COMM_H */
