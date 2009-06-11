#ifdef HAVE_MPI



/* True C header file! */


#ifndef CEPETRA_MPICOMMDATA_H
#define CEPETRA_MPICOMMDATA_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_MpiCommData_ID_t Epetra_MpiCommData_Cast(
  CTrilinos_Object_ID_t id );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_MPICOMMDATA_H */

#endif /* HAVE_MPI */


