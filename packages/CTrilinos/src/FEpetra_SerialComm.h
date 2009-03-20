/* True C header file! */
#ifndef FEPETRA_SERIALCOMM_H
#define FEPETRA_SERIALCOMM_H

#ifdef __cplusplus
extern "C" {
#endif
  typedef int SerialCommID;

  SerialCommID Epetra_SerialComm_Create();

  void Epetra_SerialComm_Destroy( SerialCommID id);
#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FEPETRA_SERIALCOMM_H */
