
/* True C header file! */


#ifndef CEPETRA_DISTRIBUTOR_H
#define CEPETRA_DISTRIBUTOR_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Distributor_ID_t Epetra_Distributor_Cast(
  CTrilinos_Object_ID_t id );

CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID );

void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID );

int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID,  const int * NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs );

int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID,  const int * NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs );

int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );

int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int * len_import_objs, char ** import_objs );

int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID );

int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );

int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );

int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );

int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID,  char * export_objs, 
  int obj_size, int ** sizes, int * len_import_objs, 
  char ** import_objs );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DISTRIBUTOR_H */

