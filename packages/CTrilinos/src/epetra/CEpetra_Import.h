
/* True C header file! */


#ifndef CEPETRA_IMPORT_H
#define CEPETRA_IMPORT_H


#include "CEpetra_BlockMap.h"
#include "CEpetra_Distributor.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Import_ID_t Epetra_Import_Cast(
  CTrilinos_Object_ID_t id );

CT_Epetra_Import_ID_t Epetra_Import_Create ( 
  CT_Epetra_BlockMap_ID_t TargetMapID, 
  CT_Epetra_BlockMap_ID_t SourceMapID );

CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( 
  CT_Epetra_Import_ID_t ImporterID );

void Epetra_Import_Destroy ( CT_Epetra_Import_ID_t * selfID );

int Epetra_Import_NumSameIDs ( CT_Epetra_Import_ID_t selfID );

int Epetra_Import_NumPermuteIDs ( CT_Epetra_Import_ID_t selfID );

int * Epetra_Import_PermuteFromLIDs ( CT_Epetra_Import_ID_t selfID );

int * Epetra_Import_PermuteToLIDs ( CT_Epetra_Import_ID_t selfID );

int Epetra_Import_NumRemoteIDs ( CT_Epetra_Import_ID_t selfID );

int * Epetra_Import_RemoteLIDs ( CT_Epetra_Import_ID_t selfID );

int Epetra_Import_NumExportIDs ( CT_Epetra_Import_ID_t selfID );

int * Epetra_Import_ExportLIDs ( CT_Epetra_Import_ID_t selfID );

int * Epetra_Import_ExportPIDs ( CT_Epetra_Import_ID_t selfID );

int Epetra_Import_NumSend ( CT_Epetra_Import_ID_t selfID );

int Epetra_Import_NumRecv ( CT_Epetra_Import_ID_t selfID );

CT_Epetra_BlockMap_ID_t Epetra_Import_SourceMap ( 
  CT_Epetra_Import_ID_t selfID );

CT_Epetra_BlockMap_ID_t Epetra_Import_TargetMap ( 
  CT_Epetra_Import_ID_t selfID );

CT_Epetra_Distributor_ID_t Epetra_Import_Distributor ( 
  CT_Epetra_Import_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_IMPORT_H */

