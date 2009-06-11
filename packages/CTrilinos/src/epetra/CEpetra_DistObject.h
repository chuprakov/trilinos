
/* True C header file! */


#ifndef CEPETRA_DISTOBJECT_H
#define CEPETRA_DISTOBJECT_H


#include "CEpetra_SrcDistObject.h"
#include "CEpetra_Import.h"
#include "CEpetra_OffsetIndex.h"
#include "Epetra_CombineMode.h"
#include "CEpetra_Export.h"
#include "CEpetra_BlockMap.h"
#include "CEpetra_Comm.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_DistObject_ID_t Epetra_DistObject_Cast(
  CTrilinos_Object_ID_t id );

void Epetra_DistObject_Destroy ( 
  CT_Epetra_DistObject_ID_t * selfID );

int Epetra_DistObject_Import ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

int Epetra_DistObject_Import_UsingExporter ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

int Epetra_DistObject_Export_UsingImporter ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

int Epetra_DistObject_Export ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( 
  CT_Epetra_DistObject_ID_t selfID );

CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( 
  CT_Epetra_DistObject_ID_t selfID );

boolean Epetra_DistObject_DistributedGlobal ( 
  CT_Epetra_DistObject_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DISTOBJECT_H */

