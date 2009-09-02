#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_EXPORT_H
#define CEPETRA_EXPORT_H


#include "CEpetra_BlockMap.h"
#include "CEpetra_Distributor.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_Export_Cast() and Epetra_Export_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_Export table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_Export_ID_t Epetra_Export_Cast ( CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_Export_Abstract ( 
  CT_Epetra_Export_ID_t id );

/* Original C++ prototype:
   Epetra_Export( const Epetra_BlockMap & SourceMap, const Epetra_BlockMap & TargetMap );
*/
CT_Epetra_Export_ID_t Epetra_Export_Create ( 
  CT_Epetra_BlockMap_ID_t SourceMapID, 
  CT_Epetra_BlockMap_ID_t TargetMapID );

/* Original C++ prototype:
   Epetra_Export(const Epetra_Export& Exporter);
*/
CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( 
  CT_Epetra_Export_ID_t ExporterID );

/* Original C++ prototype:
   virtual ~Epetra_Export(void);
*/
void Epetra_Export_Destroy ( CT_Epetra_Export_ID_t * selfID );

/* Original C++ prototype:
   int NumSameIDs() const;
*/
int Epetra_Export_NumSameIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int NumPermuteIDs() const;
*/
int Epetra_Export_NumPermuteIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int * PermuteFromLIDs () const;
*/
int * Epetra_Export_PermuteFromLIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int * PermuteToLIDs () const;
*/
int * Epetra_Export_PermuteToLIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int NumRemoteIDs() const;
*/
int Epetra_Export_NumRemoteIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int * RemoteLIDs() const;
*/
int * Epetra_Export_RemoteLIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int NumExportIDs () const;
*/
int Epetra_Export_NumExportIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int * ExportLIDs () const;
*/
int * Epetra_Export_ExportLIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int * ExportPIDs () const;
*/
int * Epetra_Export_ExportPIDs ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int NumSend() const;
*/
int Epetra_Export_NumSend ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   int NumRecv() const;
*/
int Epetra_Export_NumRecv ( CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   const Epetra_BlockMap & SourceMap() const;
*/
CT_Epetra_BlockMap_ID_t Epetra_Export_SourceMap ( 
  CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   const Epetra_BlockMap & TargetMap() const;
*/
CT_Epetra_BlockMap_ID_t Epetra_Export_TargetMap ( 
  CT_Epetra_Export_ID_t selfID );

/* Original C++ prototype:
   Epetra_Distributor & Distributor() const;
*/
CT_Epetra_Distributor_ID_t Epetra_Export_Distributor ( 
  CT_Epetra_Export_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_EXPORT_H */

