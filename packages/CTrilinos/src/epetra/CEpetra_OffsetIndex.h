#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_OFFSETINDEX_H
#define CEPETRA_OFFSETINDEX_H


#include "CEpetra_CrsGraph.h"
#include "CEpetra_Import.h"
#include "CEpetra_Export.h"
#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


/* Functions Epetra_OffsetIndex_Cast() and Epetra_OffsetIndex_Abstract()
   are used for casting CTrilinos objects from one type to another.
   The former function performs a dynamic cast on the underlying object
   and stores an RCP to it in the Epetra_OffsetIndex table, while
   the latter only converts the type of the struct that references the
   object so that an object of any type can be passed to the former
   function (use the _Abstract() function corresponding to the type
   of the object that will be casted, not the type to which it will
   be casted).
*/

CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Cast ( 
  CTrilinos_Object_ID_t id );

CTrilinos_Object_ID_t Epetra_OffsetIndex_Abstract ( 
  CT_Epetra_OffsetIndex_ID_t id );

/* Original C++ prototype:
   Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph, const Epetra_CrsGraph & TargetGraph, Epetra_Import & Importer );
*/
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromImporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Import_ID_t ImporterID );

/* Original C++ prototype:
   Epetra_OffsetIndex( const Epetra_CrsGraph & SourceGraph, const Epetra_CrsGraph & TargetGraph, Epetra_Export & Exporter );
*/
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Create_FromExporter ( 
  CT_Epetra_CrsGraph_ID_t SourceGraphID, 
  CT_Epetra_CrsGraph_ID_t TargetGraphID, 
  CT_Epetra_Export_ID_t ExporterID );

/* Original C++ prototype:
   Epetra_OffsetIndex(const Epetra_OffsetIndex & Indexor);
*/
CT_Epetra_OffsetIndex_ID_t Epetra_OffsetIndex_Duplicate ( 
  CT_Epetra_OffsetIndex_ID_t IndexorID );

/* Original C++ prototype:
   virtual ~Epetra_OffsetIndex(void);
*/
void Epetra_OffsetIndex_Destroy ( 
  CT_Epetra_OffsetIndex_ID_t * selfID );

/* Original C++ prototype:
   int ** SameOffsets() const;
*/
int ** Epetra_OffsetIndex_SameOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );

/* Original C++ prototype:
   int ** PermuteOffsets() const;
*/
int ** Epetra_OffsetIndex_PermuteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );

/* Original C++ prototype:
   int ** RemoteOffsets() const;
*/
int ** Epetra_OffsetIndex_RemoteOffsets ( 
  CT_Epetra_OffsetIndex_ID_t selfID );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_OFFSETINDEX_H */

