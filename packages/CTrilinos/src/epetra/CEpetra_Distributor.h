#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_DISTRIBUTOR_H
#define CEPETRA_DISTRIBUTOR_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_Distributor_ID_t Epetra_Distributor_Cast ( 
  CTrilinos_Object_ID_t id );

/* Original C++ prototype:
   virtual Epetra_Distributor * Clone() = 0;
*/
CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID );

/* Original C++ prototype:
   virtual ~Epetra_Distributor();
*/
void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID );

/* Original C++ prototype:
   virtual int CreateFromSends( const int & NumExportIDs, const int * ExportPIDs, bool Deterministic, int & NumRemoteIDs ) = 0;
*/
int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs );

/* Original C++ prototype:
   virtual int CreateFromRecvs( const int & NumRemoteIDs, const int * RemoteGIDs, const int * RemotePIDs, bool Deterministic, int & NumExportIDs, int *& ExportGIDs, int *& ExportPIDs) = 0;
*/
int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs );

/* Original C++ prototype:
   virtual int Do( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs) = 0;
*/
int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );

/* Original C++ prototype:
   virtual int DoReverse( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs ) = 0;
*/
int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );

/* Original C++ prototype:
   virtual int DoPosts( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs ) = 0;
*/
int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );

/* Original C++ prototype:
   virtual int DoWaits() = 0;
*/
int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID );

/* Original C++ prototype:
   virtual int DoReversePosts( char * export_objs, int obj_size, int & len_import_objs, char *& import_objs) = 0;
*/
int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs );

/* Original C++ prototype:
   virtual int DoReverseWaits() = 0;
*/
int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID );

/* Original C++ prototype:
   virtual int Do( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
*/
int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );

/* Original C++ prototype:
   virtual int DoReverse( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
*/
int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );

/* Original C++ prototype:
   virtual int DoPosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
*/
int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );

/* Original C++ prototype:
   virtual int DoReversePosts( char * export_objs, int obj_size, int *& sizes, int & len_import_objs, char *& import_objs) = 0;
*/
int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_DISTRIBUTOR_H */

