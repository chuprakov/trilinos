#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_Distributor.h"
#include "Epetra_Distributor.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Distributor>& tableOfDistributors()
{
    static Table<Epetra_Distributor> loc_tableOfDistributors(CT_Epetra_Distributor_ID, "CT_Epetra_Distributor_ID");
    return loc_tableOfDistributors;
}


} // namespace


//
// Definitions from CEpetra_Distributor.h
//


extern "C" {


CT_Epetra_Distributor_ID_t Epetra_Distributor_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfDistributors(), id);
}

CT_Epetra_Distributor_ID_t Epetra_Distributor_Clone ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::storeDistributor(
        CEpetra::getDistributor(selfID)->Clone());
}

void Epetra_Distributor_Destroy ( 
  CT_Epetra_Distributor_ID_t * selfID )
{
    tableOfDistributors().remove(selfID);
}

int Epetra_Distributor_CreateFromSends ( 
  CT_Epetra_Distributor_ID_t selfID, int NumExportIDs, 
  const int * ExportPIDs, boolean Deterministic, 
  int * NumRemoteIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromSends(
        NumExportIDs, ExportPIDs, Deterministic, *NumRemoteIDs);
}

int Epetra_Distributor_CreateFromRecvs ( 
  CT_Epetra_Distributor_ID_t selfID, int NumRemoteIDs, 
  const int * RemoteGIDs, const int * RemotePIDs, 
  boolean Deterministic, int * NumExportIDs, int ** ExportGIDs, 
  int ** ExportPIDs )
{
    return CEpetra::getDistributor(selfID)->CreateFromRecvs(
        NumRemoteIDs, RemoteGIDs, RemotePIDs, Deterministic, *NumExportIDs, *ExportGIDs, *ExportPIDs);
}

int Epetra_Distributor_Do ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(
        export_objs, obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(
        export_objs, obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(
        export_objs, obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoWaits ( CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoWaits();
}

int Epetra_Distributor_DoReversePosts ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(
        export_objs, obj_size, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverseWaits ( 
  CT_Epetra_Distributor_ID_t selfID )
{
    return CEpetra::getDistributor(selfID)->DoReverseWaits();
}

int Epetra_Distributor_Do_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->Do(
        export_objs, obj_size, *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReverse_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReverse(
        export_objs, obj_size, *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoPosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoPosts(
        export_objs, obj_size, *sizes, *len_import_objs, *import_objs);
}

int Epetra_Distributor_DoReversePosts_VarLen ( 
  CT_Epetra_Distributor_ID_t selfID, char * export_objs, int obj_size, 
  int ** sizes, int * len_import_objs, char ** import_objs )
{
    return CEpetra::getDistributor(selfID)->DoReversePosts(
        export_objs, obj_size, *sizes, *len_import_objs, *import_objs);
}


} // extern "C"


//
// Definitions from CEpetra_Distributor_Cpp.hpp
//


const Teuchos::RCP<Epetra_Distributor>
CEpetra::getDistributor( CT_Epetra_Distributor_ID_t id )
{
    return tableOfDistributors().get(id);
}

CT_Epetra_Distributor_ID_t
CEpetra::storeDistributor( const Epetra_Distributor *pobj )
{
    return tableOfDistributors().storeCopy(pobj);
}

void
CEpetra::purgeDistributorTable(  )
{
    tableOfDistributors().purge();
}



