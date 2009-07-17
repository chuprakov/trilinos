#include "CTrilinos_config.h"

#include "Epetra_CombineMode.h"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_DistObject.h"
#include "Epetra_DistObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_DistObject>& tableOfDistObjects()
{
    static Table<Epetra_DistObject>
        loc_tableOfDistObjects(CT_Epetra_DistObject_ID, "CT_Epetra_DistObject_ID");
    return loc_tableOfDistObjects;
}


} // namespace


//
// Definitions from CEpetra_DistObject.h
//


extern "C" {


CT_Epetra_DistObject_ID_t Epetra_DistObject_Cast ( 
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        CTrilinos::cast(tableOfDistObjects(), id));
}

CTrilinos_Object_ID_t Epetra_DistObject_Abstract ( 
  CT_Epetra_DistObject_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id);
}

void Epetra_DistObject_Destroy ( 
  CT_Epetra_DistObject_ID_t * selfID )
{
    CTrilinos_Object_ID_t id =
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(*selfID);
    tableOfDistObjects().remove(&id);
    *selfID = CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(id);
}

int Epetra_DistObject_Import ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<Epetra_SrcDistObject> 
        pA = CEpetra::getSrcDistObject(AID);
    const Teuchos::RCP<Epetra_Import> 
        pImporter = CEpetra::getImport(ImporterID);
    const Teuchos::RCP<Epetra_OffsetIndex> 
        pIndexor = CEpetra::getOffsetIndex(IndexorID);

    return CEpetra::getDistObject(selfID)->Import(
        *pA, *pImporter, CombineMode, &(*pIndexor));
}

int Epetra_DistObject_Import_UsingExporter ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<Epetra_SrcDistObject> 
        pA = CEpetra::getSrcDistObject(AID);
    const Teuchos::RCP<Epetra_Export> 
        pExporter = CEpetra::getExport(ExporterID);
    const Teuchos::RCP<Epetra_OffsetIndex> 
        pIndexor = CEpetra::getOffsetIndex(IndexorID);

    return CEpetra::getDistObject(selfID)->Import(
        *pA, *pExporter, CombineMode, &(*pIndexor));
}

int Epetra_DistObject_Export_UsingImporter ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Import_ID_t ImporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<Epetra_SrcDistObject> 
        pA = CEpetra::getSrcDistObject(AID);
    const Teuchos::RCP<Epetra_Import> 
        pImporter = CEpetra::getImport(ImporterID);
    const Teuchos::RCP<Epetra_OffsetIndex> 
        pIndexor = CEpetra::getOffsetIndex(IndexorID);

    return CEpetra::getDistObject(selfID)->Export(
        *pA, *pImporter, CombineMode, &(*pIndexor));
}

int Epetra_DistObject_Export ( 
  CT_Epetra_DistObject_ID_t selfID, CT_Epetra_SrcDistObject_ID_t AID, 
  CT_Epetra_Export_ID_t ExporterID, Epetra_CombineMode CombineMode, 
  CT_Epetra_OffsetIndex_ID_t IndexorID )
{
    const Teuchos::RCP<Epetra_SrcDistObject> 
        pA = CEpetra::getSrcDistObject(AID);
    const Teuchos::RCP<Epetra_Export> 
        pExporter = CEpetra::getExport(ExporterID);
    const Teuchos::RCP<Epetra_OffsetIndex> 
        pIndexor = CEpetra::getOffsetIndex(IndexorID);

    return CEpetra::getDistObject(selfID)->Export(
        *pA, *pExporter, CombineMode, &(*pIndexor));
}

CT_Epetra_BlockMap_ID_t Epetra_DistObject_Map ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeBlockMap(
        &( CEpetra::getDistObject(selfID)->Map() ));
}

CT_Epetra_Comm_ID_t Epetra_DistObject_Comm ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::storeComm(
        &( CEpetra::getDistObject(selfID)->Comm() ));
}

boolean Epetra_DistObject_DistributedGlobal ( 
  CT_Epetra_DistObject_ID_t selfID )
{
    return CEpetra::getDistObject(selfID)->DistributedGlobal();
}


} // extern "C"


//
// Definitions from CEpetra_DistObject_Cpp.hpp
//


const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CT_Epetra_DistObject_ID_t id )
{
    return tableOfDistObjects().get(
        CTrilinos::abstractType<CT_Epetra_DistObject_ID_t>(id));
}

const Teuchos::RCP<Epetra_DistObject>
CEpetra::getDistObject( CTrilinos_Object_ID_t id )
{
    return tableOfDistObjects().get(id);
}

CT_Epetra_DistObject_ID_t
CEpetra::storeDistObject( const Epetra_DistObject *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_DistObject_ID_t>(
        tableOfDistObjects().storeCopy(pobj));
}

void
CEpetra::purgeDistObjectTable(  )
{
    tableOfDistObjects().purge();
}



