#include "CEpetra_BlockMapData_Cpp.hpp"
#include "CEpetra_BlockMapData.h"
#include "Epetra_BlockMapData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_BlockMapData>& tableOfBlockMapDatas()
{
    static Table<Epetra_BlockMapData> loc_tableOfBlockMapDatas(CT_Epetra_BlockMapData_ID, "CT_Epetra_BlockMapData_ID");
    return loc_tableOfBlockMapDatas;
}


} // namespace


//
// Definitions from CEpetra_BlockMapData.h
//


extern "C" {


CT_Epetra_BlockMapData_ID_t Epetra_BlockMapData_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfBlockMapDatas(), id);
}


} // extern "C"


//
// Definitions from CEpetra_BlockMapData_Cpp.hpp
//


const Teuchos::RCP<Epetra_BlockMapData>
CEpetra::getBlockMapData( CT_Epetra_BlockMapData_ID_t id )
{
    return tableOfBlockMapDatas().get(id);
}

CT_Epetra_BlockMapData_ID_t
CEpetra::storeBlockMapData( const Epetra_BlockMapData *pobj )
{
    return tableOfBlockMapDatas().store(pobj);
}

void
CEpetra::purgeBlockMapDataTable(  )
{
    tableOfBlockMapDatas().purge();
}



