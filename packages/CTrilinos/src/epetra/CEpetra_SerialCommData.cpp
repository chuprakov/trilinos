#include "CEpetra_SerialCommData_Cpp.hpp"
#include "CEpetra_SerialCommData.h"
#include "Epetra_SerialCommData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_SerialCommData>& tableOfSerialCommDatas()
{
    static Table<Epetra_SerialCommData> loc_tableOfSerialCommDatas(CT_Epetra_SerialCommData_ID, "CT_Epetra_SerialCommData_ID");
    return loc_tableOfSerialCommDatas;
}


} // namespace


//
// Definitions from CEpetra_SerialCommData.h
//


extern "C" {


CT_Epetra_SerialCommData_ID_t Epetra_SerialCommData_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfSerialCommDatas(), id);
}


} // extern "C"


//
// Definitions from CEpetra_SerialCommData_Cpp.hpp
//


const Teuchos::RCP<Epetra_SerialCommData>
CEpetra::getSerialCommData( CT_Epetra_SerialCommData_ID_t id )
{
    return tableOfSerialCommDatas().get(id);
}

CT_Epetra_SerialCommData_ID_t
CEpetra::storeSerialCommData( const Epetra_SerialCommData *pobj )
{
    return tableOfSerialCommDatas().store(pobj);
}

void
CEpetra::purgeSerialCommDataTable(  )
{
    tableOfSerialCommDatas().purge();
}



