#include "CEpetra_Data_Cpp.hpp"
#include "CEpetra_Data.h"
#include "Epetra_Data.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Data>& tableOfDatas()
{
    static Table<Epetra_Data> loc_tableOfDatas(CT_Epetra_Data_ID, "CT_Epetra_Data_ID");
    return loc_tableOfDatas;
}


} // namespace


//
// Definitions from CEpetra_Data.h
//


extern "C" {


CT_Epetra_Data_ID_t Epetra_Data_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfDatas(), id);
}


} // extern "C"


//
// Definitions from CEpetra_Data_Cpp.hpp
//


const Teuchos::RCP<Epetra_Data>
CEpetra::getData( CT_Epetra_Data_ID_t id )
{
    return tableOfDatas().get(id);
}

CT_Epetra_Data_ID_t
CEpetra::storeData( const Epetra_Data *pobj )
{
    return tableOfDatas().store(pobj);
}

void
CEpetra::purgeDataTable(  )
{
    tableOfDatas().purge();
}



