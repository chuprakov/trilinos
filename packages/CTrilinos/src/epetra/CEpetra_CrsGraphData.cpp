#include "CEpetra_CrsGraphData_Cpp.hpp"
#include "CEpetra_CrsGraphData.h"
#include "Epetra_CrsGraphData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_CrsGraphData>& tableOfCrsGraphDatas()
{
    static Table<Epetra_CrsGraphData> loc_tableOfCrsGraphDatas(CT_Epetra_CrsGraphData_ID, "CT_Epetra_CrsGraphData_ID");
    return loc_tableOfCrsGraphDatas;
}


} // namespace


//
// Definitions from CEpetra_CrsGraphData.h
//


extern "C" {


CT_Epetra_CrsGraphData_ID_t Epetra_CrsGraphData_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfCrsGraphDatas(), id);
}


} // extern "C"


//
// Definitions from CEpetra_CrsGraphData_Cpp.hpp
//


const Teuchos::RCP<Epetra_CrsGraphData>
CEpetra::getCrsGraphData( CT_Epetra_CrsGraphData_ID_t id )
{
    return tableOfCrsGraphDatas().get(id);
}

CT_Epetra_CrsGraphData_ID_t
CEpetra::storeCrsGraphData( const Epetra_CrsGraphData *pobj )
{
    return tableOfCrsGraphDatas().store(pobj);
}

void
CEpetra::purgeCrsGraphDataTable(  )
{
    tableOfCrsGraphDatas().purge();
}



