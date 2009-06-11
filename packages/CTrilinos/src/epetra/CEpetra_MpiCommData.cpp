#ifdef HAVE_MPI


#include "CEpetra_MpiCommData_Cpp.hpp"
#include "CEpetra_MpiCommData.h"
#include "Epetra_MpiCommData.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_MpiCommData>& tableOfMpiCommDatas()
{
    static Table<Epetra_MpiCommData> loc_tableOfMpiCommDatas(CT_Epetra_MpiCommData_ID, "CT_Epetra_MpiCommData_ID");
    return loc_tableOfMpiCommDatas;
}


} // namespace


//
// Definitions from CEpetra_MpiCommData.h
//


extern "C" {


CT_Epetra_MpiCommData_ID_t Epetra_MpiCommData_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfMpiCommDatas(), id);
}


} // extern "C"


//
// Definitions from CEpetra_MpiCommData_Cpp.hpp
//


const Teuchos::RCP<Epetra_MpiCommData>
CEpetra::getMpiCommData( CT_Epetra_MpiCommData_ID_t id )
{
    return tableOfMpiCommDatas().get(id);
}

CT_Epetra_MpiCommData_ID_t
CEpetra::storeMpiCommData( const Epetra_MpiCommData *pobj )
{
    return tableOfMpiCommDatas().store(pobj);
}

void
CEpetra::purgeMpiCommDataTable(  )
{
    tableOfMpiCommDatas().purge();
}



#endif /* HAVE_MPI */


