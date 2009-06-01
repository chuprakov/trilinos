#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Comm.h"
#include "Epetra_Comm.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Comm>& tableOfComms()
{
    static Table<Epetra_Comm> loc_tableOfComms(CT_Epetra_Comm_ID, "CT_Epetra_Comm_ID");
    return loc_tableOfComms;
}


} // namespace


//
// Definitions from CEpetra_Comm.h
//


extern "C" {


CT_Epetra_Comm_ID_t Epetra_Comm_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfComms(), id);
}

void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID )
{
    tableOfComms().remove(selfID);
}


} // extern "C"


//
// Definitions from CEpetra_Comm_Cpp.hpp
//


const Teuchos::RCP<Epetra_Comm>
CEpetra::getComm( CT_Epetra_Comm_ID_t id )
{
    return tableOfComms().get(id);
}

void
CEpetra::purgeCommTable(  )
{
    tableOfComms().purge();
}



