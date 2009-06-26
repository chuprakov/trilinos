#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_Flops.h"
#include "Epetra_Flops.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_Flops>& tableOfFlopss()
{
    static Table<Epetra_Flops> loc_tableOfFlopss(CT_Epetra_Flops_ID, "CT_Epetra_Flops_ID");
    return loc_tableOfFlopss;
}


} // namespace


//
// Definitions from CEpetra_Flops.h
//


extern "C" {


CT_Epetra_Flops_ID_t Epetra_Flops_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfFlopss(), id);
}

CT_Epetra_Flops_ID_t Epetra_Flops_Create (  )
{
    return tableOfFlopss().store(new Epetra_Flops());

}

CT_Epetra_Flops_ID_t Epetra_Flops_Duplicate ( 
  CT_Epetra_Flops_ID_t Flops_inID )
{
    const Teuchos::RCP<Epetra_Flops> 
        pFlops_in = CEpetra::getFlops(Flops_inID);

    return tableOfFlopss().store(new Epetra_Flops(*pFlops_in));

}

double Epetra_Flops_Flops ( CT_Epetra_Flops_ID_t selfID )
{
    return CEpetra::getFlops(selfID)->Flops();
}

void Epetra_Flops_ResetFlops ( CT_Epetra_Flops_ID_t selfID )
{
    CEpetra::getFlops(selfID)->ResetFlops();
}

void Epetra_Flops_Destroy ( CT_Epetra_Flops_ID_t * selfID )
{
    tableOfFlopss().remove(selfID);
}

void Epetra_Flops_Assign ( 
  CT_Epetra_Flops_ID_t selfID, CT_Epetra_Flops_ID_t srcID )
{
    Epetra_Flops& self = *( CEpetra::getFlops(selfID) );

    const Teuchos::RCP<Epetra_Flops> 
        psrc = CEpetra::getFlops(srcID);

    self = *psrc;
}


} // extern "C"


//
// Definitions from CEpetra_Flops_Cpp.hpp
//


const Teuchos::RCP<Epetra_Flops>
CEpetra::getFlops( CT_Epetra_Flops_ID_t id )
{
    return tableOfFlopss().get(id);
}

CT_Epetra_Flops_ID_t
CEpetra::storeFlops( const Epetra_Flops *pobj )
{
    return tableOfFlopss().storeCopy(pobj);
}

void
CEpetra::purgeFlopsTable(  )
{
    tableOfFlopss().purge();
}



