#include "CTrilinos_config.h"

#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_CompObject.h"
#include "Epetra_CompObject.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


Table<Epetra_CompObject>& tableOfCompObjects()
{
    static Table<Epetra_CompObject> loc_tableOfCompObjects(CT_Epetra_CompObject_ID, "CT_Epetra_CompObject_ID");
    return loc_tableOfCompObjects;
}


} // namespace


//
// Definitions from CEpetra_CompObject.h
//


extern "C" {


CT_Epetra_CompObject_ID_t Epetra_CompObject_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfCompObjects(), id);
}

CT_Epetra_CompObject_ID_t Epetra_CompObject_Create (  )
{
    return tableOfCompObjects().store(new Epetra_CompObject());

}

CT_Epetra_CompObject_ID_t Epetra_CompObject_Duplicate ( 
  CT_Epetra_CompObject_ID_t SourceID )
{
    const Teuchos::RCP<Epetra_CompObject> 
        pSource = CEpetra::getCompObject(SourceID);

    return tableOfCompObjects().store(new Epetra_CompObject(*pSource));

}

void Epetra_CompObject_Destroy ( 
  CT_Epetra_CompObject_ID_t * selfID )
{
    tableOfCompObjects().remove(selfID);
}

void Epetra_CompObject_SetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_Flops_ID_t FlopCounter_inID )
{
    const Teuchos::RCP<Epetra_Flops> 
        pFlopCounter_in = CEpetra::getFlops(FlopCounter_inID);

    CEpetra::getCompObject(selfID)->SetFlopCounter(*pFlopCounter_in);
}

void Epetra_CompObject_SetFlopCounter_Matching ( 
  CT_Epetra_CompObject_ID_t selfID, 
  CT_Epetra_CompObject_ID_t CompObjectID )
{
    const Teuchos::RCP<Epetra_CompObject> 
        pCompObject = CEpetra::getCompObject(CompObjectID);

    CEpetra::getCompObject(selfID)->SetFlopCounter(*pCompObject);
}

void Epetra_CompObject_UnsetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    CEpetra::getCompObject(selfID)->UnsetFlopCounter();
}

CT_Epetra_Flops_ID_t Epetra_CompObject_GetFlopCounter ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    return CEpetra::storeFlops(
        CEpetra::getCompObject(selfID)->GetFlopCounter());
}

void Epetra_CompObject_ResetFlops ( 
  CT_Epetra_CompObject_ID_t selfID )
{
    CEpetra::getCompObject(selfID)->ResetFlops();
}

double Epetra_CompObject_Flops ( CT_Epetra_CompObject_ID_t selfID )
{
    return CEpetra::getCompObject(selfID)->Flops();
}

void Epetra_CompObject_UpdateFlops_Int ( 
  CT_Epetra_CompObject_ID_t selfID, int Flops_in )
{
    CEpetra::getCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Long ( 
  CT_Epetra_CompObject_ID_t selfID, long int Flops_in )
{
    CEpetra::getCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Double ( 
  CT_Epetra_CompObject_ID_t selfID, double Flops_in )
{
    CEpetra::getCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_UpdateFlops_Float ( 
  CT_Epetra_CompObject_ID_t selfID, float Flops_in )
{
    CEpetra::getCompObject(selfID)->UpdateFlops(Flops_in);
}

void Epetra_CompObject_Assign ( 
  CT_Epetra_CompObject_ID_t selfID, CT_Epetra_CompObject_ID_t srcID )
{
    Epetra_CompObject& self = *( CEpetra::getCompObject(selfID) );

    const Teuchos::RCP<Epetra_CompObject> 
        psrc = CEpetra::getCompObject(srcID);

    self = *psrc;
}


} // extern "C"


//
// Definitions from CEpetra_CompObject_Cpp.hpp
//


const Teuchos::RCP<Epetra_CompObject>
CEpetra::getCompObject( CT_Epetra_CompObject_ID_t id )
{
    return tableOfCompObjects().get(id);
}

CT_Epetra_CompObject_ID_t
CEpetra::storeCompObject( const Epetra_CompObject *pobj )
{
    return tableOfCompObjects().storeCopy(pobj);
}

void
CEpetra::purgeCompObjectTable(  )
{
    tableOfCompObjects().purge();
}



