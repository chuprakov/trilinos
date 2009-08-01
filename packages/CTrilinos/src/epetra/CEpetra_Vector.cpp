#include "CTrilinos_config.h"

#include "Epetra_DataAccess.h"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Vector.h"
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_Vector */
Table<Epetra_Vector>& tableOfVectors()
{
    static Table<Epetra_Vector>
        loc_tableOfVectors(CT_Epetra_Vector_ID, "CT_Epetra_Vector_ID", false);
    return loc_tableOfVectors;
}

/* table to hold objects of type const Epetra_Vector */
Table<const Epetra_Vector>& tableOfConstVectors()
{
    static Table<const Epetra_Vector>
        loc_tableOfConstVectors(CT_Epetra_Vector_ID, "CT_Epetra_Vector_ID", true);
    return loc_tableOfConstVectors;
}


} // namespace


//
// Definitions from CEpetra_Vector.h
//


extern "C" {


CT_Epetra_Vector_ID_t Epetra_Vector_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstVectors(), id);
    } else {
        newid = CTrilinos::cast(tableOfVectors(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_Vector_Abstract ( 
  CT_Epetra_Vector_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id);
}

CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        tableOfVectors().store(new Epetra_Vector(
        *CEpetra::getBlockMap(MapID), zeroOut)));
}

CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        tableOfVectors().store(new Epetra_Vector(
        *CEpetra::getVector(SourceID))));
}

CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, double * V )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        tableOfVectors().store(new Epetra_Vector(
        CV, *CEpetra::getBlockMap(MapID), V)));
}

CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
        tableOfVectors().store(new Epetra_Vector(
        CV, *CEpetra::getMultiVector(SourceID), Index)));
}

void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstVectors().remove(&aid);
    } else {
        tableOfVectors().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(aid);
}

int Epetra_Vector_ReplaceGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceGlobalValues(
        NumEntries, Values, Indices);
}

int Epetra_Vector_ReplaceMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceMyValues(
        NumEntries, Values, Indices);
}

int Epetra_Vector_SumIntoGlobalValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoGlobalValues(
        NumEntries, Values, Indices);
}

int Epetra_Vector_SumIntoMyValues ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, double * Values, 
  int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoMyValues(
        NumEntries, Values, Indices);
}

int Epetra_Vector_ReplaceGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceGlobalValues(
        NumEntries, BlockOffset, Values, Indices);
}

int Epetra_Vector_ReplaceMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->ReplaceMyValues(
        NumEntries, BlockOffset, Values, Indices);
}

int Epetra_Vector_SumIntoGlobalValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoGlobalValues(
        NumEntries, BlockOffset, Values, Indices);
}

int Epetra_Vector_SumIntoMyValues_BlockPos ( 
  CT_Epetra_Vector_ID_t selfID, int NumEntries, int BlockOffset, 
  double * Values, int * Indices )
{
    return CEpetra::getVector(selfID)->SumIntoMyValues(
        NumEntries, BlockOffset, Values, Indices);
}

int Epetra_Vector_ExtractCopy ( 
  CT_Epetra_Vector_ID_t selfID, double * V )
{
    return CEpetra::getConstVector(selfID)->ExtractCopy(V);
}

int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V )
{
    return CEpetra::getConstVector(selfID)->ExtractView(V);
}

double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index )
{
    const Epetra_Vector& self = *( CEpetra::getConstVector(selfID) );

    return self[index];
}


} // extern "C"


//
// Definitions from CEpetra_Vector_Cpp.hpp
//


/* get Epetra_Vector from non-const table using CT_Epetra_Vector_ID */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CT_Epetra_Vector_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id);
    return tableOfVectors().get(aid);
}

/* get Epetra_Vector from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CTrilinos_Object_ID_t id )
{
    return tableOfVectors().get(id);
}

/* get const Epetra_Vector from either the const or non-const table
 * using CT_Epetra_Vector_ID */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CT_Epetra_Vector_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_Vector_ID_t>(id);
    if (id.is_const) {
        return tableOfConstVectors().get(aid);
    } else {
        return tableOfVectors().get(aid);
    }
}

/* get const Epetra_Vector from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_Vector>
CEpetra::getConstVector( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstVectors().get(id);
    } else {
        return tableOfVectors().get(id);
    }
}

/* store Epetra_Vector in non-const table */
CT_Epetra_Vector_ID_t
CEpetra::storeVector( Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
            tableOfVectors().storeCopy(pobj));
}

/* store const Epetra_Vector in const table */
CT_Epetra_Vector_ID_t
CEpetra::storeConstVector( const Epetra_Vector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_Vector_ID_t>(
            tableOfConstVectors().storeCopy(pobj));
}

/* dump contents of Epetra_Vector and const Epetra_Vector tables */
void
CEpetra::purgeVectorTables(  )
{
    tableOfVectors().purge();
    tableOfConstVectors().purge();
}



