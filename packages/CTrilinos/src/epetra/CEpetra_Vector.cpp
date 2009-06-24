#include "CEpetra_BlockMap_Cpp.hpp"
#include "Epetra_DataAccess.h"
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


Table<Epetra_Vector>& tableOfVectors()
{
    static Table<Epetra_Vector> loc_tableOfVectors(CT_Epetra_Vector_ID, "CT_Epetra_Vector_ID");
    return loc_tableOfVectors;
}


} // namespace


//
// Definitions from CEpetra_Vector.h
//


extern "C" {


CT_Epetra_Vector_ID_t Epetra_Vector_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfVectors(), id);
}

CT_Epetra_Vector_ID_t Epetra_Vector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return tableOfVectors().store(new Epetra_Vector(*pMap, zeroOut));

}

CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( 
  CT_Epetra_Vector_ID_t SourceID )
{
    const Teuchos::RCP<Epetra_Vector> 
        pSource = CEpetra::getVector(SourceID);

    return tableOfVectors().store(new Epetra_Vector(*pSource));

}

CT_Epetra_Vector_ID_t Epetra_Vector_Create_FromArray ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, double * V )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return tableOfVectors().store(new Epetra_Vector(CV, *pMap, V));

}

CT_Epetra_Vector_ID_t Epetra_Vector_FromSource ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int Index )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pSource = CEpetra::getMultiVector(SourceID);

    return tableOfVectors().store(new Epetra_Vector(
        CV, *pSource, Index));

}

void Epetra_Vector_Destroy ( CT_Epetra_Vector_ID_t * selfID )
{
    tableOfVectors().remove(selfID);
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
    return CEpetra::getVector(selfID)->ExtractCopy(V);
}

int Epetra_Vector_ExtractView ( 
  CT_Epetra_Vector_ID_t selfID, double ** V )
{
    return CEpetra::getVector(selfID)->ExtractView(V);
}

double Epetra_Vector_getElement ( 
  CT_Epetra_Vector_ID_t selfID, int index )
{
    Epetra_Vector& self = *( CEpetra::getVector(selfID) );

    return self[index];
}


} // extern "C"


//
// Definitions from CEpetra_Vector_Cpp.hpp
//


const Teuchos::RCP<Epetra_Vector>
CEpetra::getVector( CT_Epetra_Vector_ID_t id )
{
    return tableOfVectors().get(id);
}

CT_Epetra_Vector_ID_t
CEpetra::storeVector( const Epetra_Vector *pobj )
{
    return tableOfVectors().store(pobj);
}

void
CEpetra::purgeVectorTable(  )
{
    tableOfVectors().purge();
}



