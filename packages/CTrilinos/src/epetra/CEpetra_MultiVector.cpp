#include "CTrilinos_config.h"

#include "Epetra_DataAccess.h"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_MultiVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_MultiVector */
Table<Epetra_MultiVector>& tableOfMultiVectors()
{
    static Table<Epetra_MultiVector>
        loc_tableOfMultiVectors(CT_Epetra_MultiVector_ID, "CT_Epetra_MultiVector_ID", false);
    return loc_tableOfMultiVectors;
}

/* table to hold objects of type const Epetra_MultiVector */
Table<const Epetra_MultiVector>& tableOfConstMultiVectors()
{
    static Table<const Epetra_MultiVector>
        loc_tableOfConstMultiVectors(CT_Epetra_MultiVector_ID, "CT_Epetra_MultiVector_ID", true);
    return loc_tableOfConstMultiVectors;
}


} // namespace


//
// Definitions from CEpetra_MultiVector.h
//


extern "C" {


CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Cast ( 
  CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstMultiVectors(), id);
    } else {
        newid = CTrilinos::cast(tableOfMultiVectors(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_MultiVector_Abstract ( 
  CT_Epetra_MultiVector_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id);
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store(new Epetra_MultiVector(
        *CEpetra::getBlockMap(MapID), NumVectors, zeroOut)));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( 
  CT_Epetra_MultiVector_ID_t SourceID )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store(new Epetra_MultiVector(
        *CEpetra::getMultiVector(SourceID))));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, double * A, 
  int MyLDA, int NumVectors )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *CEpetra::getBlockMap(MapID), A, MyLDA, NumVectors)));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, 
  double ** ArrayOfPointers, int NumVectors )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *CEpetra::getBlockMap(MapID), ArrayOfPointers, NumVectors)));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int * Indices, int NumVectors )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *CEpetra::getMultiVector(SourceID), Indices, NumVectors)));
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int StartIndex, int NumVectors )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
        tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *CEpetra::getMultiVector(SourceID), StartIndex, NumVectors)));
}

void Epetra_MultiVector_Destroy ( 
  CT_Epetra_MultiVector_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstMultiVectors().remove(&aid);
    } else {
        tableOfMultiVectors().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(aid);
}

int Epetra_MultiVector_ReplaceGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceGlobalValue(
        GlobalRow, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceGlobalValue(
        GlobalBlockRow, BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoGlobalValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoGlobalValue(
        GlobalRow, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoGlobalValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int GlobalBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoGlobalValue(
        GlobalBlockRow, BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceMyValue(
        MyRow, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_ReplaceMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->ReplaceMyValue(
        MyBlockRow, BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoMyValue ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyRow, int VectorIndex, 
  double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoMyValue(
        MyRow, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_SumIntoMyValue_BlockPos ( 
  CT_Epetra_MultiVector_ID_t selfID, int MyBlockRow, 
  int BlockRowOffset, int VectorIndex, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->SumIntoMyValue(
        MyBlockRow, BlockRowOffset, VectorIndex, ScalarValue);
}

int Epetra_MultiVector_PutScalar ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarConstant )
{
    return CEpetra::getMultiVector(selfID)->PutScalar(ScalarConstant);
}

int Epetra_MultiVector_Random ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->Random();
}

int Epetra_MultiVector_ExtractCopy_Fill2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double * A, int MyLDA )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractCopy(A, MyLDA);
}

int Epetra_MultiVector_ExtractCopy_FillAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** ArrayOfPointers )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractCopy(
        ArrayOfPointers);
}

int Epetra_MultiVector_ExtractView_Set2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** A, int * MyLDA )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractView(A, MyLDA);
}

int Epetra_MultiVector_ExtractView_SetAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double *** ArrayOfPointers )
{
    return CEpetra::getConstMultiVector(selfID)->ExtractView(
        ArrayOfPointers);
}

int Epetra_MultiVector_Dot ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID, 
  double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->Dot(
        *CEpetra::getMultiVector(AID), Result);
}

int Epetra_MultiVector_Abs ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    return CEpetra::getMultiVector(selfID)->Abs(
        *CEpetra::getMultiVector(AID));
}

int Epetra_MultiVector_Reciprocal ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    return CEpetra::getMultiVector(selfID)->Reciprocal(
        *CEpetra::getMultiVector(AID));
}

int Epetra_MultiVector_Scale_Self ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarValue )
{
    return CEpetra::getMultiVector(selfID)->Scale(ScalarValue);
}

int Epetra_MultiVector_Scale ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID )
{
    return CEpetra::getMultiVector(selfID)->Scale(
        ScalarA, *CEpetra::getMultiVector(AID));
}

int Epetra_MultiVector_Update_WithA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarThis )
{
    return CEpetra::getMultiVector(selfID)->Update(
        ScalarA, *CEpetra::getMultiVector(AID), ScalarThis);
}

int Epetra_MultiVector_Update_WithAB ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarB, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    return CEpetra::getMultiVector(selfID)->Update(
        ScalarA, *CEpetra::getMultiVector(AID), ScalarB, *CEpetra::getMultiVector(BID), ScalarThis);
}

int Epetra_MultiVector_Norm1 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->Norm1(Result);
}

int Epetra_MultiVector_Norm2 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->Norm2(Result);
}

int Epetra_MultiVector_NormInf ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->NormInf(Result);
}

int Epetra_MultiVector_NormWeighted ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t WeightsID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->NormWeighted(
        *CEpetra::getMultiVector(WeightsID), Result);
}

int Epetra_MultiVector_MinValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->MinValue(Result);
}

int Epetra_MultiVector_MaxValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->MaxValue(Result);
}

int Epetra_MultiVector_MeanValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getConstMultiVector(selfID)->MeanValue(Result);
}

int Epetra_MultiVector_Multiply_Matrix ( 
  CT_Epetra_MultiVector_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    return CEpetra::getMultiVector(selfID)->Multiply(
        TransA, TransB, ScalarAB, *CEpetra::getMultiVector(AID), *CEpetra::getMultiVector(BID), ScalarThis);
}

int Epetra_MultiVector_Multiply_ByEl ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    return CEpetra::getMultiVector(selfID)->Multiply(
        ScalarAB, *CEpetra::getMultiVector(AID), *CEpetra::getMultiVector(BID), ScalarThis);
}

int Epetra_MultiVector_ReciprocalMultiply ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    return CEpetra::getMultiVector(selfID)->ReciprocalMultiply(
        ScalarAB, *CEpetra::getMultiVector(AID), *CEpetra::getMultiVector(BID), ScalarThis);
}

int Epetra_MultiVector_SetSeed ( 
  CT_Epetra_MultiVector_ID_t selfID, unsigned int Seed_in )
{
    return CEpetra::getMultiVector(selfID)->SetSeed(Seed_in);
}

unsigned int Epetra_MultiVector_Seed ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->Seed();
}

void Epetra_MultiVector_Assign ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t SourceID )
{
    Epetra_MultiVector& self = *( CEpetra::getMultiVector(selfID) );

    self = *CEpetra::getMultiVector(SourceID);
}

double * Epetra_MultiVector_getArray ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    const Epetra_MultiVector& self = *( CEpetra::getConstMultiVector(selfID) );

    return self[i];
}

CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    const Epetra_MultiVector& self = *( CEpetra::getConstMultiVector(selfID) );

    return CEpetra::storeConstVector(self(i));
}

int Epetra_MultiVector_NumVectors ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->NumVectors();
}

int Epetra_MultiVector_MyLength ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->MyLength();
}

int Epetra_MultiVector_GlobalLength ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->GlobalLength();
}

int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getConstMultiVector(selfID)->Stride();
}

boolean Epetra_MultiVector_ConstantStride ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return         CEpetra::getConstMultiVector(selfID)->ConstantStride();
}

int Epetra_MultiVector_ReplaceMap ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID )
{
    return CEpetra::getMultiVector(selfID)->ReplaceMap(
        *CEpetra::getBlockMap(mapID));
}


} // extern "C"


//
// Definitions from CEpetra_MultiVector_Cpp.hpp
//


/* get Epetra_MultiVector from non-const table using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id);
    return tableOfMultiVectors().get(aid);
}

/* get Epetra_MultiVector from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CTrilinos_Object_ID_t id )
{
    return tableOfMultiVectors().get(id);
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CT_Epetra_MultiVector_ID */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_MultiVector_ID_t>(id);
    if (id.is_const) {
        return tableOfConstMultiVectors().get(aid);
    } else {
        return tableOfMultiVectors().get(aid);
    }
}

/* get const Epetra_MultiVector from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_MultiVector>
CEpetra::getConstMultiVector( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstMultiVectors().get(id);
    } else {
        return tableOfMultiVectors().get(id);
    }
}

/* store Epetra_MultiVector in non-const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeMultiVector( Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
            tableOfMultiVectors().storeCopy(pobj));
}

/* store const Epetra_MultiVector in const table */
CT_Epetra_MultiVector_ID_t
CEpetra::storeConstMultiVector( const Epetra_MultiVector *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_MultiVector_ID_t>(
            tableOfConstMultiVectors().storeCopy(pobj));
}

/* dump contents of Epetra_MultiVector and const Epetra_MultiVector tables */
void
CEpetra::purgeMultiVectorTables(  )
{
    tableOfMultiVectors().purge();
    tableOfConstMultiVectors().purge();
}



