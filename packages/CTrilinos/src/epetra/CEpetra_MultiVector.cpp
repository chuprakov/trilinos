#include "CTrilinos_config.h"

#include "CEpetra_BlockMap_Cpp.hpp"
#include "Epetra_DataAccess.h"
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


Table<Epetra_MultiVector>& tableOfMultiVectors()
{
    static Table<Epetra_MultiVector> loc_tableOfMultiVectors(CT_Epetra_MultiVector_ID, "CT_Epetra_MultiVector_ID");
    return loc_tableOfMultiVectors;
}


} // namespace


//
// Definitions from CEpetra_MultiVector.h
//


extern "C" {


CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Cast (
  CTrilinos_Object_ID_t id )
{
    return CTrilinos::cast(tableOfMultiVectors(), id);
}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( 
  CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return tableOfMultiVectors().store(new Epetra_MultiVector(
        *pMap, NumVectors, zeroOut));

}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( 
  CT_Epetra_MultiVector_ID_t SourceID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pSource = CEpetra::getMultiVector(SourceID);

    return tableOfMultiVectors().store(new Epetra_MultiVector(
        *pSource));

}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_From2DA ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, double * A, 
  int MyLDA, int NumVectors )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *pMap, A, MyLDA, NumVectors));

}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromAOP ( 
  Epetra_DataAccess CV, CT_Epetra_BlockMap_ID_t MapID, 
  double ** ArrayOfPointers, int NumVectors )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pMap = CEpetra::getBlockMap(MapID);

    return tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *pMap, ArrayOfPointers, NumVectors));

}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromList ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int * Indices, int NumVectors )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pSource = CEpetra::getMultiVector(SourceID);

    return tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *pSource, Indices, NumVectors));

}

CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create_FromRange ( 
  Epetra_DataAccess CV, CT_Epetra_MultiVector_ID_t SourceID, 
  int StartIndex, int NumVectors )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pSource = CEpetra::getMultiVector(SourceID);

    return tableOfMultiVectors().store(new Epetra_MultiVector(
        CV, *pSource, StartIndex, NumVectors));

}

void Epetra_MultiVector_Destroy ( 
  CT_Epetra_MultiVector_ID_t * selfID )
{
    tableOfMultiVectors().remove(selfID);
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
    return CEpetra::getMultiVector(selfID)->ExtractCopy(A, MyLDA);
}

int Epetra_MultiVector_ExtractCopy_FillAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** ArrayOfPointers )
{
    return CEpetra::getMultiVector(selfID)->ExtractCopy(
        ArrayOfPointers);
}

int Epetra_MultiVector_ExtractView_Set2DA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ** A, int * MyLDA )
{
    return CEpetra::getMultiVector(selfID)->ExtractView(A, MyLDA);
}

int Epetra_MultiVector_ExtractView_SetAOP ( 
  CT_Epetra_MultiVector_ID_t selfID, double *** ArrayOfPointers )
{
    return CEpetra::getMultiVector(selfID)->ExtractView(
        ArrayOfPointers);
}

int Epetra_MultiVector_Dot ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID, 
  double * Result )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);

    return CEpetra::getMultiVector(selfID)->Dot(*pA, Result);
}

int Epetra_MultiVector_Abs ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);

    return CEpetra::getMultiVector(selfID)->Abs(*pA);
}

int Epetra_MultiVector_Reciprocal ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_MultiVector_ID_t AID )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);

    return CEpetra::getMultiVector(selfID)->Reciprocal(*pA);
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
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);

    return CEpetra::getMultiVector(selfID)->Scale(ScalarA, *pA);
}

int Epetra_MultiVector_Update_WithA ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarThis )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);

    return CEpetra::getMultiVector(selfID)->Update(
        ScalarA, *pA, ScalarThis);
}

int Epetra_MultiVector_Update_WithAB ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarA, 
  CT_Epetra_MultiVector_ID_t AID, double ScalarB, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pB = CEpetra::getMultiVector(BID);

    return CEpetra::getMultiVector(selfID)->Update(
        ScalarA, *pA, ScalarB, *pB, ScalarThis);
}

int Epetra_MultiVector_Norm1 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getMultiVector(selfID)->Norm1(Result);
}

int Epetra_MultiVector_Norm2 ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getMultiVector(selfID)->Norm2(Result);
}

int Epetra_MultiVector_NormInf ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getMultiVector(selfID)->NormInf(Result);
}

int Epetra_MultiVector_NormWeighted ( 
  CT_Epetra_MultiVector_ID_t selfID, 
  CT_Epetra_MultiVector_ID_t WeightsID, double * Result )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pWeights = CEpetra::getMultiVector(WeightsID);

    return CEpetra::getMultiVector(selfID)->NormWeighted(
        *pWeights, Result);
}

int Epetra_MultiVector_MinValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getMultiVector(selfID)->MinValue(Result);
}

int Epetra_MultiVector_MaxValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getMultiVector(selfID)->MaxValue(Result);
}

int Epetra_MultiVector_MeanValue ( 
  CT_Epetra_MultiVector_ID_t selfID, double * Result )
{
    return CEpetra::getMultiVector(selfID)->MeanValue(Result);
}

int Epetra_MultiVector_Multiply_Matrix ( 
  CT_Epetra_MultiVector_ID_t selfID, char TransA, char TransB, 
  double ScalarAB, CT_Epetra_MultiVector_ID_t AID, 
  CT_Epetra_MultiVector_ID_t BID, double ScalarThis )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pB = CEpetra::getMultiVector(BID);

    return CEpetra::getMultiVector(selfID)->Multiply(
        TransA, TransB, ScalarAB, *pA, *pB, ScalarThis);
}

int Epetra_MultiVector_Multiply_ByEl ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pB = CEpetra::getMultiVector(BID);

    return CEpetra::getMultiVector(selfID)->Multiply(
        ScalarAB, *pA, *pB, ScalarThis);
}

int Epetra_MultiVector_ReciprocalMultiply ( 
  CT_Epetra_MultiVector_ID_t selfID, double ScalarAB, 
  CT_Epetra_MultiVector_ID_t AID, CT_Epetra_MultiVector_ID_t BID, 
  double ScalarThis )
{
    const Teuchos::RCP<Epetra_MultiVector> 
        pA = CEpetra::getMultiVector(AID);
    const Teuchos::RCP<Epetra_MultiVector> 
        pB = CEpetra::getMultiVector(BID);

    return CEpetra::getMultiVector(selfID)->ReciprocalMultiply(
        ScalarAB, *pA, *pB, ScalarThis);
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

    const Teuchos::RCP<Epetra_MultiVector> 
        pSource = CEpetra::getMultiVector(SourceID);

    self = *pSource;
}

double * Epetra_MultiVector_getArray ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    Epetra_MultiVector& self = *( CEpetra::getMultiVector(selfID) );

    return self[i];
}

CT_Epetra_Vector_ID_t Epetra_MultiVector_getVector ( 
  CT_Epetra_MultiVector_ID_t selfID, int i )
{
    Epetra_MultiVector& self = *( CEpetra::getMultiVector(selfID) );

    return CEpetra::storeVector(self(i));
}

int Epetra_MultiVector_NumVectors ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->NumVectors();
}

int Epetra_MultiVector_MyLength ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->MyLength();
}

int Epetra_MultiVector_GlobalLength ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->GlobalLength();
}

int Epetra_MultiVector_Stride ( CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->Stride();
}

boolean Epetra_MultiVector_ConstantStride ( 
  CT_Epetra_MultiVector_ID_t selfID )
{
    return CEpetra::getMultiVector(selfID)->ConstantStride();
}

int Epetra_MultiVector_ReplaceMap ( 
  CT_Epetra_MultiVector_ID_t selfID, CT_Epetra_BlockMap_ID_t mapID )
{
    const Teuchos::RCP<Epetra_BlockMap> 
        pmap = CEpetra::getBlockMap(mapID);

    return CEpetra::getMultiVector(selfID)->ReplaceMap(*pmap);
}


} // extern "C"


//
// Definitions from CEpetra_MultiVector_Cpp.hpp
//


const Teuchos::RCP<Epetra_MultiVector>
CEpetra::getMultiVector( CT_Epetra_MultiVector_ID_t id )
{
    return tableOfMultiVectors().get(id);
}

CT_Epetra_MultiVector_ID_t
CEpetra::storeMultiVector( const Epetra_MultiVector *pobj )
{
    return tableOfMultiVectors().storeCopy(pobj);
}

void
CEpetra::purgeMultiVectorTable(  )
{
    tableOfMultiVectors().purge();
}



