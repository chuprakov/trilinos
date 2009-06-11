#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_SerialCommData_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_Data_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_BlockMapData_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "CEpetra_MpiCommData_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_CrsGraphData_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"

#include "CTrilinos_enums.h"
#include "CTrilinos_Table.hpp"


#ifndef CTRILINOS_UTILS_HPP
#define CTRILINOS_UTILS_HPP


namespace CTrilinos {


template <class T>
CTrilinos_Object_ID_t
cast( Table<T> &destTable, CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    newid.type = CT_Invalid_ID;
    newid.index = -1;

    switch (id.type) {
    case CT_Epetra_Distributor_ID:
        newid = destTable.cast(CEpetra::getDistributor(id));
        break;
    case CT_Epetra_SerialCommData_ID:
        newid = destTable.cast(CEpetra::getSerialCommData(id));
        break;
    case CT_Epetra_SerialComm_ID:
        newid = destTable.cast(CEpetra::getSerialComm(id));
        break;
    case CT_Epetra_BLAS_ID:
        newid = destTable.cast(CEpetra::getBLAS(id));
        break;
    case CT_Epetra_Comm_ID:
        newid = destTable.cast(CEpetra::getComm(id));
        break;
    case CT_Epetra_Operator_ID:
        newid = destTable.cast(CEpetra::getOperator(id));
        break;
    case CT_Epetra_MultiVector_ID:
        newid = destTable.cast(CEpetra::getMultiVector(id));
        break;
    case CT_Epetra_OffsetIndex_ID:
        newid = destTable.cast(CEpetra::getOffsetIndex(id));
        break;
    case CT_Epetra_Object_ID:
        newid = destTable.cast(CEpetra::getObject(id));
        break;
    case CT_Epetra_Data_ID:
        newid = destTable.cast(CEpetra::getData(id));
        break;
    case CT_Epetra_RowMatrix_ID:
        newid = destTable.cast(CEpetra::getRowMatrix(id));
        break;
    case CT_Epetra_CompObject_ID:
        newid = destTable.cast(CEpetra::getCompObject(id));
        break;
    case CT_Epetra_Directory_ID:
        newid = destTable.cast(CEpetra::getDirectory(id));
        break;
    case CT_Epetra_BlockMapData_ID:
        newid = destTable.cast(CEpetra::getBlockMapData(id));
        break;
    case CT_Epetra_Flops_ID:
        newid = destTable.cast(CEpetra::getFlops(id));
        break;
    case CT_Epetra_SrcDistObject_ID:
        newid = destTable.cast(CEpetra::getSrcDistObject(id));
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        newid = destTable.cast(CEpetra::getMpiComm(id));
        break;
#endif /* HAVE_MPI */
#ifdef HAVE_MPI
    case CT_Epetra_MpiCommData_ID:
        newid = destTable.cast(CEpetra::getMpiCommData(id));
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        newid = destTable.cast(CEpetra::getCrsMatrix(id));
        break;
    case CT_Epetra_CrsGraph_ID:
        newid = destTable.cast(CEpetra::getCrsGraph(id));
        break;
    case CT_Epetra_DistObject_ID:
        newid = destTable.cast(CEpetra::getDistObject(id));
        break;
    case CT_Epetra_Export_ID:
        newid = destTable.cast(CEpetra::getExport(id));
        break;
    case CT_Epetra_Vector_ID:
        newid = destTable.cast(CEpetra::getVector(id));
        break;
    case CT_Epetra_Map_ID:
        newid = destTable.cast(CEpetra::getMap(id));
        break;
    case CT_Epetra_CrsGraphData_ID:
        newid = destTable.cast(CEpetra::getCrsGraphData(id));
        break;
    case CT_Epetra_BlockMap_ID:
        newid = destTable.cast(CEpetra::getBlockMap(id));
        break;
    case CT_Epetra_Import_ID:
        newid = destTable.cast(CEpetra::getImport(id));
        break;
    default:
        break;
    }

    return newid;
}


} // namespace CTrilinos


#endif

