#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"

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
    case CT_Epetra_Comm_ID:
        newid = destTable.cast(CEpetra::getComm(id));
        break;
    case CT_Epetra_Map_ID:
        newid = destTable.cast(CEpetra::getMap(id));
        break;
    case CT_Epetra_MultiVector_ID:
        newid = destTable.cast(CEpetra::getMultiVector(id));
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        newid = destTable.cast(CEpetra::getMpiComm(id));
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_SerialComm_ID:
        newid = destTable.cast(CEpetra::getSerialComm(id));
        break;
    case CT_Epetra_BlockMap_ID:
        newid = destTable.cast(CEpetra::getBlockMap(id));
        break;
    case CT_Epetra_Vector_ID:
        newid = destTable.cast(CEpetra::getVector(id));
        break;
    default:
        break;
    }

    return newid;
}


} // namespace CTrilinos


#endif

