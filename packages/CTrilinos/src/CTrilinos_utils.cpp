#include "CTrilinos_utils.hpp"


namespace CTrilinos {


/* convert struct from specific type to generic CTrilinos_Object_ID_t
 * but keep the content in tact */
const char *str2cptr( std::string s )
{
    char *pc = NULL;
    pc = new char [s.size()+1];
    strcpy(pc, s.c_str());
    return pc;
}


/* convert struct from generic CTrilinos_Object_ID_t to specific type
 * but keep the content in tact */
std::string cptr2str( const char *pc )
{
    std::string s(pc);
    return s;
}


/* isSameObject(generic_id, generic_id) */
bool
isSameObject( CTrilinos_Object_ID_t id1, CTrilinos_Object_ID_t id2 )
{
    bool shares = false;

    switch (id1.type) {
    case CT_Epetra_Distributor_ID:
        shares = isSameObject(CEpetra::getDistributor(id1), id2);
        break;
    case CT_Epetra_SerialComm_ID:
        shares = isSameObject(CEpetra::getSerialComm(id1), id2);
        break;
    case CT_Epetra_BLAS_ID:
        shares = isSameObject(CEpetra::getBLAS(id1), id2);
        break;
    case CT_Epetra_Comm_ID:
        shares = isSameObject(CEpetra::getComm(id1), id2);
        break;
    case CT_Epetra_Operator_ID:
        shares = isSameObject(CEpetra::getOperator(id1), id2);
        break;
    case CT_Epetra_MultiVector_ID:
        shares = isSameObject(CEpetra::getMultiVector(id1), id2);
        break;
    case CT_Epetra_OffsetIndex_ID:
        shares = isSameObject(CEpetra::getOffsetIndex(id1), id2);
        break;
    case CT_Epetra_Object_ID:
        shares = isSameObject(CEpetra::getObject(id1), id2);
        break;
    case CT_Epetra_RowMatrix_ID:
        shares = isSameObject(CEpetra::getRowMatrix(id1), id2);
        break;
    case CT_Epetra_CompObject_ID:
        shares = isSameObject(CEpetra::getCompObject(id1), id2);
        break;
    case CT_Epetra_Directory_ID:
        shares = isSameObject(CEpetra::getDirectory(id1), id2);
        break;
    case CT_Epetra_Flops_ID:
        shares = isSameObject(CEpetra::getFlops(id1), id2);
        break;
    case CT_Epetra_SrcDistObject_ID:
        shares = isSameObject(CEpetra::getSrcDistObject(id1), id2);
        break;
#ifdef HAVE_MPI
    case CT_Epetra_MpiComm_ID:
        shares = isSameObject(CEpetra::getMpiComm(id1), id2);
        break;
#endif /* HAVE_MPI */
    case CT_Epetra_CrsMatrix_ID:
        shares = isSameObject(CEpetra::getCrsMatrix(id1), id2);
        break;
    case CT_Epetra_CrsGraph_ID:
        shares = isSameObject(CEpetra::getCrsGraph(id1), id2);
        break;
    case CT_Epetra_DistObject_ID:
        shares = isSameObject(CEpetra::getDistObject(id1), id2);
        break;
    case CT_Epetra_Vector_ID:
        shares = isSameObject(CEpetra::getVector(id1), id2);
        break;
    case CT_Epetra_Export_ID:
        shares = isSameObject(CEpetra::getExport(id1), id2);
        break;
    case CT_Epetra_Map_ID:
        shares = isSameObject(CEpetra::getMap(id1), id2);
        break;
    case CT_Epetra_BlockMap_ID:
        shares = isSameObject(CEpetra::getBlockMap(id1), id2);
        break;
    case CT_Epetra_Import_ID:
        shares = isSameObject(CEpetra::getImport(id1), id2);
        break;
    default:
        break;
    }

    return shares;
}


} // namespace CTrilinos

