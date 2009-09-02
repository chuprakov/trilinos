#include "CTrilinos_utils.hpp"
#include "CTrilinos_utils_templ.hpp"


#include "CEpetra_Distributor_Cpp.hpp"
#include "CEpetra_SerialComm_Cpp.hpp"
#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_Comm_Cpp.hpp"
#include "CEpetra_Operator_Cpp.hpp"
#include "CEpetra_MultiVector_Cpp.hpp"
#include "CEpetra_OffsetIndex_Cpp.hpp"
#include "CEpetra_Object_Cpp.hpp"
#include "CEpetra_RowMatrix_Cpp.hpp"
#include "CEpetra_CompObject_Cpp.hpp"
#include "CEpetra_Directory_Cpp.hpp"
#include "CEpetra_Flops_Cpp.hpp"
#include "CEpetra_SrcDistObject_Cpp.hpp"
#include "CEpetra_MpiComm_Cpp.hpp"
#include "CEpetra_CrsMatrix_Cpp.hpp"
#include "CEpetra_CrsGraph_Cpp.hpp"
#include "CEpetra_DistObject_Cpp.hpp"
#include "CEpetra_Vector_Cpp.hpp"
#include "CEpetra_Export_Cpp.hpp"
#include "CEpetra_Map_Cpp.hpp"
#include "CEpetra_BlockMap_Cpp.hpp"
#include "CEpetra_Import_Cpp.hpp"
#include "CEpetra_Time_Cpp.hpp"
#include "CEpetra_JadMatrix_Cpp.hpp"
#include "CEpetra_LinearProblem_Cpp.hpp"
#include "CEpetra_LAPACK_Cpp.hpp"
#include "CTeuchos_CommandLineProcessor_Cpp.hpp"
#include "CTeuchos_ParameterList_Cpp.hpp"
#include "CTeuchos_ParameterEntry_Cpp.hpp"
#include "CTeuchos_any_Cpp.hpp"
#include "CAmesos_BaseSolver_Cpp.hpp"
#include "CAmesos_Cpp.hpp"


namespace CTrilinos {


void pass_string_out( const std::string * const s, char *c[] )
{
    *c = new char[s->size()+1];
    strcpy(*c, s->c_str());
}


void pass_string_in( const char * const c[], std::string *s )
{
    s = new std::string(*c);
}


/* stringify the enum name */
std::string
enum2str( CTrilinos_Type_ID_t ty )
{
    std::string s;
    switch (ty) {
    case CT_Epetra_Distributor_ID:
        s.assign("CT_Epetra_Distributor_ID");
        break;
    case CT_Epetra_SerialComm_ID:
        s.assign("CT_Epetra_SerialComm_ID");
        break;
    case CT_Epetra_BLAS_ID:
        s.assign("CT_Epetra_BLAS_ID");
        break;
    case CT_Epetra_Comm_ID:
        s.assign("CT_Epetra_Comm_ID");
        break;
    case CT_Epetra_Operator_ID:
        s.assign("CT_Epetra_Operator_ID");
        break;
    case CT_Epetra_MultiVector_ID:
        s.assign("CT_Epetra_MultiVector_ID");
        break;
    case CT_Epetra_OffsetIndex_ID:
        s.assign("CT_Epetra_OffsetIndex_ID");
        break;
    case CT_Epetra_Object_ID:
        s.assign("CT_Epetra_Object_ID");
        break;
    case CT_Epetra_RowMatrix_ID:
        s.assign("CT_Epetra_RowMatrix_ID");
        break;
    case CT_Epetra_CompObject_ID:
        s.assign("CT_Epetra_CompObject_ID");
        break;
    case CT_Epetra_Directory_ID:
        s.assign("CT_Epetra_Directory_ID");
        break;
    case CT_Epetra_Flops_ID:
        s.assign("CT_Epetra_Flops_ID");
        break;
    case CT_Epetra_SrcDistObject_ID:
        s.assign("CT_Epetra_SrcDistObject_ID");
        break;
    case CT_Epetra_MpiComm_ID:
        s.assign("CT_Epetra_MpiComm_ID");
        break;
    case CT_Epetra_CrsMatrix_ID:
        s.assign("CT_Epetra_CrsMatrix_ID");
        break;
    case CT_Epetra_CrsGraph_ID:
        s.assign("CT_Epetra_CrsGraph_ID");
        break;
    case CT_Epetra_DistObject_ID:
        s.assign("CT_Epetra_DistObject_ID");
        break;
    case CT_Epetra_Vector_ID:
        s.assign("CT_Epetra_Vector_ID");
        break;
    case CT_Epetra_Export_ID:
        s.assign("CT_Epetra_Export_ID");
        break;
    case CT_Epetra_Map_ID:
        s.assign("CT_Epetra_Map_ID");
        break;
    case CT_Epetra_BlockMap_ID:
        s.assign("CT_Epetra_BlockMap_ID");
        break;
    case CT_Epetra_Import_ID:
        s.assign("CT_Epetra_Import_ID");
        break;
    case CT_Epetra_Time_ID:
        s.assign("CT_Epetra_Time_ID");
        break;
    case CT_Epetra_JadMatrix_ID:
        s.assign("CT_Epetra_JadMatrix_ID");
        break;
    case CT_Epetra_LinearProblem_ID:
        s.assign("CT_Epetra_LinearProblem_ID");
        break;
    case CT_Epetra_LAPACK_ID:
        s.assign("CT_Epetra_LAPACK_ID");
        break;
    case CT_Teuchos_CommandLineProcessor_ID:
        s.assign("CT_Teuchos_CommandLineProcessor_ID");
        break;
    case CT_Teuchos_ParameterList_ID:
        s.assign("CT_Teuchos_ParameterList_ID");
        break;
    case CT_Teuchos_ParameterEntry_ID:
        s.assign("CT_Teuchos_ParameterEntry_ID");
        break;
    case CT_Teuchos_any_ID:
        s.assign("CT_Teuchos_any_ID");
        break;
    case CT_Amesos_BaseSolver_ID:
        s.assign("CT_Amesos_BaseSolver_ID");
        break;
    case CT_Amesos_ID:
        s.assign("CT_Amesos_ID");
        break;
    default:
        s.assign("(unrecognized)");
        break;
    }

    return s;
}


/* isSameObject(generic_id, generic_id) */
bool
isSameObject( CTrilinos_Object_ID_t id1, CTrilinos_Object_ID_t id2 )
{
    bool shares = false;

    if (id1.is_const) {
        switch (id1.type) {
        case CT_Epetra_Distributor_ID:
            shares = isSameObject(CEpetra::getConstDistributor(id1), id2);
            break;
        case CT_Epetra_SerialComm_ID:
            shares = isSameObject(CEpetra::getConstSerialComm(id1), id2);
            break;
        case CT_Epetra_BLAS_ID:
            shares = isSameObject(CEpetra::getConstBLAS(id1), id2);
            break;
        case CT_Epetra_Comm_ID:
            shares = isSameObject(CEpetra::getConstComm(id1), id2);
            break;
        case CT_Epetra_Operator_ID:
            shares = isSameObject(CEpetra::getConstOperator(id1), id2);
            break;
        case CT_Epetra_MultiVector_ID:
            shares = isSameObject(CEpetra::getConstMultiVector(id1), id2);
            break;
        case CT_Epetra_OffsetIndex_ID:
            shares = isSameObject(CEpetra::getConstOffsetIndex(id1), id2);
            break;
        case CT_Epetra_Object_ID:
            shares = isSameObject(CEpetra::getConstObject(id1), id2);
            break;
        case CT_Epetra_RowMatrix_ID:
            shares = isSameObject(CEpetra::getConstRowMatrix(id1), id2);
            break;
        case CT_Epetra_CompObject_ID:
            shares = isSameObject(CEpetra::getConstCompObject(id1), id2);
            break;
        case CT_Epetra_Directory_ID:
            shares = isSameObject(CEpetra::getConstDirectory(id1), id2);
            break;
        case CT_Epetra_Flops_ID:
            shares = isSameObject(CEpetra::getConstFlops(id1), id2);
            break;
        case CT_Epetra_SrcDistObject_ID:
            shares = isSameObject(CEpetra::getConstSrcDistObject(id1), id2);
            break;
#ifdef HAVE_MPI
        case CT_Epetra_MpiComm_ID:
            shares = isSameObject(CEpetra::getConstMpiComm(id1), id2);
            break;
#endif /* HAVE_MPI */
        case CT_Epetra_CrsMatrix_ID:
            shares = isSameObject(CEpetra::getConstCrsMatrix(id1), id2);
            break;
        case CT_Epetra_CrsGraph_ID:
            shares = isSameObject(CEpetra::getConstCrsGraph(id1), id2);
            break;
        case CT_Epetra_DistObject_ID:
            shares = isSameObject(CEpetra::getConstDistObject(id1), id2);
            break;
        case CT_Epetra_Vector_ID:
            shares = isSameObject(CEpetra::getConstVector(id1), id2);
            break;
        case CT_Epetra_Export_ID:
            shares = isSameObject(CEpetra::getConstExport(id1), id2);
            break;
        case CT_Epetra_Map_ID:
            shares = isSameObject(CEpetra::getConstMap(id1), id2);
            break;
        case CT_Epetra_BlockMap_ID:
            shares = isSameObject(CEpetra::getConstBlockMap(id1), id2);
            break;
        case CT_Epetra_Import_ID:
            shares = isSameObject(CEpetra::getConstImport(id1), id2);
            break;
        case CT_Epetra_Time_ID:
            shares = isSameObject(CEpetra::getConstTime(id1), id2);
            break;
        case CT_Epetra_JadMatrix_ID:
            shares = isSameObject(CEpetra::getConstJadMatrix(id1), id2);
            break;
        case CT_Epetra_LinearProblem_ID:
            shares = isSameObject(CEpetra::getConstLinearProblem(id1), id2);
            break;
        case CT_Epetra_LAPACK_ID:
            shares = isSameObject(CEpetra::getConstLAPACK(id1), id2);
            break;
        case CT_Teuchos_CommandLineProcessor_ID:
            shares = isSameObject(CTeuchos::getConstCommandLineProcessor(id1), id2);
            break;
        case CT_Teuchos_ParameterList_ID:
            shares = isSameObject(CTeuchos::getConstParameterList(id1), id2);
            break;
        case CT_Teuchos_ParameterEntry_ID:
            shares = isSameObject(CTeuchos::getConstParameterEntry(id1), id2);
            break;
        case CT_Teuchos_any_ID:
            shares = isSameObject(CTeuchos::getConstany(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            shares = isSameObject(CAmesos::getConstBaseSolver(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_ID:
            shares = isSameObject(CAmesos::getConstAmesos(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        default:
            break;
        }
    } else {
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
        case CT_Epetra_Time_ID:
            shares = isSameObject(CEpetra::getTime(id1), id2);
            break;
        case CT_Epetra_JadMatrix_ID:
            shares = isSameObject(CEpetra::getJadMatrix(id1), id2);
            break;
        case CT_Epetra_LinearProblem_ID:
            shares = isSameObject(CEpetra::getLinearProblem(id1), id2);
            break;
        case CT_Epetra_LAPACK_ID:
            shares = isSameObject(CEpetra::getLAPACK(id1), id2);
            break;
        case CT_Teuchos_CommandLineProcessor_ID:
            shares = isSameObject(CTeuchos::getCommandLineProcessor(id1), id2);
            break;
        case CT_Teuchos_ParameterList_ID:
            shares = isSameObject(CTeuchos::getParameterList(id1), id2);
            break;
        case CT_Teuchos_ParameterEntry_ID:
            shares = isSameObject(CTeuchos::getParameterEntry(id1), id2);
            break;
        case CT_Teuchos_any_ID:
            shares = isSameObject(CTeuchos::getany(id1), id2);
            break;
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_BaseSolver_ID:
            shares = isSameObject(CAmesos::getBaseSolver(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
#ifdef HAVE_CTRILINOS_AMESOS
        case CT_Amesos_ID:
            shares = isSameObject(CAmesos::getAmesos(id1), id2);
            break;
#endif /* HAVE_CTRILINOS_AMESOS */
        default:
            break;
        }
    }

    return shares;
}


} // namespace CTrilinos

