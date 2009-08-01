#include "CTrilinos_config.h"

#include "CEpetra_BLAS_Cpp.hpp"
#include "CEpetra_BLAS.h"
#include "Epetra_BLAS.h"
#include "Teuchos_RCP.hpp"
#include "CTrilinos_enums.h"
#include "CTrilinos_utils.hpp"
#include "CTrilinos_Table.hpp"


namespace {


using Teuchos::RCP;
using CTrilinos::Table;


/* table to hold objects of type Epetra_BLAS */
Table<Epetra_BLAS>& tableOfBLASs()
{
    static Table<Epetra_BLAS>
        loc_tableOfBLASs(CT_Epetra_BLAS_ID, "CT_Epetra_BLAS_ID", false);
    return loc_tableOfBLASs;
}

/* table to hold objects of type const Epetra_BLAS */
Table<const Epetra_BLAS>& tableOfConstBLASs()
{
    static Table<const Epetra_BLAS>
        loc_tableOfConstBLASs(CT_Epetra_BLAS_ID, "CT_Epetra_BLAS_ID", true);
    return loc_tableOfConstBLASs;
}


} // namespace


//
// Definitions from CEpetra_BLAS.h
//


extern "C" {


CT_Epetra_BLAS_ID_t Epetra_BLAS_Cast ( CTrilinos_Object_ID_t id )
{
    CTrilinos_Object_ID_t newid;
    if (id.is_const) {
        newid = CTrilinos::cast(tableOfConstBLASs(), id);
    } else {
        newid = CTrilinos::cast(tableOfBLASs(), id);
    }
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(newid);
}

CTrilinos_Object_ID_t Epetra_BLAS_Abstract ( 
  CT_Epetra_BLAS_ID_t id )
{
    return CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id);
}

CT_Epetra_BLAS_ID_t Epetra_BLAS_Create (  )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tableOfBLASs().store(new Epetra_BLAS()));
}

CT_Epetra_BLAS_ID_t Epetra_BLAS_Duplicate ( 
  CT_Epetra_BLAS_ID_t BLASID )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tableOfBLASs().store(new Epetra_BLAS(
        *CEpetra::getBLAS(BLASID))));
}

void Epetra_BLAS_Destroy ( CT_Epetra_BLAS_ID_t * selfID )
{
    CTrilinos_Object_ID_t aid
        = CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(*selfID);
    if (selfID->is_const) {
        tableOfConstBLASs().remove(&aid);
    } else {
        tableOfBLASs().remove(&aid);
    }
    *selfID = CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(aid);
}

float Epetra_BLAS_ASUM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->ASUM(N, X, INCX);
}

double Epetra_BLAS_ASUM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->ASUM(N, X, INCX);
}

float Epetra_BLAS_DOT_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const float * Y, const int INCX, const int INCY )
{
    return CEpetra::getConstBLAS(selfID)->DOT(N, X, Y, INCX, INCY);
}

double Epetra_BLAS_DOT_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const double * Y, const int INCX, const int INCY )
{
    return CEpetra::getConstBLAS(selfID)->DOT(N, X, Y, INCX, INCY);
}

float Epetra_BLAS_NRM2_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->NRM2(N, X, INCX);
}

double Epetra_BLAS_NRM2_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->NRM2(N, X, INCX);
}

void Epetra_BLAS_SCAL_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  float * X, const int INCX )
{
    CEpetra::getConstBLAS(selfID)->SCAL(N, ALPHA, X, INCX);
}

void Epetra_BLAS_SCAL_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  double * X, const int INCX )
{
    CEpetra::getConstBLAS(selfID)->SCAL(N, ALPHA, X, INCX);
}

void Epetra_BLAS_COPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, float * Y, 
  const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->COPY(N, X, Y, INCX, INCY);
}

void Epetra_BLAS_COPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  double * Y, const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->COPY(N, X, Y, INCX, INCY);
}

int Epetra_BLAS_IAMAX_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->IAMAX(N, X, INCX);
}

int Epetra_BLAS_IAMAX_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getConstBLAS(selfID)->IAMAX(N, X, INCX);
}

void Epetra_BLAS_AXPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  const float * X, float * Y, const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->AXPY(N, ALPHA, X, Y, INCX, INCY);
}

void Epetra_BLAS_AXPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  const double * X, double * Y, const int INCX, const int INCY )
{
    CEpetra::getConstBLAS(selfID)->AXPY(N, ALPHA, X, Y, INCX, INCY);
}

void Epetra_BLAS_GEMV_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const float ALPHA, const float * A, const int LDA, 
  const float * X, const float BETA, float * Y, const int INCX, 
  const int INCY )
{
    CEpetra::getConstBLAS(selfID)->GEMV(
        TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, INCX, INCY);
}

void Epetra_BLAS_GEMV_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const double ALPHA, const double * A, const int LDA, 
  const double * X, const double BETA, double * Y, const int INCX, 
  const int INCY )
{
    CEpetra::getConstBLAS(selfID)->GEMV(
        TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, INCX, INCY);
}

void Epetra_BLAS_GEMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const float ALPHA, 
  const float * A, const int LDA, const float * B, const int LDB, 
  const float BETA, float * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->GEMM(
        TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_GEMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const double ALPHA, 
  const double * A, const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->GEMM(
        TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_SYMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const float ALPHA, const float * A, 
  const int LDA, const float * B, const int LDB, const float BETA, 
  float * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->SYMM(
        SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_SYMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const double ALPHA, const double * A, 
  const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC )
{
    CEpetra::getConstBLAS(selfID)->SYMM(
        SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_TRMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const float ALPHA, const float * A, const int LDA, float * B, 
  const int LDB )
{
    CEpetra::getConstBLAS(selfID)->TRMM(
        SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
}

void Epetra_BLAS_TRMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const double ALPHA, const double * A, const int LDA, double * B, 
  const int LDB )
{
    CEpetra::getConstBLAS(selfID)->TRMM(
        SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
}


} // extern "C"


//
// Definitions from CEpetra_BLAS_Cpp.hpp
//


/* get Epetra_BLAS from non-const table using CT_Epetra_BLAS_ID */
const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CT_Epetra_BLAS_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id);
    return tableOfBLASs().get(aid);
}

/* get Epetra_BLAS from non-const table using CTrilinos_Object_ID_t */
const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CTrilinos_Object_ID_t id )
{
    return tableOfBLASs().get(id);
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CT_Epetra_BLAS_ID */
const Teuchos::RCP<const Epetra_BLAS>
CEpetra::getConstBLAS( CT_Epetra_BLAS_ID_t id )
{
    CTrilinos_Object_ID_t aid
            = CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id);
    if (id.is_const) {
        return tableOfConstBLASs().get(aid);
    } else {
        return tableOfBLASs().get(aid);
    }
}

/* get const Epetra_BLAS from either the const or non-const table
 * using CTrilinos_Object_ID_t */
const Teuchos::RCP<const Epetra_BLAS>
CEpetra::getConstBLAS( CTrilinos_Object_ID_t id )
{
    if (id.is_const) {
        return tableOfConstBLASs().get(id);
    } else {
        return tableOfBLASs().get(id);
    }
}

/* store Epetra_BLAS in non-const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeBLAS( Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
            tableOfBLASs().storeCopy(pobj));
}

/* store const Epetra_BLAS in const table */
CT_Epetra_BLAS_ID_t
CEpetra::storeConstBLAS( const Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
            tableOfConstBLASs().storeCopy(pobj));
}

/* dump contents of Epetra_BLAS and const Epetra_BLAS tables */
void
CEpetra::purgeBLASTables(  )
{
    tableOfBLASs().purge();
    tableOfConstBLASs().purge();
}



