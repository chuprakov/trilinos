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


Table<Epetra_BLAS>& tableOfBLASs()
{
    static Table<Epetra_BLAS>
        loc_tableOfBLASs(CT_Epetra_BLAS_ID, "CT_Epetra_BLAS_ID");
    return loc_tableOfBLASs;
}


} // namespace


//
// Definitions from CEpetra_BLAS.h
//


extern "C" {


CT_Epetra_BLAS_ID_t Epetra_BLAS_Cast ( CTrilinos_Object_ID_t id )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        CTrilinos::cast(tableOfBLASs(), id));
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
    const Teuchos::RCP<Epetra_BLAS> 
        pBLAS = CEpetra::getBLAS(BLASID);

    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tableOfBLASs().store(new Epetra_BLAS(
        *pBLAS)));

}

void Epetra_BLAS_Destroy ( CT_Epetra_BLAS_ID_t * selfID )
{
    CTrilinos_Object_ID_t id =
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(*selfID);
    tableOfBLASs().remove(&id);
    *selfID = CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(id);
}

float Epetra_BLAS_ASUM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getBLAS(selfID)->ASUM(N, X, INCX);
}

double Epetra_BLAS_ASUM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getBLAS(selfID)->ASUM(N, X, INCX);
}

float Epetra_BLAS_DOT_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const float * Y, const int INCX, const int INCY )
{
    return CEpetra::getBLAS(selfID)->DOT(N, X, Y, INCX, INCY);
}

double Epetra_BLAS_DOT_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const double * Y, const int INCX, const int INCY )
{
    return CEpetra::getBLAS(selfID)->DOT(N, X, Y, INCX, INCY);
}

float Epetra_BLAS_NRM2_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getBLAS(selfID)->NRM2(N, X, INCX);
}

double Epetra_BLAS_NRM2_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getBLAS(selfID)->NRM2(N, X, INCX);
}

void Epetra_BLAS_SCAL_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  float * X, const int INCX )
{
    CEpetra::getBLAS(selfID)->SCAL(N, ALPHA, X, INCX);
}

void Epetra_BLAS_SCAL_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  double * X, const int INCX )
{
    CEpetra::getBLAS(selfID)->SCAL(N, ALPHA, X, INCX);
}

void Epetra_BLAS_COPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, float * Y, 
  const int INCX, const int INCY )
{
    CEpetra::getBLAS(selfID)->COPY(N, X, Y, INCX, INCY);
}

void Epetra_BLAS_COPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  double * Y, const int INCX, const int INCY )
{
    CEpetra::getBLAS(selfID)->COPY(N, X, Y, INCX, INCY);
}

int Epetra_BLAS_IAMAX_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX )
{
    return CEpetra::getBLAS(selfID)->IAMAX(N, X, INCX);
}

int Epetra_BLAS_IAMAX_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX )
{
    return CEpetra::getBLAS(selfID)->IAMAX(N, X, INCX);
}

void Epetra_BLAS_AXPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  const float * X, float * Y, const int INCX, const int INCY )
{
    CEpetra::getBLAS(selfID)->AXPY(N, ALPHA, X, Y, INCX, INCY);
}

void Epetra_BLAS_AXPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  const double * X, double * Y, const int INCX, const int INCY )
{
    CEpetra::getBLAS(selfID)->AXPY(N, ALPHA, X, Y, INCX, INCY);
}

void Epetra_BLAS_GEMV_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const float ALPHA, const float * A, const int LDA, 
  const float * X, const float BETA, float * Y, const int INCX, 
  const int INCY )
{
    CEpetra::getBLAS(selfID)->GEMV(
        TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, INCX, INCY);
}

void Epetra_BLAS_GEMV_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const double ALPHA, const double * A, const int LDA, 
  const double * X, const double BETA, double * Y, const int INCX, 
  const int INCY )
{
    CEpetra::getBLAS(selfID)->GEMV(
        TRANS, M, N, ALPHA, A, LDA, X, BETA, Y, INCX, INCY);
}

void Epetra_BLAS_GEMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const float ALPHA, 
  const float * A, const int LDA, const float * B, const int LDB, 
  const float BETA, float * C, const int LDC )
{
    CEpetra::getBLAS(selfID)->GEMM(
        TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_GEMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const double ALPHA, 
  const double * A, const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC )
{
    CEpetra::getBLAS(selfID)->GEMM(
        TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_SYMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const float ALPHA, const float * A, 
  const int LDA, const float * B, const int LDB, const float BETA, 
  float * C, const int LDC )
{
    CEpetra::getBLAS(selfID)->SYMM(
        SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_SYMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const double ALPHA, const double * A, 
  const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC )
{
    CEpetra::getBLAS(selfID)->SYMM(
        SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB, BETA, C, LDC);
}

void Epetra_BLAS_TRMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const float ALPHA, const float * A, const int LDA, float * B, 
  const int LDB )
{
    CEpetra::getBLAS(selfID)->TRMM(
        SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
}

void Epetra_BLAS_TRMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const double ALPHA, const double * A, const int LDA, double * B, 
  const int LDB )
{
    CEpetra::getBLAS(selfID)->TRMM(
        SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA, B, LDB);
}


} // extern "C"


//
// Definitions from CEpetra_BLAS_Cpp.hpp
//


const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CT_Epetra_BLAS_ID_t id )
{
    return tableOfBLASs().get(
        CTrilinos::abstractType<CT_Epetra_BLAS_ID_t>(id));
}

const Teuchos::RCP<Epetra_BLAS>
CEpetra::getBLAS( CTrilinos_Object_ID_t id )
{
    return tableOfBLASs().get(id);
}

CT_Epetra_BLAS_ID_t
CEpetra::storeBLAS( const Epetra_BLAS *pobj )
{
    return CTrilinos::concreteType<CT_Epetra_BLAS_ID_t>(
        tableOfBLASs().storeCopy(pobj));
}

void
CEpetra::purgeBLASTable(  )
{
    tableOfBLASs().purge();
}



