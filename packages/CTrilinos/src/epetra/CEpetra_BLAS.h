#include "CTrilinos_config.h"


/* True C header file! */


#ifndef CEPETRA_BLAS_H
#define CEPETRA_BLAS_H


#include "CTrilinos_enums.h"


#ifdef __cplusplus
extern "C" {
#endif


CT_Epetra_BLAS_ID_t Epetra_BLAS_Cast ( CTrilinos_Object_ID_t id );

/* Original C++ prototype:
   Epetra_BLAS(void);
*/
CT_Epetra_BLAS_ID_t Epetra_BLAS_Create (  );

/* Original C++ prototype:
   Epetra_BLAS(const Epetra_BLAS& BLAS);
*/
CT_Epetra_BLAS_ID_t Epetra_BLAS_Duplicate ( 
  CT_Epetra_BLAS_ID_t BLASID );

/* Original C++ prototype:
   virtual ~Epetra_BLAS(void);
*/
void Epetra_BLAS_Destroy ( CT_Epetra_BLAS_ID_t * selfID );

/* Original C++ prototype:
   float ASUM(const int N, const float * X, const int INCX = 1) const;
*/
float Epetra_BLAS_ASUM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX );

/* Original C++ prototype:
   double ASUM(const int N, const double * X, const int INCX = 1) const;
*/
double Epetra_BLAS_ASUM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX );

/* Original C++ prototype:
   float DOT(const int N, const float * X, const float * Y, const int INCX = 1, const int INCY = 1) const;
*/
float Epetra_BLAS_DOT_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const float * Y, const int INCX, const int INCY );

/* Original C++ prototype:
   double DOT(const int N, const double * X, const double * Y, const int INCX = 1, const int INCY = 1) const;
*/
double Epetra_BLAS_DOT_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const double * Y, const int INCX, const int INCY );

/* Original C++ prototype:
   float NRM2(const int N, const float * X, const int INCX = 1) const;
*/
float Epetra_BLAS_NRM2_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX );

/* Original C++ prototype:
   double NRM2(const int N, const double * X, const int INCX = 1) const;
*/
double Epetra_BLAS_NRM2_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX );

/* Original C++ prototype:
   void SCAL( const int N, const float ALPHA, float * X, const int INCX = 1) const;
*/
void Epetra_BLAS_SCAL_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  float * X, const int INCX );

/* Original C++ prototype:
   void SCAL( const int N, const double ALPHA, double * X, const int INCX = 1) const;
*/
void Epetra_BLAS_SCAL_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  double * X, const int INCX );

/* Original C++ prototype:
   void COPY( const int N, const float * X, float * Y, const int INCX = 1, const int INCY = 1) const;
*/
void Epetra_BLAS_COPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, float * Y, 
  const int INCX, const int INCY );

/* Original C++ prototype:
   void COPY( const int N, const double * X, double * Y, const int INCX = 1, const int INCY = 1) const;
*/
void Epetra_BLAS_COPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  double * Y, const int INCX, const int INCY );

/* Original C++ prototype:
   int IAMAX( const int N, const float * X, const int INCX = 1) const;
*/
int Epetra_BLAS_IAMAX_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float * X, 
  const int INCX );

/* Original C++ prototype:
   int IAMAX( const int N, const double * X, const int INCX = 1) const;
*/
int Epetra_BLAS_IAMAX_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double * X, 
  const int INCX );

/* Original C++ prototype:
   void AXPY( const int N, const float ALPHA, const float * X, float * Y, const int INCX = 1, const int INCY = 1) const;
*/
void Epetra_BLAS_AXPY_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const float ALPHA, 
  const float * X, float * Y, const int INCX, const int INCY );

/* Original C++ prototype:
   void AXPY( const int N, const double ALPHA, const double * X, double * Y, const int INCX = 1, const int INCY = 1) const;
*/
void Epetra_BLAS_AXPY_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const int N, const double ALPHA, 
  const double * X, double * Y, const int INCX, const int INCY );

/* Original C++ prototype:
   void GEMV(const char TRANS, const int M, const int N, const float ALPHA, const float * A, const int LDA, const float * X, const float BETA, float * Y, const int INCX = 1, const int INCY = 1) const;
*/
void Epetra_BLAS_GEMV_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const float ALPHA, const float * A, const int LDA, 
  const float * X, const float BETA, float * Y, const int INCX, 
  const int INCY );

/* Original C++ prototype:
   void GEMV(const char TRANS, const int M, const int N, const double ALPHA, const double * A, const int LDA, const double * X, const double BETA, double * Y, const int INCX = 1, const int INCY = 1) const;
*/
void Epetra_BLAS_GEMV_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANS, const int M, 
  const int N, const double ALPHA, const double * A, const int LDA, 
  const double * X, const double BETA, double * Y, const int INCX, 
  const int INCY );

/* Original C++ prototype:
   void GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K, const float ALPHA, const float * A, const int LDA, const float * B, const int LDB, const float BETA, float * C, const int LDC) const;
*/
void Epetra_BLAS_GEMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const float ALPHA, 
  const float * A, const int LDA, const float * B, const int LDB, 
  const float BETA, float * C, const int LDC );

/* Original C++ prototype:
   void GEMM(const char TRANSA, const char TRANSB, const int M, const int N, const int K, const double ALPHA, const double * A, const int LDA, const double * B, const int LDB, const double BETA, double * C, const int LDC) const;
*/
void Epetra_BLAS_GEMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char TRANSA, const char TRANSB, 
  const int M, const int N, const int K, const double ALPHA, 
  const double * A, const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC );

/* Original C++ prototype:
   void SYMM(const char SIDE, const char UPLO, const int M, const int N, const float ALPHA, const float * A, const int LDA, const float * B, const int LDB, const float BETA, float * C, const int LDC) const;
*/
void Epetra_BLAS_SYMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const float ALPHA, const float * A, 
  const int LDA, const float * B, const int LDB, const float BETA, 
  float * C, const int LDC );

/* Original C++ prototype:
   void SYMM(const char SIDE, const char UPLO, const int M, const int N, const double ALPHA, const double * A, const int LDA, const double * B, const int LDB, const double BETA, double * C, const int LDC) const;
*/
void Epetra_BLAS_SYMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const int M, const int N, const double ALPHA, const double * A, 
  const int LDA, const double * B, const int LDB, 
  const double BETA, double * C, const int LDC );

/* Original C++ prototype:
   void TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, const int N, const float ALPHA, const float * A, const int LDA, float * B, const int LDB) const;
*/
void Epetra_BLAS_TRMM_Float ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const float ALPHA, const float * A, const int LDA, float * B, 
  const int LDB );

/* Original C++ prototype:
   void TRMM(const char SIDE, const char UPLO, const char TRANSA, const char DIAG, const int M, const int N, const double ALPHA, const double * A, const int LDA, double * B, const int LDB) const;
*/
void Epetra_BLAS_TRMM_Double ( 
  CT_Epetra_BLAS_ID_t selfID, const char SIDE, const char UPLO, 
  const char TRANSA, const char DIAG, const int M, const int N, 
  const double ALPHA, const double * A, const int LDA, double * B, 
  const int LDB );


#ifdef __cplusplus
} /* extern "C" */
#endif


#endif /* CEPETRA_BLAS_H */

