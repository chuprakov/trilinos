// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_UTILITIES_DECL_HPP
#define MUELU_UTILITIES_DECL_HPP

#include <unistd.h> //necessary for "sleep" function in debugging methods

#include "MueLu_ConfigDefs.hpp"

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_Matrix_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>
#include <Xpetra_MatrixFactory_fwd.hpp>
#include <Xpetra_CrsMatrixWrap_fwd.hpp>
#include <Xpetra_CrsMatrix_fwd.hpp>
#include <Xpetra_BlockedCrsMatrix_fwd.hpp>
#include <Xpetra_Vector_fwd.hpp>
#include <Xpetra_VectorFactory_fwd.hpp>
#include <Xpetra_MultiVector_fwd.hpp>
#include <Xpetra_MultiVectorFactory_fwd.hpp>
#include <Xpetra_ImportFactory_fwd.hpp>

#ifdef HAVE_MUELU_EPETRA
namespace Xpetra {
  class EpetraCrsMatrix; // TODO: replace by include of _fwd.hpp
}

// needed because of inlined function
//TODO: remove inline function?
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#endif

#include "MueLu_Exceptions.hpp"

#ifdef HAVE_MUELU_EPETRAEXT
class Epetra_CrsMatrix;
class Epetra_MultiVector;
#endif

#ifdef HAVE_MUELU_TPETRA
#include <Xpetra_TpetraMultiVector_fwd.hpp>
#include <Xpetra_TpetraCrsMatrix_fwd.hpp>

namespace Tpetra {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>                    class MultiVector;
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps> class CrsMatrix;
}

#endif

// MPI helper
#define sumAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));
#define minAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MIN, in, Teuchos::outArg(out));
#define maxAll(rcpComm, in, out)                                        \
  Teuchos::reduceAll(*rcpComm, Teuchos::REDUCE_MAX, in, Teuchos::outArg(out));

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  //defined after Utils class
  template<typename SC,typename LO,typename GO,typename NO, typename LMO>
  RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB);
#endif

  /*!
    @class Utils
    @brief MueLu utility class.

    This class provides a number of static helper methods. Some are temporary and will eventually
    go away, while others should be moved to Xpetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = KokkosClassic::DefaultNode::DefaultNodeType,
            class LocalMatOps   = typename KokkosClassic::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Utils {
#undef MUELU_UTILITIES_SHORT
#include "MueLu_UseShortNames.hpp"

  public:
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

#ifdef HAVE_MUELU_EPETRA
    //! Helper utility to pull out the underlying Epetra objects from an Xpetra object
    // @{
    static RCP<const Epetra_MultiVector>                    MV2EpetraMV(RCP<MultiVector> const Vec);
    static RCP<      Epetra_MultiVector>                    MV2NonConstEpetraMV(RCP<MultiVector> Vec);

    static const Epetra_MultiVector&                        MV2EpetraMV(const MultiVector& Vec);
    static       Epetra_MultiVector&                        MV2NonConstEpetraMV(MultiVector& Vec);

    static RCP<const Epetra_CrsMatrix>                      Op2EpetraCrs(RCP<const Matrix> Op);
    static RCP<      Epetra_CrsMatrix>                      Op2NonConstEpetraCrs(RCP<Matrix> Op);

    static const Epetra_CrsMatrix&                          Op2EpetraCrs(const Matrix& Op);
    static       Epetra_CrsMatrix&                          Op2NonConstEpetraCrs(Matrix& Op);

    static const Epetra_Map&                                Map2EpetraMap(const Map& map);
    // @}
#endif

#ifdef HAVE_MUELU_TPETRA
    //! Helper utility to pull out the underlying Tpetra objects from an Xpetra object
    // @{
    static RCP<const Tpetra::MultiVector<SC,LO,GO,NO> >     MV2TpetraMV(RCP<MultiVector> const Vec);
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV(RCP<MultiVector> Vec);
    static RCP<      Tpetra::MultiVector<SC,LO,GO,NO> >     MV2NonConstTpetraMV2(MultiVector& Vec);

    static const Tpetra::MultiVector<SC,LO,GO,NO>&          MV2TpetraMV(const MultiVector& Vec);
    static       Tpetra::MultiVector<SC,LO,GO,NO>&          MV2NonConstTpetraMV(MultiVector& Vec);

    static RCP<const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >   Op2TpetraCrs(RCP<const Matrix> Op);
    static RCP<      Tpetra::CrsMatrix<SC,LO,GO,NO,LMO> >   Op2NonConstTpetraCrs(RCP<Matrix> Op);

    static const Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>&        Op2TpetraCrs(const Matrix& Op);
    static       Tpetra::CrsMatrix<SC,LO,GO,NO,LMO>&        Op2NonConstTpetraCrs(Matrix& Op);

    static const RCP<const Tpetra::Map<LO, GO, NO> >        Map2TpetraMap(const Map& map);
#endif

    /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    */
    static RCP<Matrix> Multiply(const Matrix & A,
                                bool transposeA,
                                const Matrix & B,
                                bool transposeB,
                                //Teuchos::FancyOStream &fos = *(Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout))),
                                Teuchos::FancyOStream &fos,
                                bool callFillCompleteOnResult = true,
                                bool doOptimizeStorage        = true,
                                bool allowMLMultiply          = true) {
      return Utils<SC,LO,GO,NO,LMO>::Multiply(A, transposeA, B, transposeB, Teuchos::null, fos, callFillCompleteOnResult, doOptimizeStorage, allowMLMultiply);
    }


    /*! @brief Helper function to do matrix-matrix multiply

    Returns C = AB.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B

    @param C_in advanced usage. Use Teuchos::null by default.
           When C_in is available, its memory is reused to build
           This is useful in the case of multiple solve: if the pattern of C does not change, we can keep the memory and pattern of previous C matrix (C_in)
           C_in is modified in place and is not valid after the call.
    @param callFillCompleteOnResult if true, the resulting matrix should be fillComplete'd
    @param doOptimizedStorage if true, optimize storage
    @param allowMLMultiply    if true, allow usage of ML's matrix matrix multiply

           ML MxM multiply does not reuse the pattern of C_in. If a C_in matrix is provided, then it is ignored.
           This can create a memory penalty if both the useless C_in and the new C are present in memory at the same time
           => Do not enable ML MxM at the same time as the option "reuse pattern".
    */
    static RCP<Matrix> Multiply(const Matrix& A,
                                bool transposeA,
                                const Matrix& B,
                                bool transposeB,
                                RCP<Matrix> C_in,
                                //Teuchos::FancyOStream &fos = *(Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout)))
                                Teuchos::FancyOStream &fos,
                                bool callFillCompleteOnResult = true,
                                bool doOptimizeStorage        = true,
                                bool allowMLMultiply          = true
                                );

#ifdef HAVE_MUELU_EPETRAEXT
    // Michael Gee's MLMultiply
    static RCP<Epetra_CrsMatrix> MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
                                                     const Epetra_CrsMatrix& epB,
                                                     Teuchos::FancyOStream& fos);
#endif //ifdef HAVE_MUELU_EPETRAEXT

    /*! @brief Helper function to do matrix-matrix multiply "in-place"

    Returns RCP to non-constant Xpetra::BlockedCrsMatrix.

    @param A left matrix
    @param transposeA if true, use the transpose of A
    @param B right matrix
    @param transposeB if true, use the transpose of B
    @param doOptimizeStorage if true, the resulting matrix should be fillComplete'd
    */
    static RCP<BlockedCrsMatrix> TwoMatrixMultiplyBlock(BlockedCrsMatrix& A, bool transposeA,
                                                        BlockedCrsMatrix& B, bool transposeB,
                                                        bool doFillComplete    = true,
                                                        bool doOptimizeStorage = true);

    /*! @brief Extract Matrix Diagonal

    Returns Matrix diagonal in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetMatrixDiagonal(const Matrix& A);

    /*! @brief Extract Matrix Diagonal of lumped matrix

    Returns Matrix diagonal of lumped matrix in ArrayRCP.

    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static Teuchos::ArrayRCP<SC> GetLumpedMatrixDiagonal(const Matrix& A);

    /*! @brief Extract Overlapped Matrix Diagonal

    Returns overlapped Matrix diagonal in ArrayRCP.

    The local overlapped diagonal has an entry for each index in A's column map.
    NOTE -- it's assumed that A has been fillComplete'd.
    */
    static RCP<Vector> GetMatrixOverlappedDiagonal(const Matrix& A);

    /*! @brief Left scale matrix by an arbitrary vector.

    Algorithmically, this left scales a matrix by a diagonal matrix.
    The inverse of a diagonal matrix can also be applied.

    @param Op matrix to be scaled
    @param scalingVector vector that represents diagonal matrix
    @doInverse Indicates whether the inverse of the diagonal matrix should be applied.  (Default is to use inverse.)
    */
    static void ScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doInverse = true);


    // TODO: should NOT return an Array. Definition must be changed to:
    // - ArrayRCP<> ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS)
    // or
    // - void ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS, Array &)
    static Teuchos::Array<Magnitude> ResidualNorm(const Matrix& Op, const MultiVector& X, const MultiVector& RHS);

    static RCP<MultiVector> Residual(const Matrix& Op, const MultiVector& X, const MultiVector& RHS);

    // NOTE:
    // A better place for the Read/Write function is probably Xpetra

    //! Read/Write methods
    //@{
    /*! @brief Save map to file. */
    static void Write(const std::string& fileName, const Map& M);

    /*! @brief Save vector to file in Matrix Market format.  */
    static void Write(const std::string& fileName, const MultiVector& Vec);

    /*! @brief Save matrix to file in Matrix Market format. */
    static void Write(const std::string& fileName, const Matrix& Op);

    //! @brief Read matrix from file in Matrix Market or binary format.
    static Teuchos::RCP<Matrix> Read(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm, bool binary = false);

    /*! @brief Read matrix from file in Matrix Market or binary format.

        If only rowMap is specified, then it is used for the domainMap and rangeMap, as well.
    */
    static Teuchos::RCP<Matrix> Read(const std::string&   filename,
                                     const RCP<const Map> rowMap,
                                           RCP<const Map> colMap           = Teuchos::null,
                                     const RCP<const Map> domainMap        = Teuchos::null,
                                     const RCP<const Map> rangeMap         = Teuchos::null,
                                     const bool           callFillComplete = true,
                                     const bool           tolerant         = false,
                                     const bool           debug            = false);
    //@}

    static void PauseForDebugger();

    /*! @brief Simple transpose for Tpetra::CrsMatrix types

        Note:  This is very inefficient, as it inserts one entry at a time.
    */

    /*! @brief Power method.

    @param A matrix
    @param scaleByDiag if true, estimate the largest eigenvalue of \f$ D^; A \f$.
    @param niters maximum number of iterations
    @param tolerance stopping tolerance
    @verbose if true, print iteration information

    (Shamelessly grabbed from tpetra/examples.)
    */
    static Scalar PowerMethod(const Matrix& A, bool scaleByDiag = true,
                              LO niters = 10, Magnitude tolerance = 1e-2, bool verbose = false, unsigned int seed = 123);

    static void MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse = true,
                                 bool doFillComplete = true, bool doOptimizeStorage = true);

    static void MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                                        bool doFillComplete, bool doOptimizeStorage);

    static RCP<Teuchos::FancyOStream> MakeFancy(std::ostream& os);

    /*! @brief Squared distance between two rows in a multivector

       Used for coordinate vectors.
    */
    static typename Teuchos::ScalarTraits<Scalar>::magnitudeType Distance2(const MultiVector& v, LocalOrdinal i0, LocalOrdinal i1);

    /*! @brief Detect Dirichlet rows

        @param[in] A matrix
        @param[in] tol If a row entry's magnitude is less than or equal to this tolerance, the entry is treated as zero.

        @return boolean array.  The ith entry is true iff row i is a Dirichlet row.
    */
    static Teuchos::ArrayRCP<const bool> DetectDirichletRows(const Matrix& A, const Magnitude& tol = Teuchos::ScalarTraits<SC>::zero());

    /*! @brief print matrix info
    */
    static std::string PrintMatrixInfo(const Matrix& A, const std::string& msgTag, RCP<const Teuchos::ParameterList> params = Teuchos::null);

  }; // class Utils

#ifdef HAVE_MUELU_EPETRA
  //This non-member templated function exists so that the matrix-matrix multiply will compile if Epetra, Tpetra, and ML are enabled.
  template<class SC,class LO,class GO,class NO, class LMO>
  RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB) {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
    return Teuchos::null;
  }

  typedef KokkosClassic::DefaultNode::DefaultNodeType KDNT;
  typedef KokkosClassic::DefaultKernels<void,int,KokkosClassic::DefaultNode::DefaultNodeType>::SparseOps KDKSO;

  //specialization for the case of ScalarType=double and LocalOrdinal=GlobalOrdinal=int
  template<>
  inline RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> &epAB) {
    RCP<Xpetra::EpetraCrsMatrix> tmpC1 = rcp(new Xpetra::EpetraCrsMatrix(epAB));
    RCP<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> > tmpC2 = rcp_implicit_cast<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> >(tmpC1);
    RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > tmpC3 = rcp(new Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO>(tmpC2));
    return tmpC3;
  }
#endif

  //! Little helper function to convert non-string types to strings
  template<class T>
  std::string toString(const T& what) {
    std::ostringstream buf;
    buf << what;
    return buf.str();
  }

  /*!
    @class Utils2
    @brief MueLu utility class.

    Separate class for utilities that need a specialization for Epetra.
  */
  template <class Scalar,
            class LocalOrdinal  = int,
            class GlobalOrdinal = LocalOrdinal,
            class Node          = KokkosClassic::DefaultNode::DefaultNodeType,
            class LocalMatOps   = typename KokkosClassic::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps > //TODO: or BlockSparseOp ?
  class Utils2 {

#include "MueLu_UseShortNames.hpp"

  public:

    /*! @brief Transpose a Xpetra::Matrix

    Note: Currently, an error is thrown if the matrix isn't a Tpetra::CrsMatrix or Epetra_CrsMatrix.
    In principle, however, we could allow any Epetra_RowMatrix because the Epetra transposer does.
    */
    static RCP<Matrix> Transpose(Matrix& Op, bool optimizeTranspose = false);

    //! Scale an Epetra matrix.
    static void MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doFillComplete, bool doOptimizeStorage);

    /*! @brief Helper function to calculate B = alpha*A + beta*B.

    @param A      left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha  scalar multiplier for A
    @param B      right matrix operand
    @param beta   scalar multiplier for B

    @return sum in B.

    Note that B does not have to be fill-completed.
    */
    static void TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta);

    /*! @brief Helper function to calculate C = alpha*A + beta*B.

    @param A          left matrix operand
    @param transposeA indicate whether to use transpose of A
    @param alpha      scalar multiplier for A, defaults to 1.0
    @param B          right matrix operand
    @param transposeB indicate whether to use transpose of B
    @param beta       scalar multiplier for B, defaults to 1.0
    @param C          resulting sum

    It is up to the caller to ensure that the resulting matrix sum is fillComplete'd.
    */
    static void TwoMatrixAdd(const Matrix& A, bool transposeA, const SC& alpha,
                             const Matrix& B, bool transposeB, const SC& beta,
                             RCP<Matrix>& C,  Teuchos::FancyOStream &fos, bool AHasFixedNnzPerRow = false);

    static RCP<MultiVector> ReadMultiVector (const std::string& fileName, const RCP<const Map>& map);
    static RCP<const Map>   ReadMap         (const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm);
  }; // class Utils2

  // specialization Utils2 for SC=double, LO=GO=int
  template<>
  class Utils2<double,int,int> {
    typedef double                                                     SC;
    typedef int                                                        LO;
    typedef int                                                        GO;
    typedef KokkosClassic::DefaultNode::DefaultNodeType                NO;
    typedef KokkosClassic::DefaultKernels<double,int,NO>::SparseOps    LMO;
    typedef Xpetra::Map<int,int,NO>                                    Map;
    typedef Xpetra::Matrix<double,int,int,NO,LMO>                      Matrix;
    typedef Xpetra::MultiVector<double,int,int,NO>                     MultiVector;

  public:

    static RCP<Matrix>      Transpose               (Matrix& Op, bool optimizeTranspose = false);
    static void             MyOldScaleMatrix_Epetra (Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doFillComplete, bool doOptimizeStorage);
    static void             TwoMatrixAdd            (const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta);
    static void             TwoMatrixAdd            (const Matrix& A, bool transposeA, SC alpha,
                                                     const Matrix& B, bool transposeB, SC beta,
                                                     RCP<Matrix>& C,  Teuchos::FancyOStream & fos, bool AHasFixedNnzPerRow = false);
    static RCP<MultiVector> ReadMultiVector         (const std::string& fileName, const RCP<const Map>& map);
    static RCP<const Map>   ReadMap                 (const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm);
  };


} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DECL_HPP
