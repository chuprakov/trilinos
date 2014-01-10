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
#ifndef MUELU_UTILITIES_DEF_HPP
#define MUELU_UTILITIES_DEF_HPP

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ParameterList.hpp>

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_EPETRA
# ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
# endif
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_MultiVectorIn.h>
#include <EpetraExt_BlockMapIn.h>
#include <Xpetra_EpetraUtils.hpp>
#include <Xpetra_EpetraMultiVector.hpp>
#include <EpetraExt_BlockMapOut.h>
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <TpetraExt_MatrixMatrix.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif // HAVE_MUELU_TPETRA

#ifdef HAVE_MUELU_EPETRA
#include <Xpetra_EpetraMap.hpp>
#endif //ifdef HAVE_MUELU_EPETRA

#include <Xpetra_Map.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include "Xpetra_DefaultPlatform.hpp"

#include <XpetraExt_MatrixMatrix.hpp>  // standard matrix matrix routines

#include <MueLu_Utilities_decl.hpp>

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_ML)
#include <ml_operator.h>
#include <ml_epetra_utils.h>
#endif

namespace MueLu {

#ifdef HAVE_MUELU_EPETRA
  using Xpetra::EpetraCrsMatrix;   // TODO: mv in Xpetra_UseShortNamesScalar
  using Xpetra::EpetraMultiVector;
#endif

#ifdef HAVE_MUELU_EPETRA
  //defined after Utils class
  template<typename SC,typename LO,typename GO,typename NO, typename LMO>
  RCP<Xpetra::CrsMatrixWrap<SC,LO,GO,NO,LMO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB);
#endif

#ifdef HAVE_MUELU_EPETRA
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2EpetraMV(const RCP<MultiVector> Vec) {
    RCP<const EpetraMultiVector > tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstEpetraMV(RCP<MultiVector> Vec) {
    RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed");
    return tmpVec->getEpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Epetra_MultiVector& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstEpetraMV(MultiVector &Vec) {
    const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(Vec);
    return *(tmpVec.getEpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Epetra_MultiVector& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2EpetraMV(const MultiVector& Vec) {
    const EpetraMultiVector& tmpVec = dynamic_cast<const EpetraMultiVector&>(Vec);
    return *(tmpVec.getEpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2EpetraCrs(RCP<const Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrix();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(RCP<Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
    return tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Epetra_CrsMatrix& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2EpetraCrs(const Matrix& Op) {
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<const EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrix();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Epetra_CrsMatrix& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(Matrix& Op) {
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        EpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<EpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getEpetra_CrsMatrixNonConst();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Epetra_Map& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Map2EpetraMap(const Map& map) {
    RCP<const Xpetra::EpetraMap> xeMap = rcp_dynamic_cast<const Xpetra::EpetraMap>(rcpFromRef(map));
    if (xeMap == Teuchos::null)
      throw Exceptions::BadCast("Utils::Map2EpetraMap : Cast from Xpetra::Map to Xpetra::EpetraMap failed");
    return xeMap->getEpetra_Map();
  }
#endif

#ifdef HAVE_MUELU_TPETRA
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2TpetraMV(RCP<MultiVector> const Vec) {
    RCP<const TpetraMultiVector > tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV(RCP<MultiVector> Vec) {
    RCP<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed");
    return tmpVec->getTpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV(MultiVector& Vec) {
    const TpetraMultiVector& tmpVec = dynamic_cast<const TpetraMultiVector&>(Vec);
    return *(tmpVec.getTpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV2(MultiVector &Vec) {
    const TpetraMultiVector& tmpVec = dynamic_cast<const TpetraMultiVector&>(Vec);
    return tmpVec.getTpetra_MultiVector();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2TpetraMV(const MultiVector& Vec) {
    const TpetraMultiVector& tmpVec = dynamic_cast<const TpetraMultiVector&>(Vec);
    return *(tmpVec.getTpetra_MultiVector());
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(RCP<const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Op) {
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrix();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstTpetraCrs(RCP<Matrix> Op) {
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(crsOp->getCrsMatrix());
    if (tmp_ECrsMtx == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
    return tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(const Matrix& Op) {
    try {
      const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
      try {
        const TpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<const TpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrix();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstTpetraCrs(Matrix& Op) {
    try {
      CrsMatrixWrap& crsOp = dynamic_cast<CrsMatrixWrap&>(Op);
      try {
        TpetraCrsMatrix& tmp_ECrsMtx = dynamic_cast<TpetraCrsMatrix&>(*crsOp.getCrsMatrix());
        return *tmp_ECrsMtx.getTpetra_CrsMatrixNonConst();
      } catch (std::bad_cast) {
        throw Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed");
      }
    } catch (std::bad_cast) {
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  const RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal,Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Map2TpetraMap(const Map& map) {
    const RCP<const TpetraMap>& tmp_TMap = rcp_dynamic_cast<const TpetraMap>(rcpFromRef(map));
    if (tmp_TMap == Teuchos::null)
      throw Exceptions::BadCast("Utils::Map2TpetraMap : Cast from Xpetra::Map to Xpetra::TpetraMap failed");
    return tmp_TMap->getTpetra_Map();
  }
#endif






  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Jacobi(Scalar omega,
									const Vector& Dinv,
									const Matrix& A,
									const Matrix& B,
									RCP<Matrix> C_in,
									Teuchos::FancyOStream &fos) {
    // Sanity checks
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    // Default case: Xpetra Jacobi
    RCP<Matrix> C = C_in;
    if (C == Teuchos::null)
      C = MatrixFactory::Build(B.getRowMap(),Teuchos::OrdinalTraits<LO>::zero());

    Xpetra::MatrixMatrix::Jacobi(omega, Dinv, A, B, *C, true,true);
    C->CreateView("stridedMaps", rcpFromRef(A),false, rcpFromRef(B), false);
    return C;
  } //Jacobi


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Multiply(const Matrix& A, bool transposeA,
                                                                          const Matrix& B, bool transposeB,
                                                                          RCP<Matrix> C_in,
                                                                          Teuchos::FancyOStream &fos,
                                                                          bool doFillComplete,
                                                                          bool doOptimizeStorage) {


    // Preconditions
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    // Optimization using ML Multiply when available
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_ML_MMM)
    if (B.getDomainMap()->lib() == Xpetra::UseEpetra && !transposeA && !transposeB) {
      RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(rcpFromRef(A));
      RCP<const Epetra_CrsMatrix> epB = Op2EpetraCrs(rcpFromRef(B));
      RCP<Epetra_CrsMatrix>       epC = MLTwoMatrixMultiply(*epA, *epB, fos);

      RCP<Matrix> C = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(epC);
      if (doFillComplete) {
        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage", doOptimizeStorage);
        C->fillComplete(B.getDomainMap(), A.getRangeMap(), params);
      }

      C->CreateView("stridedMaps", rcpFromRef(A), transposeA, rcpFromRef(B), transposeB);
      return C;
    }
#endif // EPETRA + EPETRAEXT + ML

    // Default case: Xpetra Multiply
    RCP<Matrix> C = C_in;

    if (C == Teuchos::null) {
      double nnzPerRow = Teuchos::as<double>(0);

      if (A.getDomainMap()->lib() == Xpetra::UseTpetra) {
        // For now, follow what ML and Epetra do.
        GO numRowsA = A.getGlobalNumRows();
        GO numRowsB = B.getGlobalNumRows();
        nnzPerRow = sqrt(Teuchos::as<double>(A.getGlobalNumEntries())/numRowsA) +
                    sqrt(Teuchos::as<double>(B.getGlobalNumEntries())/numRowsB) - 1;
        nnzPerRow *=  nnzPerRow;
        double totalNnz = nnzPerRow * A.getGlobalNumRows() * 0.75 + 100;
        double minNnz = Teuchos::as<double>(1.2 * A.getGlobalNumEntries());
        if (totalNnz < minNnz)
          totalNnz = minNnz;
        nnzPerRow = totalNnz / A.getGlobalNumRows();

        fos << "Utils::Multiply : Estimate for nnz per row of product matrix = " << Teuchos::as<LO>(nnzPerRow) << std::endl;
      }

      if (transposeA) C = MatrixFactory::Build(A.getDomainMap(), Teuchos::as<LO>(nnzPerRow));
      else            C = MatrixFactory::Build(A.getRowMap(),    Teuchos::as<LO>(nnzPerRow));

    } else {
      C->resumeFill(); // why this is not done inside of Tpetra MxM?
      fos << "Reuse C pattern" << std::endl;
    }

    Xpetra::MatrixMatrix::Multiply(A, transposeA, B, transposeB, *C, doFillComplete, doOptimizeStorage);

    // Fill strided maps information
    // This is necessary since the ML matrix matrix multiplication routine has no handling for this
    // TODO: move this call to MLMultiply...
    C->CreateView("stridedMaps", rcpFromRef(A), transposeA, rcpFromRef(B), transposeB);

    return C;
  } //Multiply()

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA, const Epetra_CrsMatrix& epB, Teuchos::FancyOStream &fos) {
#if defined(HAVE_MUELU_ML_MMM)
    ML_Comm* comm;
    ML_Comm_Create(&comm);
    fos << "****** USING ML's MATRIX MATRIX MULTIPLY (LNM version) ******" << std::endl;
#ifdef HAVE_MPI
    // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
    const Epetra_MpiComm * Mcomm = dynamic_cast<const Epetra_MpiComm*>(&(epA.Comm()));
    if (Mcomm)
      ML_Comm_Set_UsrComm(comm,Mcomm->GetMpiComm());
# endif
    //in order to use ML, there must be no indices missing from the matrix column maps.
    EpetraExt::CrsMatrix_SolverMap Atransform;
    EpetraExt::CrsMatrix_SolverMap Btransform;
    const Epetra_CrsMatrix& A = Atransform(const_cast<Epetra_CrsMatrix&>(epA));
    const Epetra_CrsMatrix& B = Btransform(const_cast<Epetra_CrsMatrix&>(epB));

    if (!A.Filled())    throw Exceptions::RuntimeError("A has to be FillCompleted");
    if (!B.Filled())    throw Exceptions::RuntimeError("B has to be FillCompleted");

    // create ML operators from EpetraCrsMatrix
    ML_Operator* ml_As = ML_Operator_Create(comm);
    ML_Operator* ml_Bs = ML_Operator_Create(comm);
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&A),ml_As); // Should we test if the lightweight wrapper is actually used or if WrapEpetraCrsMatrix fall backs to the heavy one?
    ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&B),ml_Bs);
    ML_Operator* ml_AtimesB = ML_Operator_Create(comm);
    {
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("ML_2matmult kernel"));
      ML_2matmult(ml_As,ml_Bs,ml_AtimesB,ML_CSR_MATRIX); // do NOT use ML_EpetraCRS_MATRIX!!!
    }
    ML_Operator_Destroy(&ml_As);
    ML_Operator_Destroy(&ml_Bs);

    // For ml_AtimesB we have to reconstruct the column map in global indexing,
    // The following is going down to the salt-mines of ML ...
    // note: we use integers, since ML only works for Epetra...
    int N_local = ml_AtimesB->invec_leng;
    ML_CommInfoOP* getrow_comm = ml_AtimesB->getrow->pre_comm;
    if (!getrow_comm)   throw(Exceptions::RuntimeError("ML_Operator does not have a CommInfo"));
    ML_Comm* comm_AB = ml_AtimesB->comm;   // get comm object
    if (N_local != B.DomainMap().NumMyElements())
      throw(Exceptions::RuntimeError("Mismatch in local row dimension between ML and Epetra"));
    int N_rcvd = 0;
    int N_send = 0;
    int flag   = 0;
    for (int i = 0; i < getrow_comm->N_neighbors; i++) {
      N_rcvd += (getrow_comm->neighbors)[i].N_rcv;
      N_send += (getrow_comm->neighbors)[i].N_send;
      if ( ((getrow_comm->neighbors)[i].N_rcv    != 0) &&
           ((getrow_comm->neighbors)[i].rcv_list != NULL) ) flag = 1;
    }
    // For some unknown reason, ML likes to have stuff one larger than
    // neccessary...
    std::vector<double> dtemp(N_local + N_rcvd + 1); // "double" vector for comm function
    std::vector<int>    cmap (N_local + N_rcvd + 1); // vector for gids
    for (int i = 0; i < N_local; ++i) {
      cmap[i]  = B.DomainMap().GID(i);
      dtemp[i] = (double) cmap[i];
    }
    ML_cheap_exchange_bdry(&dtemp[0],getrow_comm,N_local,N_send,comm_AB); // do communication
    if (flag) { // process received data
      int count = N_local;
      const int neighbors = getrow_comm->N_neighbors;
      for (int i = 0; i < neighbors; i++) {
        const int nrcv = getrow_comm->neighbors[i].N_rcv;
        for (int j = 0; j < nrcv; j++)
          cmap[getrow_comm->neighbors[i].rcv_list[j]] = (int) dtemp[count++];
      }
    } else {
      for (int i = 0; i < N_local+N_rcvd; ++i)
        cmap[i] = (int)dtemp[i];
    }
    dtemp.clear();  // free double array

    // we can now determine a matching column map for the result
    Epetra_Map gcmap(-1, N_local+N_rcvd, &cmap[0], B.ColMap().IndexBase(), A.Comm());

    int     allocated = 0;
    int     rowlength;
    double* val       = NULL;
    int*    bindx     = NULL;

    const int myrowlength    = A.RowMap().NumMyElements();
    const Epetra_Map& rowmap = A.RowMap();

    // Determine the maximum bandwith for the result matrix.
    // replaces the old, very(!) memory-consuming guess:
    // int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
    int educatedguess = 0;
    for (int i = 0; i < myrowlength; ++i) {
      // get local row
      ML_get_matrix_row(ml_AtimesB, 1, &i, &allocated, &bindx, &val, &rowlength, 0);
      if (rowlength>educatedguess)
        educatedguess = rowlength;
    }

    // allocate our result matrix and fill it
    RCP<Epetra_CrsMatrix> result = rcp(new Epetra_CrsMatrix(::Copy, A.RangeMap(), gcmap, educatedguess, false));

    std::vector<int> gcid(educatedguess);
    for (int i = 0; i < myrowlength; ++i) {
      const int grid = rowmap.GID(i);
      // get local row
      ML_get_matrix_row(ml_AtimesB, 1, &i, &allocated, &bindx, &val, &rowlength, 0);
      if (!rowlength) continue;
      if ((int)gcid.size() < rowlength) gcid.resize(rowlength);
      for (int j = 0; j < rowlength; ++j) {
        gcid[j] = gcmap.GID(bindx[j]);
        if (gcid[j] < 0)
          throw Exceptions::RuntimeError("Error: cannot find gcid!");
      }
      int err = result->InsertGlobalValues(grid, rowlength, val, &gcid[0]);
      if (err != 0 && err != 1) {
        std::ostringstream errStr;
        errStr << "Epetra_CrsMatrix::InsertGlobalValues returned err=" << err;
        throw Exceptions::RuntimeError(errStr.str());
      }
    }
    // free memory
    if (bindx) ML_free(bindx);
    if (val)   ML_free(val);
    ML_Operator_Destroy(&ml_AtimesB);
    ML_Comm_Destroy(&comm);

    return result;
#else // no MUELU_ML
    TEUCHOS_TEST_FOR_EXCEPTION(true, Xpetra::Exceptions::RuntimeError,
                               "HAVE_MUELU_ML compiler flag not set. no ML multiply available.");
    return Teuchos::null;
#endif
  }
#endif //ifdef HAVE_MUELU_EPETRAEXT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::
  TwoMatrixMultiplyBlock(BlockedCrsMatrix& A, bool transposeA,
                         BlockedCrsMatrix& B, bool transposeB,
                         bool doFillComplete,
                         bool doOptimizeStorage) {
    if (transposeA || transposeB)
      throw Exceptions::RuntimeError("TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true");

    // Preconditions
    if (!A.isFillComplete())
      throw Exceptions::RuntimeError("A is not fill-completed");
    if (!B.isFillComplete())
      throw Exceptions::RuntimeError("B is not fill-completed");

    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgmapextractor = A.getRangeMapExtractor();
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domapextractor = B.getDomainMapExtractor();

    RCP<BlockedCrsMatrix> C = rcp(new BlockedCrsMatrix(rgmapextractor, domapextractor, 33 /* TODO fix me */));

    for (size_t i = 0; i < A.Rows(); ++i) { // loop over all block rows of A
      for (size_t j = 0; j < B.Cols(); ++j) { // loop over all block columns of B
        // empty CrsMatrixWrap
        RCP<Matrix> Cij = MatrixFactory::Build(A.getRangeMap(i), 33 /* TODO fix me */);

        for (size_t l = 0; l < B.Rows(); ++l) { // loop for calculating entry C_{ij}
          RCP<CrsMatrix> crmat1 = A.getMatrix(i,l);
          RCP<CrsMatrix> crmat2 = B.getMatrix(l,j);
          RCP<CrsMatrixWrap> crop1 = rcp(new CrsMatrixWrap(crmat1));
          RCP<CrsMatrixWrap> crop2 = rcp(new CrsMatrixWrap(crmat2));

          RCP<Matrix> temp = MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Multiply(*crop1, false, *crop2, false,
          *(Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout))));

          // sum up
          MueLu::Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(*temp, false, 1.0, *Cij, 1.0);
        }

        Cij->fillComplete(B.getDomainMap(j), A.getRangeMap(i));

        RCP<CrsMatrixWrap> crsCij = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(Cij);
        TEUCHOS_TEST_FOR_EXCEPTION( Cij==Teuchos::null, Xpetra::Exceptions::BadCast,
                                    "MatrixFactory failed in generating a CrsMatrixWrap." );

        RCP<CrsMatrix> crsMatCij = crsCij->getCrsMatrix();
        C->setMatrix(i, j, crsMatCij);
      }
    }

    if (doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

    return C;
  } // TwoMatrixMultiplyBlock

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<Scalar> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMatrixDiagonal(const Matrix& A) {
    size_t numRows = A.getRowMap()->getNodeNumElements();
    Teuchos::ArrayRCP<SC> diag(numRows);

    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i = 0; i < numRows; ++i) {
      A.getLocalRowView(i, cols, vals);

      LO j = 0;
      for (; j < cols.size(); ++j) {
        if (Teuchos::as<size_t>(cols[j]) == i) {
          diag[i] = vals[j];
          break;
        }
      }
      if (j == cols.size()) {
        // Diagonal entry is absent
        diag[i] = Teuchos::ScalarTraits<SC>::zero();
      }
    }

    return diag;
  } //GetMatrixDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMatrixDiagonalInverse(const Matrix& A,Magnitude tol) {
    RCP<const Map> rowMap = A.getRowMap();
    RCP<Vector> diag      = VectorFactory::Build(rowMap);
    ArrayRCP<SC> diagVals = diag->getDataNonConst(0);

    size_t numRows = rowMap->getNodeNumElements();

    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i = 0; i < numRows; ++i) {
      A.getLocalRowView(i, cols, vals);

      LO j = 0;
      for (; j < cols.size(); ++j) {
        if (Teuchos::as<size_t>(cols[j]) == i) {
	  if(Teuchos::ScalarTraits<SC>::magnitude(vals[j]) > tol)
	    diagVals[i] = Teuchos::ScalarTraits<SC>::one() / vals[j];
	  else
	    diagVals[i]=Teuchos::ScalarTraits<SC>::zero();
	  break;
	}
      }
      if (j == cols.size()) {
        // Diagonal entry is absent
        diagVals[i]=Teuchos::ScalarTraits<SC>::zero();
      }
    }
    diagVals=null;

    return diag;
  } //GetMatrixDiagonalInverse



  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<Scalar> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetLumpedMatrixDiagonal(const Matrix &A) {
    size_t numRows = A.getRowMap()->getNodeNumElements();
    Teuchos::ArrayRCP<SC> diag(numRows);

    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i = 0; i < numRows; ++i) {
      A.getLocalRowView(i, cols, vals);

      diag[i] = Teuchos::ScalarTraits<Scalar>::zero();
      for (LO j = 0; j < cols.size(); ++j) {
        diag[i] += Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
      }
    }

    return diag;
  } //GetMatrixDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMatrixOverlappedDiagonal(const Matrix& A) {
    RCP<const Map> rowMap = A.getRowMap(), colMap = A.getColMap();
    RCP<Vector>    localDiag     = VectorFactory::Build(rowMap);
    ArrayRCP<SC>   localDiagVals = localDiag->getDataNonConst(0);

    Teuchos::ArrayRCP<SC> diagVals = GetMatrixDiagonal(A);
    for (LO i = 0; i < localDiagVals.size(); i++)
      localDiagVals[i] = diagVals[i];
    localDiagVals = diagVals = null;

    // TODO there's a problem with the importer from the underlying Tpetra::CrsGraph
    // TODO so right now construct an importer.
    // diagonal->doImport(*localDiag, *(A.getCrsGraph()->getImporter()), Xpetra::INSERT);
    RCP<const Import> importer = ImportFactory::Build(rowMap, colMap);

    RCP<Vector> diagonal = VectorFactory::Build(colMap);
    diagonal->doImport(*localDiag, *importer, Xpetra::INSERT);

    return diagonal;
  } //GetMatrixOverlappedDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector, bool doInverse) {
#ifdef HAVE_MUELU_TPETRA
    try {
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& tpOp = Op2NonConstTpetraCrs(Op);

      Tpetra::Vector<SC,LO,GO,NO> x(tpOp.getRowMap(), scalingVector());
      if (doInverse){
        Tpetra::Vector<SC,LO,GO,NO> xi(tpOp.getRowMap());
        xi.reciprocal(x);
        tpOp.leftScale(xi);

      } else {
        tpOp.leftScale(x);
      }
    } catch(...) {
      throw Exceptions::RuntimeError("Matrix scaling has not been implemented Epetra");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling has not been implemented Epetra");
#endif // HAVE_MUELU_TPETRA
  } //ScaleMatrix()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ResidualNorm(const Matrix& Op, const MultiVector& X, const MultiVector& RHS) {
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
    const size_t numVecs = X.getNumVectors();

    RCP<MultiVector> RES = Residual(Op, X, RHS);
    Teuchos::Array<Magnitude> norms(numVecs);
    RES->norm2(norms);

    return norms;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Residual(const Matrix& Op, const MultiVector& X, const MultiVector& RHS) {
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
    const size_t numVecs = X.getNumVectors();

    SC one = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();

    RCP<MultiVector> RES = MultiVectorFactory::Build(Op.getRangeMap(), numVecs, false); // no need to initialize to zero
    Op.apply(X, *RES, Teuchos::NO_TRANS, one, zero);
    RES->update(one, RHS, negone);

    return RES;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write(const std::string& fileName, const Matrix& Op) {
    std::string mapfile = "rowmap_" + fileName;
    Write(mapfile, *(Op.getRowMap()));

    const CrsMatrixWrap& crsOp = dynamic_cast<const CrsMatrixWrap&>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp.getCrsMatrix();
#if defined(HAVE_MUELU_EPETRA)
    const RCP<const EpetraCrsMatrix>& tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx != Teuchos::null) {
#if defined(HAVE_MUELU_EPETRAEXT)
      RCP<const Epetra_CrsMatrix> A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      int rv = EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(), *A);
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::RowMatrixToMatrixMarketFile return value of " + toString(rv));
#else
      throw Exceptions::RuntimeError("Compiled without EpetraExt");
#endif
      return;
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraCrsMatrix>& tmp_TCrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_TCrsMtx != Teuchos::null) {
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >::writeSparseFile(fileName, A);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing");
  } //Write

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ReadMultiVector(const std::string& fileName, const RCP<const Map>& map){
    Xpetra::UnderlyingLib lib = map->lib();

    if (lib == Xpetra::UseEpetra) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");

    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>                          reader_type;
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>                            map_type;
      typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>            multivector_type;

      RCP<const map_type>   temp = toTpetra(map);
      RCP<multivector_type> TMV  = reader_type::readDenseFile(fileName,map->getComm(),map->getNode(),temp);
      RCP<MultiVector>      rmv  = Xpetra::toXpetra(TMV);
      return rmv;
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    return Teuchos::null;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >
  Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ReadMap(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm) {
    if (lib == Xpetra::UseEpetra) {
      TEUCHOS_TEST_FOR_EXCEPTION(true, ::Xpetra::Exceptions::BadCast, "Epetra can only be used with Scalar=double and Ordinal=int");
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>                          reader_type;

      RCP<Node> node = rcp(new Node());

      RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > tMap = reader_type::readMapFile(fileName, comm, node);
      if (tMap.is_null())
        throw Exceptions::RuntimeError("The Tpetra::Map returned from readSparseFile() is null.");

      return Xpetra::toXpetra(tMap);
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    return Teuchos::null;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Read(const std::string& fileName, Xpetra::UnderlyingLib lib, const RCP<const Teuchos::Comm<int> >& comm, bool binary) {

    if (binary == false) {
      // Matrix Market file format (ASCII)
      if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
        Epetra_CrsMatrix *eA;
        const RCP<const Epetra_Comm> epcomm = Xpetra::toEpetra(comm);
        int rv = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), *epcomm, eA);
        if (rv != 0)
          throw Exceptions::RuntimeError("EpetraExt::MatrixMarketFileToCrsMatrix return value of " + toString(rv));

        RCP<Epetra_CrsMatrix> tmpA = rcp(eA);
        RCP<Matrix>           A    = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(tmpA);

        return A;
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
      } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
//        typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Kokkos::SerialNode, typename KokkosClassic::DefaultKernels<void, LocalOrdinal, Kokkos::SerialNode>::SparseOps> sparse_matrix_type;
	typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> sparse_matrix_type;

        typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

        //RCP<Node> node = Xpetra::DefaultPlatform::getDefaultPlatform().getNode();
	Teuchos::ParameterList pl = Teuchos::ParameterList();
	RCP<Node> node = rcp(new Node(pl));
        bool callFillComplete = true;

        RCP<sparse_matrix_type> tA = reader_type::readSparseFile(fileName, comm, node, callFillComplete);

        if (tA.is_null())
          throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

        RCP<TpetraCrsMatrix> tmpA1 = rcp(new TpetraCrsMatrix(tA));
        RCP<CrsMatrix>       tmpA2 = rcp_implicit_cast<CrsMatrix>(tmpA1);
        RCP<Matrix>          A     = rcp(new CrsMatrixWrap(tmpA2));

        return A;
#else
        throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
      } else {
        throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
      }
    } else {
      // Custom file format (binary)
      std::ifstream ifs(fileName.c_str(), std::ios::binary);
      TEUCHOS_TEST_FOR_EXCEPTION(!ifs.good(), Exceptions::RuntimeError, "Can not read \"" << fileName << "\"");
      int m, n, nnz;
      ifs.read(reinterpret_cast<char*>(&m),   sizeof(m));
      ifs.read(reinterpret_cast<char*>(&n),   sizeof(n));
      ifs.read(reinterpret_cast<char*>(&nnz), sizeof(nnz));

      GO indexBase = 0;
      RCP<Map>    map = MapFactory   ::Build(lib, m, indexBase, comm);
      RCP<Matrix> A   = MatrixFactory::Build(map, 1);

      TEUCHOS_TEST_FOR_EXCEPTION(sizeof(int) != sizeof(GO), Exceptions::RuntimeError, "Incompatible sizes");

      Teuchos::Array<GO> inds;
      Teuchos::Array<SC> vals;
      for (int i = 0; i < m; i++) {
        int row, rownnz;
        ifs.read(reinterpret_cast<char*>(&row),    sizeof(row));
        ifs.read(reinterpret_cast<char*>(&rownnz), sizeof(rownnz));
        inds.resize(rownnz);
        vals.resize(rownnz);
        for (int j = 0; j < rownnz; j++) {
          int index;
          ifs.read(reinterpret_cast<char*>(&index), sizeof(index));
          inds[j] = Teuchos::as<GO>(index);
        }
        for (int j = 0; j < rownnz; j++) {
          double value;
          ifs.read(reinterpret_cast<char*>(&value), sizeof(value));
          vals[j] = Teuchos::as<SC>(value);
        }
        A->insertGlobalValues(row, inds, vals);
      }
      A->fillComplete();
      return A;
    }

    return Teuchos::null;

  } //Read()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Read(const std::string&   fileName,
                                                                      const RCP<const Map> rowMap,
                                                                            RCP<const Map> colMap,
                                                                      const RCP<const Map> domainMap,
                                                                      const RCP<const Map> rangeMap,
                                                                      const bool           callFillComplete,
                                                                      const bool           tolerant,
                                                                      const bool           debug
                                                                     ) {
    TEUCHOS_TEST_FOR_EXCEPTION(rowMap.is_null(), Exceptions::RuntimeError, "Utils::Read() : rowMap cannot be null");

    RCP<const Map> domain = (domainMap.is_null() ? rowMap : domainMap);
    RCP<const Map> range  = (rangeMap .is_null() ? rowMap : rangeMap);

    const Xpetra::UnderlyingLib lib = rowMap->lib();
    if (lib == Xpetra::UseEpetra) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      Epetra_CrsMatrix *eA;
      const RCP<const Epetra_Comm> epcomm = Xpetra::toEpetra(rowMap->getComm());
      const Epetra_Map& epetraRowMap    = Map2EpetraMap(*rowMap);
      const Epetra_Map& epetraDomainMap = (domainMap.is_null() ? epetraRowMap : Map2EpetraMap(*domainMap));
      const Epetra_Map& epetraRangeMap  = (rangeMap .is_null() ? epetraRowMap : Map2EpetraMap(*rangeMap));
      int rv;
      if (colMap.is_null()) {
        rv = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), epetraRowMap, epetraRangeMap, epetraDomainMap, eA);

      } else {
        const Epetra_Map& epetraColMap  = Map2EpetraMap(*colMap);
        rv = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), epetraRowMap, epetraColMap, epetraRangeMap, epetraDomainMap, eA);
      }

      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::MatrixMarketFileToCrsMatrix return value of " + toString(rv));

      RCP<Epetra_CrsMatrix> tmpA = rcp(eA);
      RCP<Matrix>           A    = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(tmpA);

      return A;
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support.");
#endif
    } else if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type>                          reader_type;
      typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>                            map_type;

      const RCP<const map_type> tpetraRowMap    = Map2TpetraMap(*rowMap);
      RCP<const map_type>       tpetraColMap    = (colMap.is_null()    ? Teuchos::null : Map2TpetraMap(*colMap));
      const RCP<const map_type> tpetraRangeMap  = (rangeMap.is_null()  ? tpetraRowMap  : Map2TpetraMap(*rangeMap));
      const RCP<const map_type> tpetraDomainMap = (domainMap.is_null() ? tpetraRowMap  : Map2TpetraMap(*domainMap));

      RCP<sparse_matrix_type> tA = reader_type::readSparseFile(fileName, tpetraRowMap, tpetraColMap, tpetraDomainMap, tpetraRangeMap,
                                                               callFillComplete, tolerant, debug);
      if (tA.is_null())
        throw Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null.");

      RCP<TpetraCrsMatrix> tmpA1 = rcp(new TpetraCrsMatrix(tA));
      RCP<CrsMatrix>       tmpA2 = rcp_implicit_cast<CrsMatrix>(tmpA1);
      RCP<Matrix>          A     = rcp(new CrsMatrixWrap(tmpA2));

      return A;
#else
      throw Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support.");
#endif
    } else {
      throw Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra.");
    }

    return Teuchos::null;
  } //Read

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write(const std::string& fileName, const MultiVector& x) {

    std::string mapfile = "map_" + fileName;
    Write(mapfile, *(x.getMap()));

    RCP<const MultiVector> tmp_Vec = rcpFromRef(x);
#ifdef HAVE_MUELU_EPETRA
    const RCP<const EpetraMultiVector>& tmp_EVec = rcp_dynamic_cast<const EpetraMultiVector>(tmp_Vec);
    if (tmp_EVec != Teuchos::null) {
#ifdef HAVE_MUELU_EPETRAEXT
      int rv = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(), *(tmp_EVec->getEpetra_MultiVector()));
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::RowMatrixToMatrixMarketFile return value of " + toString(rv));
#else
      throw Exceptions::RuntimeError("Compiled without EpetraExt");
#endif
      return;
    }
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraMultiVector> &tmp_TVec = rcp_dynamic_cast<const TpetraMultiVector>(tmp_Vec);
    if (tmp_TVec != Teuchos::null) {
      RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > TVec = tmp_TVec->getTpetra_MultiVector();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeDenseFile(fileName, TVec);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in matrix writing");

  } //Write

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write(const std::string& fileName, const Map& M) {
    RCP<const Map> tmp_Map = rcpFromRef(M);
#ifdef HAVE_MUELU_EPETRAEXT
    const RCP<const Xpetra::EpetraMap>& tmp_EMap = rcp_dynamic_cast<const Xpetra::EpetraMap>(tmp_Map);
    if (tmp_EMap != Teuchos::null) {
#ifdef HAVE_MUELU_EPETRAEXT
      int rv = EpetraExt::BlockMapToMatrixMarketFile(fileName.c_str(), tmp_EMap->getEpetra_Map());
      if (rv != 0)
        throw Exceptions::RuntimeError("EpetraExt::BlockMapToMatrixMarketFile() return value of " + toString(rv));
#else
      throw(Exceptions::RuntimeError("Compiled without EpetraExt"));
#endif
      return;
    }
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraMap> &tmp_TMap = rcp_dynamic_cast<const TpetraMap>(tmp_Map);
    if (tmp_TMap != Teuchos::null) {
      RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> > TMap = tmp_TMap->getTpetra_Map();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeMapFile(fileName, *TMap);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in matrix writing");

  } //Write

#include <unistd.h>

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PauseForDebugger() {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    int myPID = comm->getRank();
    int pid   = getpid();

    char hostname[80];
    for (int i = 0; i <comm->getSize(); i++) {
      if (i == myPID) {
        gethostname(hostname, sizeof(hostname));
        std::cout << "Host: " << hostname << "\tMPI rank: " << myPID << ",\tPID: " << pid << "\n\tattach " << pid << std::endl;
        sleep(1);
      }
    }

    if (myPID == 0) {
      std::cout << "** Enter a character to continue > " << std::endl;
      char go = ' ';
      scanf("%c", &go);
    }
    comm->barrier();
  } //PauseForDebugger

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PowerMethod(const Matrix& A, bool scaleByDiag,
                                                                                    LO niters, Magnitude tolerance, bool verbose, unsigned int seed)
  {
    if (!(A.getRangeMap()->isSameAs(*(A.getDomainMap()))))
      throw Exceptions::Incompatible("Utils::PowerMethod: operator must have domain and range maps that are equivalent.");

    // Create three vectors, fill z with random numbers
    RCP<Vector> q     = VectorFactory::Build(A.getDomainMap());
    RCP<Vector> qinit = VectorFactory::Build(A.getDomainMap());
    RCP<Vector> r     = VectorFactory::Build(A.getRangeMap());
    RCP<Vector> z     = VectorFactory::Build(A.getRangeMap());
    z->setSeed(seed);  // seed random number generator
    z->randomize(true);// use Xpetra implementation: -> same results for Epetra and Tpetra

    Teuchos::Array<Magnitude> norms(1);

    typedef Teuchos::ScalarTraits<Scalar> STS;

    const Scalar zero = STS::zero();
    const Scalar one  = STS::one();

    Scalar lambda = zero;
    Magnitude residual = STS::magnitude(zero);

    // power iteration
    RCP<Vector> diagVec, oneOverDiagonal;
    if (scaleByDiag) {
      diagVec = VectorFactory::Build(A.getRowMap());
      A.getLocalDiagCopy(*diagVec);
      oneOverDiagonal = VectorFactory::Build(A.getRowMap());
      oneOverDiagonal->reciprocal(*diagVec);
    }

    for (int iter = 0; iter < niters; ++iter) {
      z->norm2(norms);                               // Compute 2-norm of z
      q->update(one / norms[0],*z,zero);                 // Set q = z / normz
      A.apply(*q, *z);                               // Compute z = A*q
      if (scaleByDiag) z->elementWiseMultiply(one, *oneOverDiagonal, *z, zero);
      lambda = q->dot(*z);                            // Approximate maximum eigenvalue: lamba = dot(q,z)
      if (iter % 100 == 0 || iter + 1 == niters) {
        r->update(1.0, *z, -lambda, *q, zero);         // Compute A*q - lambda*q
        r->norm2(norms);
        residual = STS::magnitude(norms[0] / lambda);
        if (verbose) {
          std::cout << "Iter = " << iter
                    << "  Lambda = " << lambda
                    << "  Residual of A*q - lambda*q = " << residual
                    << std::endl;
        }
      }
      if (residual < tolerance)
        break;
    }

    return lambda;
  } //PowerMethod

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix(Matrix& Op, const Teuchos::ArrayRCP<const SC>& scalingVector, bool doInverse,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
    SC one = Teuchos::ScalarTraits<SC>::one();
    Teuchos::ArrayRCP<SC> sv(scalingVector.size());
    if (doInverse) {
      for (int i = 0; i < scalingVector.size(); ++i)
        sv[i] = one / scalingVector[i];
    } else {
      for (int i = 0; i < scalingVector.size(); ++i)
        sv[i] = scalingVector[i];
    }

    switch (Op.getRowMap()->lib()) {
      case Xpetra::UseTpetra:
        MyOldScaleMatrix_Tpetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;

      case Xpetra::UseEpetra:
        Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;

      default:
        throw Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled.");
        break;
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Tpetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
#ifdef HAVE_MUELU_TPETRA
    try {
      Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& tpOp = Op2NonConstTpetraCrs(Op);

      const RCP<const Tpetra::Map<LO,GO,NO> > rowMap    = tpOp.getRowMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > domainMap = tpOp.getDomainMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > rangeMap  = tpOp.getRangeMap();

      size_t maxRowSize = tpOp.getNodeMaxNumRowEntries();
      if (maxRowSize == Teuchos::as<size_t>(-1)) // hasn't been determined yet
        maxRowSize = 20;

      std::vector<SC> scaledVals(maxRowSize);
      if (tpOp.isFillComplete())
        tpOp.resumeFill();

      if (Op.isLocallyIndexed() == true) {
        Teuchos::ArrayView<const LO> cols;
        Teuchos::ArrayView<const SC> vals;

        for (size_t i = 0; i < rowMap->getNodeNumElements(); ++i) {
          tpOp.getLocalRowView(i, cols, vals);
          size_t nnz = tpOp.getNumEntriesInLocalRow(i);
          if (nnz > maxRowSize) {
            maxRowSize = nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j = 0; j < nnz; ++j)
            scaledVals[j] = vals[j]*scalingVector[i];

          if (nnz > 0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0], nnz);
            tpOp.replaceLocalValues(i, cols, valview);
          }
        } //for (size_t i=0; ...

      } else {
        Teuchos::ArrayView<const GO> cols;
        Teuchos::ArrayView<const SC> vals;

        for (size_t i = 0; i < rowMap->getNodeNumElements(); ++i) {
          GO gid = rowMap->getGlobalElement(i);
          tpOp.getGlobalRowView(gid, cols, vals);
          size_t nnz = tpOp.getNumEntriesInGlobalRow(gid);
          if (nnz > maxRowSize) {
            maxRowSize = nnz;
            scaledVals.resize(maxRowSize);
          }
          // FIXME FIXME FIXME FIXME FIXME FIXME
          for (size_t j = 0; j < nnz; ++j)
            scaledVals[j] = vals[j]*scalingVector[i]; //FIXME i or gid?

          if (nnz > 0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0], nnz);
            tpOp.replaceGlobalValues(gid, cols, valview);
          }
        } //for (size_t i=0; ...
      }

      if (doFillComplete) {
        if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
          throw Exceptions::RuntimeError("In Utils::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined");

        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage",    doOptimizeStorage);
        params->set("No Nonlocal Changes", true);
        Op.fillComplete(Op.getDomainMap(), Op.getRangeMap(), params);
      }
    } catch(...) {
      throw Exceptions::RuntimeError("Only Tpetra::CrsMatrix types can be scaled (Err.1)");
    }
#else
    throw Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been enabled.");
#endif
  } //MyOldScaleMatrix_Tpetra()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Teuchos::FancyOStream> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MakeFancy(std::ostream& os) {
    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(os));
    return fancy;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Distance2(const MultiVector& v, LocalOrdinal i0, LocalOrdinal i1) {
    size_t numVectors = v.getNumVectors();

    Scalar d = Teuchos::ScalarTraits<Scalar>::zero();
    for (size_t j = 0; j < numVectors; j++) {
      Teuchos::ArrayRCP<const Scalar> vv = v.getData(j);
      d += (vv[i0] - vv[i1])*(vv[i0] - vv[i1]);
    }

    return Teuchos::ScalarTraits<SC>::magnitude(d);
  }

  template <class SC, class LO, class GO, class NO, class LMO>
  ArrayRCP<const bool> Utils<SC, LO, GO, NO, LMO>::DetectDirichletRows(const Matrix& A, const typename Teuchos::ScalarTraits<SC>::magnitudeType& tol) {
    LO numRows = A.getNodeNumRows();

    typedef Teuchos::ScalarTraits<SC> STS;

    ArrayRCP<bool> boundaryNodes(numRows, true);
    for (LO row = 0; row < numRows; row++) {
      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      A.getLocalRowView(row, indices, vals);

      size_t nnz = A.getNumEntriesInLocalRow(row);
      if (nnz > 1)
        for (size_t col = 0; col < nnz; col++)
          if ( (indices[col] != row) && STS::magnitude(vals[col]) > tol) {
            boundaryNodes[row] = false;
            break;
          }
    }

    return boundaryNodes;
  }

  template <class SC, class LO, class GO, class NO, class LMO>
  std::string Utils<SC, LO, GO, NO, LMO>::PrintMatrixInfo(const Matrix& A, const std::string& msgTag, RCP<const ParameterList> params) {
    std::ostringstream ss;
    ss << msgTag << " size =  " << A.getGlobalNumRows() << " x " << A.getGlobalNumCols() << ", nnz = " << A.getGlobalNumEntries() << std::endl;

    if (params.is_null())
      return ss.str();

    if (params->isParameter("printLoadBalancingInfo") && params->get<bool>("printLoadBalancingInfo")) {
      RCP<const Teuchos::Comm<int> > comm = A.getRowMap()->getComm();
      GO numActiveProcesses = comm->getSize(), numProcessesWithData = 0;

      // aggregate data
      GO  numMyNnz = A.getNodeNumEntries(),     minNnz,     maxNnz;
      GO numMyRows = A.getNodeNumRows(),    minNumRows, maxNumRows;
      double  numMyNnz2 =  Teuchos::as<double>(numMyNnz)* numMyNnz,     sumNnz = Teuchos::as<double>(A.getGlobalNumEntries()),     sum2Nnz;
      double numMyRows2 = Teuchos::as<double>(numMyRows)*numMyRows, sumNumRows = Teuchos::as<double>(A.getGlobalNumRows()),    sum2NumRows;
      maxAll(comm,                                       numMyNnz, maxNnz);
      maxAll(comm,                                      numMyRows, maxNumRows);
      sumAll(comm, (GO)((numMyRows > 0) ?         1 :          0), numProcessesWithData);
      minAll(comm, (GO)(( numMyNnz > 0) ?  numMyNnz :     maxNnz), minNnz);
      minAll(comm, (GO)((numMyRows > 0) ? numMyRows : maxNumRows), minNumRows);
      sumAll(comm,                                      numMyNnz2, sum2Nnz);
      sumAll(comm,                                     numMyRows2, sum2NumRows);

      double avgNumRows = sumNumRows / numProcessesWithData;
      double avgNnz     = sumNnz     / numProcessesWithData;
      double devNumRows = (numProcessesWithData != 1 ? sqrt((sum2NumRows - sumNumRows*sumNumRows/numProcessesWithData)/(numProcessesWithData-1)) : 0);
      double devNnz     = (numProcessesWithData != 1 ? sqrt((sum2Nnz     -     sumNnz*    sumNnz/numProcessesWithData)/(numProcessesWithData-1)) : 0);

      char buf[256];
      ss << msgTag << " Load balancing info:" << std::endl;
      ss << msgTag << "   # active processes: "   << numActiveProcesses << ",  # processes with data = " << numProcessesWithData << std::endl;
      sprintf(buf, "avg = %.2e,  dev = %4.1f%%,  min = %+5.1f%%,  max = %+5.1f%%", avgNumRows,
              (devNumRows/avgNumRows)*100, (minNumRows/avgNumRows-1)*100, (maxNumRows/avgNumRows-1)*100);
      ss << msgTag << "   # rows per proc   : " << buf << std::endl;
      sprintf(buf, "avg = %.2e,  dev = %4.1f%%,  min = %+5.1f%%,  max = %+5.1f%%", avgNnz,
              (devNnz/avgNnz)*100, (minNnz/avgNnz-1)*100, (maxNnz/avgNnz-1)*100);
      ss << msgTag << "   #  nnz per proc   : " << buf << std::endl;
    }

    return ss.str();
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Transpose(Matrix& Op, bool optimizeTranspose) {
   Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("YY Entire Transpose"));
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    try {
      const Epetra_CrsMatrix& epetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(Op);
      (void) epetraOp; // silence "unused variable" compiler warning
    } catch (...) {
      TorE = "tpetra";
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    if (TorE == "tpetra") {
      try {
        const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& tpetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(Op);

        RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
        {
          Teuchos::TimeMonitor tmm(*Teuchos::TimeMonitor::getNewTimer("YY Tpetra Transpose Only"));
          Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> transposer(rcpFromRef(tpetraOp)); //more than meets the eye
          A = transposer.createTranspose();
        }

        RCP<TpetraCrsMatrix> AA = rcp(new TpetraCrsMatrix(A) );
        RCP<CrsMatrix> AAA = rcp_implicit_cast<CrsMatrix>(AA);
        RCP<Matrix> AAAA = rcp( new CrsMatrixWrap(AAA) );
        if (!AAAA->isFillComplete())
          AAAA->fillComplete(Op.getRangeMap(),Op.getDomainMap());

        return AAAA;

      } catch (...) {
        throw Exceptions::RuntimeError("Utils::Transpose: Can only transpose Crs matrices");
      }
    } //if
#endif

    // Epetra case
    std::cout << "Utilities::Transpose() not implemented for Epetra" << std::endl;
    return Teuchos::null;

  } // Transpose

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Epetra(Matrix& Op, const Teuchos::ArrayRCP<SC>& scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MyOldScalematrix and Epetra cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(const Matrix& A, bool transposeA, SC alpha, Matrix& B, SC beta) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (A.getRowMap()->lib() == Xpetra::UseEpetra) {
      throw Exceptions::RuntimeError("You cannot use Epetra::MatrixMatrix::Add with Scalar!=double or Ordinal!=int");

    } else if (A.getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& tpA = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(A);
            Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& tpB = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, beta);
#else
      throw Exceptions::RuntimeError("MueLu must be compiled with Tpetra.");
#endif
    }
  } //Utils2::TwoMatrixAdd()


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(const Matrix& A, bool transposeA, SC const &alpha,
                           const Matrix& B, bool transposeB, SC const &beta,
                           RCP<Matrix>& C, Teuchos::FancyOStream & fos, bool AHasFixedNnzPerRow) {
    if (!(A.getRowMap()->isSameAs(*(B.getRowMap()))))
      throw Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same.");

    if (C == Teuchos::null) {
      if (!A.isFillComplete() || !B.isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Global statistics are not available for estimates.");

      size_t maxNzInA     = A.getGlobalMaxNumRowEntries();
      size_t maxNzInB     = B.getGlobalMaxNumRowEntries();
      size_t numLocalRows = A.getNodeNumRows();

      if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
        // first check if either A or B has at most 1 nonzero per row
        // the case of both having at most 1 nz per row is handled by the ``else''
        Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);

        if ((maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = B.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;

        } else {
          for (size_t i = 0; i < numLocalRows; ++i)
            exactNnzPerRow[i] = A.getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
        }

        fos << "Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
             << ", using static profiling" << std::endl;
        C = rcp(new CrsMatrixWrap(A.getRowMap(), exactNnzPerRow, Xpetra::StaticProfile));

      } else {
        // general case
        double nnzPerRowInA = Teuchos::as<double>(A.getGlobalNumEntries()) / A.getGlobalNumRows();
        double nnzPerRowInB = Teuchos::as<double>(B.getGlobalNumEntries()) / B.getGlobalNumRows();
        LO    nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

        LO maxPossible = A.getGlobalMaxNumRowEntries() + B.getGlobalMaxNumRowEntries();
        //Use static profiling (more efficient) if the estimate is at least as big as the max
        //possible nnz's in any single row of the result.
        Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

        fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
        fos << "Utils::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
             << ", max possible nnz per row in sum = " << maxPossible
             << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
             << std::endl;
        C = rcp(new CrsMatrixWrap(A.getRowMap(), nnzToAllocate, pft));
      }
      if (transposeB)
        fos << "Utils::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
      throw Exceptions::RuntimeError("You cannot use Epetra::MatrixMatrix::Add with Scalar!=double or Ordinal!=int");

    } else if (C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& tpA = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(A);
      const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& tpB = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(B);
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >  tpC = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Add(tpA, transposeA, alpha, tpB, transposeB, beta, tpC);
#else
      throw Exceptions::RuntimeError("MueLu must be compile with Tpetra.");
#endif
    }

    ///////////////////////// EXPERIMENTAL
    if (A.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(A));
    if (B.IsView("stridedMaps")) C->CreateView("stridedMaps", rcpFromRef(B));
    ///////////////////////// EXPERIMENTAL

  } //Utils2::TwoMatrixAdd()

} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DEF_HPP

//  LocalWords:  LocalOrdinal
