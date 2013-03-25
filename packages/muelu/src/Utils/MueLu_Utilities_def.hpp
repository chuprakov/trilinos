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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
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
#include <Epetra_RowMatrixTransposer.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <Xpetra_EpetraUtils.hpp>
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
#include <TpetraExt_MatrixMatrix.hpp>
#include <Tpetra_RowMatrixTransposer.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Xpetra_TpetraMultiVector.hpp>
#include <Xpetra_TpetraCrsMatrix.hpp>
#endif // HAVE_MUELU_TPETRA

#include <Xpetra_Map.hpp>
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

#define scan(rcpComm, in, out)                                        \
  Teuchos::scan(*rcpComm, Teuchos::REDUCE_SUM, in, Teuchos::outArg(out));

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
  RCP<const Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2EpetraMV(RCP<MultiVector> const Vec) {
    //rcp<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    RCP<const EpetraMultiVector > tmpVec;
    tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed"));
    RCP<const Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
    return epVec;
  } //MV2EpetraMV

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_MultiVector> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstEpetraMV(RCP<MultiVector> Vec) {
    RCP<const EpetraMultiVector> tmpVec = rcp_dynamic_cast<EpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::EpetraMultiVector failed"));
    RCP<Epetra_MultiVector> epVec = tmpVec->getEpetra_MultiVector();
    return epVec;
  } //MV2EpetraMV

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Epetra_MultiVector& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstEpetraMV(MultiVector &Vec) {
    EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
    RCP<Epetra_MultiVector> epVec = tmpVec.getEpetra_MultiVector();
    return *epVec;
  } //MV2EpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Epetra_MultiVector const& Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2EpetraMV(MultiVector const &Vec) {
    EpetraMultiVector const &tmpVec = dynamic_cast<EpetraMultiVector const&>(Vec);
    RCP<Epetra_MultiVector const> epVec = tmpVec.getEpetra_MultiVector();
    return *epVec;
  } //MV2EpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2EpetraCrs(RCP<const Matrix> Op) {
    RCP<const Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getEpetra_CrsMatrix();
    return A;
  } //Op2EpetraCrs


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(RCP<Matrix> Op) {
    RCP<Epetra_CrsMatrix> A;
    // Get the underlying Epetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::EpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getEpetra_CrsMatrixNonConst();
    return A;
  } //Op2NonConstEpetraCrs
#endif

#ifdef HAVE_MUELU_TPETRA

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2TpetraMV(RCP<MultiVector> const Vec) {
    //rcp<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    RCP<const TpetraMultiVector > tmpVec;
    tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed"));
    RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
    return tpVec;
  } //MV2TpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV(RCP<MultiVector> Vec) {
    RCP<const TpetraMultiVector> tmpVec = rcp_dynamic_cast<TpetraMultiVector>(Vec);
    if (tmpVec == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::MultiVector to Xpetra::TpetraMultiVector failed"));
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec->getTpetra_MultiVector();
    return tpVec;
  } //MV2NonConstTpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> & Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV(MultiVector &Vec) {
    TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
    return *tpVec;
  } //MV2NonConstTpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2NonConstTpetraMV2(MultiVector &Vec) {
    TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
    RCP<Tpetra::MultiVector<SC,LO,GO,NO> > tpVec = tmpVec.getTpetra_MultiVector();
    return tpVec;
  } //MV2NonConstTpetraMV2


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>  const&
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MV2TpetraMV(MultiVector const &Vec) {
    TpetraMultiVector const &tmpVec = dynamic_cast<TpetraMultiVector const&>(Vec);
    RCP<Tpetra::MultiVector<SC,LO,GO,NO>  const> tpVec = tmpVec.getTpetra_MultiVector();
    return *tpVec;
  } //MV2TpetraMV


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(RCP<Matrix> Op) {
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getTpetra_CrsMatrix();
    return A;
  } //Op2TpetraCrs


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstTpetraCrs(RCP<Matrix> Op) {
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
    // Get the underlying Tpetra Mtx
    RCP<const CrsMatrixWrap> crsOp = rcp_dynamic_cast<const CrsMatrixWrap>(Op);
    if (crsOp == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed"));
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp->getCrsMatrix();
    const RCP<const TpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null)
      throw(Exceptions::BadCast("Cast from Xpetra::CrsMatrix to Xpetra::TpetraCrsMatrix failed"));
    A = tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();
    return A;
  } //Op2NonConstTpetraCrs

#endif

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Multiply(const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& A,
                                                                          bool transposeA,
                                                                          const Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>& B,
                                                                          bool transposeB,
                                                                          RCP< Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > C_in,
                                                                          bool doFillComplete,
                                                                          bool doOptimizeStorage) {

    RCP<Matrix> C;

    // Preconditions
    if (!A.isFillComplete())
      throw(Exceptions::RuntimeError("A is not fill-completed"));
    if (!B.isFillComplete())
      throw(Exceptions::RuntimeError("B is not fill-completed"));

    // Optimization using ML Multiply when available
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT) && defined(HAVE_MUELU_ML)
    if (B.getDomainMap()->lib() == Xpetra::UseEpetra && !transposeA && !transposeB) {

      if (!transposeA && !transposeB) {
        RCP<const Epetra_CrsMatrix> epA = Op2EpetraCrs(rcpFromRef(A)); // TODO: do conversion without RCPs.
        RCP<const Epetra_CrsMatrix> epB = Op2EpetraCrs(rcpFromRef(B));
        RCP<Epetra_CrsMatrix> epC = MLTwoMatrixMultiply(*epA, *epB);
        C = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(epC);
        if(doFillComplete) {
          RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
          params->set("Optimize Storage",doOptimizeStorage);
          C->fillComplete(B.getDomainMap(), A.getRangeMap(), params);
        }
      }

      C->CreateView("stridedMaps", rcpFromRef(A), transposeA, rcpFromRef(B), transposeB);
      return C;

    }
# endif // EPETRA + EPETRAEXT + ML

    // Default case: Xpetra Multiply

    C = C_in;

    if (C == Teuchos::null) {
      double nnzPerRow=Teuchos::as<double>(0);
      if (A.getDomainMap()->lib() == Xpetra::UseTpetra) {
        // for now, follow what ML and Epetra do.
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
        RCP<Teuchos::FancyOStream> fos = MakeFancy(std::cout);
        fos->setOutputToRootOnly(0);
        *fos << "Utils::Multiply : Estimate for nnz per row of product matrix = " << Teuchos::as<LO>(nnzPerRow) << std::endl;
      }

      if (transposeA) C = MatrixFactory::Build(A.getDomainMap(), Teuchos::as<LO>(nnzPerRow));
      else            C = MatrixFactory::Build(A.getRowMap(),    Teuchos::as<LO>(nnzPerRow));
    } else {
      C->resumeFill(); // why this is not done inside of Tpetra MxM?
      std::cout << "Reuse C pattern" << std::endl;
    }

    Xpetra::MatrixMatrix::Multiply(A, transposeA, B, transposeB, *C, doFillComplete, doOptimizeStorage);

    // fill strided maps information
    // this is necessary since the ML matrix matrix multiplication routine
    // has no handling for this
    // TODO: move this call to MLMultiply...
    C->CreateView("stridedMaps", rcpFromRef(A), transposeA, rcpFromRef(B), transposeB);

    return C;

  } // Multiply()

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MLTwoMatrixMultiply(const Epetra_CrsMatrix& epA,
      const Epetra_CrsMatrix& epB)
  {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_ML)
    ML_Comm* comm;
    ML_Comm_Create(&comm);
    if (comm->ML_mypid == 0)
      std::cout << "****** USING ML's MATRIX MATRIX MULTIPLY (LNM version) ******" << std::endl;
#           ifdef HAVE_MPI
    // ML_Comm uses MPI_COMM_WORLD, so try to use the same communicator as epA.
    const Epetra_MpiComm * Mcomm=dynamic_cast<const Epetra_MpiComm*>(&(epA.Comm()));
    if(Mcomm) ML_Comm_Set_UsrComm(comm,Mcomm->GetMpiComm());
#           endif
    //in order to use ML, there must be no indices missing from the matrix column maps.
    EpetraExt::CrsMatrix_SolverMap Atransform;
    EpetraExt::CrsMatrix_SolverMap Btransform;
    const Epetra_CrsMatrix& A = Atransform(const_cast<Epetra_CrsMatrix&>(epA));
    const Epetra_CrsMatrix& B = Btransform(const_cast<Epetra_CrsMatrix&>(epB));

    if (!A.Filled())    throw(Exceptions::RuntimeError("A has to be FillCompeleted"));
    if (!B.Filled())    throw(Exceptions::RuntimeError("B has to be FillCompeleted"));

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
    for (int i=0; i<getrow_comm->N_neighbors; i++)
      {
        N_rcvd += (getrow_comm->neighbors)[i].N_rcv;
        N_send += (getrow_comm->neighbors)[i].N_send;
        if ( ((getrow_comm->neighbors)[i].N_rcv !=0) &&
             ((getrow_comm->neighbors)[i].rcv_list != NULL) ) flag = 1;
      }
    // For some unknown reason, ML likes to have stuff one larger than
    // neccessary...
    std::vector<double> dtemp(N_local + N_rcvd + 1); // "double" vector for comm function
    std::vector<int>    cmap (N_local + N_rcvd + 1); // vector for gids
    for (int i=0; i<N_local; ++i)
      {
        cmap[i] = B.DomainMap().GID(i);
        dtemp[i] = (double) cmap[i];
      }
    ML_cheap_exchange_bdry(&dtemp[0],getrow_comm,N_local,N_send,comm_AB); // do communication
    if (flag) // process received data
      {
        int count = N_local;
        const int neighbors = getrow_comm->N_neighbors;
        for (int i=0; i< neighbors; i++)
          {
            const int nrcv = getrow_comm->neighbors[i].N_rcv;
            for (int j=0; j<nrcv; j++)
              cmap[getrow_comm->neighbors[i].rcv_list[j]] = (int) dtemp[count++];
          }
      }
    else
      for (int i=0; i<N_local+N_rcvd; ++i) cmap[i] = (int)dtemp[i];
    dtemp.clear();  // free double array

    // we can now determine a matching column map for the result
    Epetra_Map gcmap(-1,N_local+N_rcvd,&cmap[0],B.ColMap().IndexBase(),A.Comm());

    int allocated=0;
    int rowlength;
    double* val=NULL;
    int* bindx=NULL;
    const int myrowlength = A.RowMap().NumMyElements();
    const Epetra_Map& rowmap = A.RowMap();

    // determine the maximum bandwith for the result matrix.
    // replaces the old, very(!) memory-consuming guess:
    // int guessnpr = A.MaxNumEntries()*B.MaxNumEntries();
    int educatedguess = 0;
    for (int i=0; i<myrowlength; ++i)
      {
        // get local row
        ML_get_matrix_row(ml_AtimesB,1,&i,&allocated,&bindx,&val,&rowlength,0);
        if (rowlength>educatedguess) educatedguess = rowlength;
      }

    // allocate our result matrix and fill it
    RCP<Epetra_CrsMatrix> result
      = rcp(new Epetra_CrsMatrix(::Copy,A.RangeMap(),gcmap,educatedguess,false));

    std::vector<int> gcid(educatedguess);
    for (int i=0; i<myrowlength; ++i)
      {
        const int grid = rowmap.GID(i);
        // get local row
        ML_get_matrix_row(ml_AtimesB,1,&i,&allocated,&bindx,&val,&rowlength,0);
        if (!rowlength) continue;
        if ((int)gcid.size() < rowlength) gcid.resize(rowlength);
        for (int j=0; j<rowlength; ++j)
          {
            gcid[j] = gcmap.GID(bindx[j]);
            if (gcid[j]<0) throw(Exceptions::RuntimeError("Error: cannot find gcid!"));
          }
        int err = result->InsertGlobalValues(grid,rowlength,val,&gcid[0]);
        if (err!=0 && err!=1) throw(Exceptions::RuntimeError("Epetra_CrsMatrix::InsertGlobalValues returned err="+err));
      }
    // free memory
    if (bindx) ML_free(bindx);
    if (val) ML_free(val);
    ML_Operator_Destroy(&ml_AtimesB);
    ML_Comm_Destroy(&comm);

    return result;
#else // no MUELU_ML
    TEUCHOS_TEST_FOR_EXCEPTION( true, Xpetra::Exceptions::RuntimeError,
                                "HAVE_MUELU_ML compiler flag not set. no ML multiply available." );
    return Teuchos::null;
#endif
  }
#endif //ifdef HAVE_MUELU_EPETRAEXT

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixMultiplyBlock(RCP<BlockedCrsMatrix> const &A, bool transposeA,
                                                        RCP<Xpetra::BlockedCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &B, bool transposeB,
                                                        bool doFillComplete,
                                                        bool doOptimizeStorage)
  {
    if(transposeA || transposeB)
      throw(Exceptions::RuntimeError("TwoMatrixMultiply for BlockedCrsMatrix not implemented for transposeA==true or transposeB==true"));

    // todo make sure that A and B are filled and completed

    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > rgmapextractor = A->getRangeMapExtractor();
    RCP<const Xpetra::MapExtractor<Scalar, LocalOrdinal, GlobalOrdinal, Node> > domapextractor = B->getDomainMapExtractor();

    RCP<BlockedCrsMatrix> C = rcp(new BlockedCrsMatrix(rgmapextractor,
                                                           domapextractor,
                                                           33 /* TODO fix me */));

    // loop over all block rows of A
    for(size_t i=0; i<A->Rows(); ++i)
      {
        // loop over all block columns of B
        for(size_t j=0; j<B->Cols(); ++j)
          {
            // empty CrsMatrixWrap
            RCP<Matrix> Cij = MatrixFactory::Build(A->getRangeMap(i), 33 /* TODO fix me */);

            // loop for calculating entry C_{ij}
            for(size_t l=0; l<B->Rows(); ++l)
              {
                RCP<CrsMatrix> crmat1 = A->getMatrix(i,l);
                RCP<CrsMatrix> crmat2 = B->getMatrix(l,j);
                RCP<CrsMatrixWrap> crop1 = rcp(new CrsMatrixWrap(crmat1));
                RCP<CrsMatrixWrap> crop2 = rcp(new CrsMatrixWrap(crmat2));

                RCP<Matrix> temp = MueLu::Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Multiply(*crop1, false, *crop2, false);

                // sum up
                MueLu::Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(temp, false, 1.0, Cij, 1.0);
              }

            Cij->fillComplete(B->getDomainMap(j), A->getRangeMap(i));

            RCP<CrsMatrixWrap> crsCij = Teuchos::rcp_dynamic_cast<CrsMatrixWrap>(Cij);
            TEUCHOS_TEST_FOR_EXCEPTION( Cij==Teuchos::null, Xpetra::Exceptions::BadCast,
                                        "MatrixFactory failed in generating a CrsMatrixWrap." );

            RCP<CrsMatrix> crsMatCij = crsCij->getCrsMatrix();
            C->setMatrix(i,j,crsMatCij);

          }
      }

    if(doFillComplete)
      C->fillComplete();  // call default fillComplete for BlockCrsMatrixWrap objects

    return C;
  } // TwoMatrixMultiplyBlock


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MatrixPrint(RCP<Matrix> const &Op) {
    std::string label = "unlabeled operator";
    MatrixPrint(Op, label);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MatrixPrint(RCP<Matrix> const &Op, std::string const &label) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    RCP<const Epetra_CrsMatrix> epOp = Op2EpetraCrs(Op);
    int mypid = epOp->RowMap().Comm().MyPID();
    if (mypid == 0)
      std::cout << "\n===============\n" << label << "\n==============" << std::endl;

    if (mypid == 0) std::cout << "\n   -- row map -- \n" << std::endl;
    std::cout << epOp->RowMap() << std::endl;
    sleep(1);
    epOp->RowMap().Comm().Barrier();

    if (mypid == 0) std::cout << "\n   -- column map -- \n" << std::endl;
    std::cout << epOp->ColMap() << std::endl;
    sleep(1);
    epOp->RowMap().Comm().Barrier();

    std::cout << *epOp << std::endl;
#endif
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildMatrixDiagonal(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &A)
  {
    const RCP<const Map> rowmap = A->getRowMap();
    //std::vector<SC> diag(A->getNodeNumRows());
    std::vector<SC> diag(rowmap->getNodeNumElements());
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    //for (size_t i=0; i<A->getNodeNumRows(); ++i)
    for (size_t i=0; i<rowmap->getNodeNumElements(); ++i) {
      A->getLocalRowView(i,cols,vals);
      //for (Teuchos::ArrayView<const LO>::size_type j=0; j<cols.size(); j++)
      for (size_t j=0; j<Teuchos::as<size_t>(cols.size()); j++) { // TODO: cleanup LO vs size_t
        //TODO this will break down if diagonal entry is not present
        //if (!(cols[j] > i)) //JG says this will work ... maybe
        if (cols[j] == Teuchos::as<LO>(i)) {  // TODO: cleanup LO vs size_t
          diag[i] = vals[j];
          break;
        }
      }
    }

    RCP< Matrix > D = rcp( new CrsMatrixWrap(rowmap, 1) );
    std::vector<GO> diagInd(1);
    Teuchos::ArrayView<GO> iv(&diagInd[0],1);
    //for (size_t i=0; i< A->getNodeNumRows(); ++i)
    for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
      Teuchos::ArrayView<SC> av(&diag[i],1);
      diagInd[0] = rowmap->getGlobalElement(i);
      D->insertGlobalValues(i,iv,av);
    }
    D->fillComplete();
    //MatrixPrint(D);

    return D;

  } //BuildMatrixDiagonal()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<Scalar> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMatrixDiagonal(const Matrix &A)
  {
    const RCP<const Map> rowmap = A.getRowMap();
    size_t locSize = rowmap->getNodeNumElements();
    Teuchos::ArrayRCP<SC> diag(locSize);
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i=0; i<locSize; ++i) {
      A.getLocalRowView(i,cols,vals);
      for (LO j=0; j<cols.size(); ++j) {
        //TODO this will break down if diagonal entry is not present
        //if (!(cols[j] > i))   //JG says this will work ... maybe
        if (Teuchos::as<size_t>(cols[j]) == i) {
          diag[i] = vals[j];
          break;
        }
      }
    }
    //for (int i=0; i<locSize; ++i) std::cout << "diag[" << i << "] = " << diag[i] << std::endl;
    return diag;
  } //GetMatrixDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<Scalar> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetLumpedMatrixDiagonal(const Matrix &A)
  {
    const RCP<const Map> rowmap = A.getRowMap();
    size_t locSize = rowmap->getNodeNumElements();
    Teuchos::ArrayRCP<SC> diag(locSize);
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    for (size_t i=0; i<locSize; ++i) { // loop over rows
      A.getLocalRowView(i,cols,vals);
      Scalar absRowSum = Teuchos::ScalarTraits<Scalar>::zero();
      for (LO j=0; j<cols.size(); ++j) { // loop over cols
        absRowSum += Teuchos::ScalarTraits<Scalar>::magnitude(vals[j]);
      }
      diag[i] = absRowSum;
    }
    //for (int i=0; i<locSize; ++i) std::cout << "diag[" << i << "] = " << diag[i] << std::endl;
    return diag;
  } //GetMatrixDiagonal

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMatrixOverlappedDiagonal(const Matrix &A)
  {
    Teuchos::ArrayRCP<SC> diagVals = GetMatrixDiagonal(A);  //FIXME should this return a Vector instead?
    RCP<Vector> diagonal = VectorFactory::Build(A.getColMap());
    RCP<Vector> localDiag = VectorFactory::Build(A.getRowMap());
    ArrayRCP<SC> localDiagVals = localDiag->getDataNonConst(0);
    for (LO i=0; i<localDiagVals.size(); ++i)
      localDiagVals[i] = diagVals[i];
    localDiagVals = null;  //release view
    diagVals = null;
    //TODO there's a problem with the importer from the underlying Tpetra::CrsGraph
    //TODO so right now construct an importer.
    //diagonal->doImport(*localDiag,*(A.getCrsGraph()->getImporter()),Xpetra::INSERT);
    RCP<const Import> importer = ImportFactory::Build(A.getRowMap(), A.getColMap());
    diagonal->doImport(*localDiag,*importer,Xpetra::INSERT);
    return diagonal;
  } //GetMatrixOverlappedDiagonal


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ScaleMatrix(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector, bool doInverse)
  {
#ifdef HAVE_MUELU_TPETRA
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpOp;
    try {
      tpOp = Op2NonConstTpetraCrs(Op);
    }
    catch(...) {
      throw(Exceptions::RuntimeError("Sorry, haven't implemented matrix scaling for epetra"));
    }
    Tpetra::Vector<SC,LO,GO,NO> x(tpOp->getRowMap(),scalingVector());
    if(doInverse){
      Tpetra::Vector<SC,LO,GO,NO> xi(tpOp->getRowMap());
      xi.reciprocal(x);
      tpOp->leftScale(xi);
    }
    else
      tpOp->leftScale(x);
#else
    throw(Exceptions::RuntimeError("Sorry, haven't implemented matrix scaling for epetra"));
#endif // HAVE_MUELU_TPETRA
  } //ScaleMatrix()

#ifdef UNUSED // and does not work with SC=complex
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildMatrixInverseDiagonal(RCP<Matrix> const &A)
  {
    const RCP<const Map> rowmap = A->getRowMap();
    //std::vector<SC> diag(A->getNodeNumRows());
    std::vector<SC> diag(rowmap->getNodeNumElements());
    Teuchos::ArrayView<const LO> cols;
    Teuchos::ArrayView<const SC> vals;
    //for (size_t i=0; i<A->getNodeNumRows(); ++i)
    LO rowmapLocalSize = (LO) rowmap->getNodeNumElements();
    for (LO i=0; i<rowmapLocalSize; ++i) {
      A->getLocalRowView(i,cols,vals);
      for (LO j=0; j<cols.size(); ++j) {
        //TODO this will break down if diagonal entry is not present
        if (cols[j] == i) {
          diag[i] = 1 / vals[j];
          break;
        }
      }
    }

    RCP< Matrix > D = rcp( new CrsMatrixWrap(rowmap, 1) );
    std::vector<GO> diagInd(1);
    Teuchos::ArrayView<GO> iv(&diagInd[0],1);
    //for (size_t i=0; i< A->getNodeNumRows(); ++i)

    for (size_t i=0; i< rowmap->getNodeNumElements(); ++i) {
      Teuchos::ArrayView<SC> av(&diag[i],1);
      diagInd[0] = rowmap->getGlobalElement(i);
      D->insertGlobalValues(rowmap->getGlobalElement(i),iv,av); //TODO is this expensive?
    }
    D->fillComplete();

    return D;

  } //BuildMatrixInverseDiagonal()
#endif

  //   typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
  Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ResidualNorm(Matrix const &Op, MultiVector const &X, MultiVector const &RHS) {
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
    const size_t numVecs = X.getNumVectors();

    RCP<MultiVector> RES = Residual(Op, X, RHS);
    Teuchos::Array<Magnitude> norms(numVecs);
    RES->norm2(norms);

    return norms;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Residual(Matrix const &Op, MultiVector const &X, MultiVector const &RHS) {
    TEUCHOS_TEST_FOR_EXCEPTION(X.getNumVectors() != RHS.getNumVectors(), Exceptions::RuntimeError, "Number of solution vectors != number of right-hand sides")
    const size_t numVecs = X.getNumVectors();

    SC one    = Teuchos::ScalarTraits<Scalar>::one(), negone = -one, zero = Teuchos::ScalarTraits<Scalar>::zero();

    RCP<MultiVector> RES = MultiVectorFactory::Build(Op.getRangeMap(), numVecs);
    Op.apply(X, *RES, Teuchos::NO_TRANS, one, zero);
    RES->update(one, RHS, negone);

    return RES;
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write(std::string const & fileName, Matrix const & Op) {
    CrsMatrixWrap const & crsOp = dynamic_cast<CrsMatrixWrap const &>(Op);
    RCP<const CrsMatrix> tmp_CrsMtx = crsOp.getCrsMatrix();
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    const RCP<const EpetraCrsMatrix> &tmp_ECrsMtx = rcp_dynamic_cast<const EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx != Teuchos::null) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      RCP<const Epetra_CrsMatrix> A = tmp_ECrsMtx->getEpetra_CrsMatrix();
      int rv = EpetraExt::RowMatrixToMatrixMarketFile(fileName.c_str(), *A);
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "EpetraExt::RowMatrixToMatrixMarketFile return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }
#else
      throw(Exceptions::RuntimeError("Compiled without EpetraExt"));
#endif
      return;
    }
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraCrsMatrix> &tmp_TCrsMtx = rcp_dynamic_cast<const TpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_TCrsMtx != Teuchos::null) {
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A = tmp_TCrsMtx->getTpetra_CrsMatrix();
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >::writeSparseFile(fileName,A);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw(Exceptions::BadCast("Could not cast to EpetraCrsMatrix or TpetraCrsMatrix in matrix writing"));

  } //Write

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Read(std::string const & fileName, Xpetra::UnderlyingLib lib, RCP<const Teuchos::Comm<int> > const &comm)
  {

    if (lib == Xpetra::UseEpetra) {

#     if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
      
      Epetra_CrsMatrix *A;
      const RCP<const Epetra_Comm> epcomm = Xpetra::toEpetra(comm);
      int rv = EpetraExt::MatrixMarketFileToCrsMatrix( fileName.c_str(), *epcomm, A);
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "EpetraExt::MatlabFileToCrsMatrix return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }
      RCP<Epetra_CrsMatrix> tmpA = rcp(A);
      RCP<Matrix> rcpA = Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(tmpA);
      return rcpA;
#     else
      throw(Exceptions::RuntimeError("MueLu has not been compiled with Epetra and EpetraExt support."));
#     endif

    } else if (lib == Xpetra::UseTpetra) {

#ifdef HAVE_MUELU_TPETRA
      typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> sparse_matrix_type;
      typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;

      RCP<Node> node = Xpetra::DefaultPlatform::getDefaultPlatform().getNode();
      bool callFillComplete = true;
      RCP<sparse_matrix_type> tpA = reader_type::readSparseFile(fileName, comm, node, callFillComplete);
      if (tpA.is_null())
        throw(Exceptions::RuntimeError("The Tpetra::CrsMatrix returned from readSparseFile() is null."));
      RCP<TpetraCrsMatrix> tmpA1 = rcp(new TpetraCrsMatrix(tpA) );
      RCP<CrsMatrix> tmpA2 = rcp_implicit_cast<CrsMatrix>(tmpA1);
      RCP<Matrix> returnA = rcp( new CrsMatrixWrap(tmpA2) );
      return returnA;

#     else
      throw(Exceptions::RuntimeError("MueLu has not been compiled with Tpetra support."));
#     endif

    } else {
        throw(Exceptions::RuntimeError("Utils::Read : you must specify Xpetra::UseEpetra or Xpetra::UseTpetra."));
    }

    return Teuchos::null;
  } //Read()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Write(std::string const & fileName, const MultiVector& x) {
    RCP<const MultiVector> tmp_Vec = rcpFromRef(x);
#ifdef HAVE_MUELU_EPETRAEXT
    const RCP<const EpetraMultiVector> &tmp_EVec = rcp_dynamic_cast<const EpetraMultiVector>(tmp_Vec);
    if (tmp_EVec != Teuchos::null) {
#ifdef HAVE_MUELU_EPETRAEXT
      int rv = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(), *(tmp_EVec->getEpetra_MultiVector()));
      if (rv != 0) {
        std::ostringstream buf;
        buf << rv;
        std::string msg = "EpetraExt::RowMatrixToMatrixMarketFile return value of " + buf.str();
        throw(Exceptions::RuntimeError(msg));
      }
#else
      throw(Exceptions::RuntimeError("Compiled without EpetraExt"));
#endif
      return;
    }
#endif // HAVE_MUELU_EPETRAEXT

#ifdef HAVE_MUELU_TPETRA
    const RCP<const TpetraMultiVector> &tmp_TVec = rcp_dynamic_cast<const TpetraMultiVector>(tmp_Vec);
    if (tmp_TVec != Teuchos::null) {
      RCP<const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > TVec = tmp_TVec->getTpetra_MultiVector();
      Tpetra::MatrixMarket::Writer<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::writeDenseFile(fileName, TVec);
      return;
    }
#endif // HAVE_MUELU_TPETRA

    throw(Exceptions::BadCast("Could not cast to EpetraMultiVector or TpetraMultiVector in matrix writing"));

  } //Write

#include <unistd.h>


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PauseForDebugger()
  {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    int mypid = comm->getRank();

    for (int i = 0; i <comm->getSize(); i++) {
      if (i == mypid ) {
        char buf[80];
        char hostname[80];
        gethostname(hostname, sizeof(hostname));
        int pid = getpid();
        sprintf(buf, "Host: %s\tMPI rank: %d,\tPID: %d\n\tattach %d\n\tcontinue\n",
                hostname, mypid, pid, pid);
        printf("%s\n",buf);
        fflush(stdout);
        sleep(1);
      }
    }

    if (mypid == 0) {
      printf( "** Enter a character to continue > "); fflush(stdout);
      char go = ' ';
      scanf("%c",&go);
    }
    comm->barrier();
  } //PauseForDebugger

#ifdef HAVE_MUELU_TPETRA

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::simple_Transpose(RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > const &A)
  {
    LocalOrdinal N=A->getNodeNumRows();
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > AT=rcp(new Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>(A->getDomainMap(),0));
    const RCP<const Tpetra::Map<LO,GO,NO> > & rowMap=A->getRowMap();
    const RCP<const Tpetra::Map<LO,GO,NO> > & colMap=A->getColMap();

    for(LO i=0;i<N;i++){
      GO grid= rowMap->getGlobalElement(i);
      Teuchos::ArrayView<const LO> indices;
      Teuchos::ArrayView<const SC> vals;
      A->getLocalRowView(i,indices,vals);
      for(LO j=0;j<indices.size();j++){
        GO gcid=colMap->getGlobalElement(indices[j]);
        AT->insertGlobalValues(gcid,Teuchos::tuple(grid),Teuchos::tuple(vals[j]));
      }
    }
    //AT->fillComplete(A->getRangeMap(),A->getDomainMap());

    return AT;
  } //simple_Transpose
#endif // HAVE_MUELU_TPETRA

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Epetra_CrsMatrix> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::simple_EpetraTranspose(RCP<const Epetra_CrsMatrix> const &A)
  {
    int N=A->NumMyRows();
    RCP<Epetra_CrsMatrix> AT=rcp(new Epetra_CrsMatrix(Copy,A->DomainMap(),0));
    const Epetra_Map& rowMap=A->RowMap();
    const Epetra_Map& colMap=A->ColMap();

    for(int i=0;i<N;++i){
      int grid= rowMap.GID(i);
      int *indices,nnz;
      double *vals;
      A->ExtractMyRowView(i,nnz,vals,indices);
      for(int j=0;j<nnz;++j){
        int gcid=colMap.GID(indices[j]);

        //if(AT->RowMap().MyGID(gcid)) // no communication between procs!
        //{
        int rv = AT->InsertGlobalValues(gcid,1,vals+j,&grid);
        if (rv < 0) { // throw on errors, not on warnings...
          std::ostringstream buf;
          buf << rv;
          std::string msg = "Utils::simple_EpetraTranspose: Epetra_CrsMatrix::InsertGlobalValues() returned value of " + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }
        //}
      }
    }

    return AT;

  } //simple_Transpose
#endif


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PowerMethod(Matrix const &A, bool scaleByDiag,
                                                                                    LO niters, Magnitude tolerance, bool verbose, unsigned int seed)
  {
    if ( !(A.getRangeMap()->isSameAs(*(A.getDomainMap()))) ) {
      throw(Exceptions::Incompatible("Utils::PowerMethod: operator must have domain and range maps that are equivalent."));
    }
    // create three vectors, fill z with random numbers
    RCP<Vector> q = VectorFactory::Build(A.getDomainMap());
    RCP<Vector> qinit = VectorFactory::Build(A.getDomainMap());
    RCP<Vector> r = VectorFactory::Build(A.getRangeMap());
    RCP<Vector> z = VectorFactory::Build(A.getRangeMap());
    z->setSeed(seed);  // seed random number generator
    z->randomize(true);// use Xpetra implementation: -> same results for Epetra and Tpetra

    Teuchos::Array<Magnitude> norms(1);


    const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

    Scalar lambda=zero;
    Magnitude residual = 0.0;
    // power iteration
    RCP<Vector> diagVec,oneOverDiagonal;
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
      if ( iter % 100 == 0 || iter + 1 == niters ) {
        r->update(1.0, *z, -lambda, *q, zero);         // Compute A*q - lambda*q
        r->norm2(norms);
        residual = Teuchos::ScalarTraits<Scalar>::magnitude(norms[0] / lambda);
        if (verbose) {
          std::cout << "Iter = " << iter
                    << "  Lambda = " << lambda
                    << "  Residual of A*q - lambda*q = " << residual
                    << std::endl;
        }
      }
      if (residual < tolerance) {
        break;
      }
    }
    return lambda;
  } //PowerMethod


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix(RCP<Matrix> &Op, Teuchos::ArrayRCP<const SC> scalingVector, bool doInverse,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {

    Teuchos::ArrayRCP<SC> sv(scalingVector.size());
    if (doInverse) {
      for (int i=0; i<scalingVector.size(); ++i)
        sv[i] = 1.0 / scalingVector[i];
    } else {
      for (int i=0; i<scalingVector.size(); ++i)
        sv[i] = scalingVector[i];
    }

    switch (Op->getRowMap()->lib()) {
      case Xpetra::UseTpetra:
        MyOldScaleMatrix_Tpetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;
      case Xpetra::UseEpetra:
        Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Epetra(Op, sv, doFillComplete, doOptimizeStorage);
        break;
      default:
        throw(Exceptions::RuntimeError("Only Epetra and Tpetra matrices can be scaled."));
        break;
    } //switch
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Tpetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
#ifdef HAVE_MUELU_TPETRA
    RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpOp;
    try {
      tpOp = Op2NonConstTpetraCrs(Op);
    }
    catch(...) {
      throw(Exceptions::RuntimeError("Only Tpetra::CrsMatrix types can be scaled (Err.1)"));
    }

      const RCP<const Tpetra::Map<LO,GO,NO> > rowMap = tpOp->getRowMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > domainMap = tpOp->getDomainMap();
      const RCP<const Tpetra::Map<LO,GO,NO> > rangeMap = tpOp->getRangeMap();
      size_t maxRowSize = tpOp->getNodeMaxNumRowEntries();
      if (maxRowSize==(size_t)-1) //hasn't been determined yet
        maxRowSize=20;
      std::vector<SC> scaledVals(maxRowSize);
      if (tpOp->isFillComplete()) {
        tpOp->resumeFill();
      }

      if (Op->isLocallyIndexed() == true) {
        Teuchos::ArrayView<const LO> cols;
        Teuchos::ArrayView<const SC> vals;
        for (size_t i=0; i<rowMap->getNodeNumElements(); ++i) {
          tpOp->getLocalRowView(i,cols,vals);
          size_t nnz = tpOp->getNumEntriesInLocalRow(i);
          if (nnz>maxRowSize) {
            maxRowSize=nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j=0; j<nnz; ++j) {
            scaledVals[j] = vals[j]*scalingVector[i];
          }
          if (nnz>0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0],nnz);
            tpOp->replaceLocalValues(i,cols,valview);
          }
        } //for (size_t i=0; ...
      } else {
        Teuchos::ArrayView<const GO> cols;
        Teuchos::ArrayView<const SC> vals;
        for (size_t i=0; i<rowMap->getNodeNumElements(); ++i) {
          GO gid = rowMap->getGlobalElement(i);
          tpOp->getGlobalRowView(gid,cols,vals);
          size_t nnz = tpOp->getNumEntriesInGlobalRow(gid);
          if (nnz>maxRowSize) {
            maxRowSize=nnz;
            scaledVals.resize(maxRowSize);
          }
          for (size_t j=0; j<nnz; ++j) {
            scaledVals[j] = vals[j]*scalingVector[i]; //FIXME i or gid?
          }
          if (nnz>0) {
            Teuchos::ArrayView<const SC> valview(&scaledVals[0],nnz);
            tpOp->replaceGlobalValues(gid,cols,valview);
          }
        } //for (size_t i=0; ...
      }

      if (doFillComplete) {
        if (domainMap == Teuchos::null || rangeMap == Teuchos::null)
          throw(Exceptions::RuntimeError("In Utils::Scaling: cannot fillComplete because the domain and/or range map hasn't been defined"));
        RCP<Teuchos::ParameterList> params = rcp(new Teuchos::ParameterList());
        params->set("Optimize Storage",doOptimizeStorage);
        Op->fillComplete(Op->getDomainMap(),Op->getRangeMap(),params);
      }
#else
      throw(Exceptions::RuntimeError("Matrix scaling is not possible because Tpetra has not been enabled."));
#endif
  } //MyOldScaleMatrix_Tpetra()

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ArrayRCP<double> Utils<Scalar,LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::CoalesceCoordinates(Teuchos::ArrayRCP<double> coord, LocalOrdinal blksize) {
    if (blksize == 1)
      return coord;

    ArrayRCP<double> coalesceCoord(coord.size()/blksize); //TODO: how to avoid automatic initialization of the vector? using arcp()?

    for(int i=0; i<coord.size(); ++i) {
#define myDEBUG
#ifdef myDEBUG //FIXME-> HAVE_MUELU_DEBUG
      for(int j=1; j < blksize; ++j) {
        TEUCHOS_TEST_FOR_EXCEPTION(coord[i*blksize + j] != coord[i*blksize], Exceptions::RuntimeError, "MueLu::ZoltanInterface: coalesceCoord problem");
      }
#endif
      coalesceCoord[i] = coalesceCoord[i*blksize];
    }

    //std::cout << coord << std::endl;
    //std::cout << coalesceCoord << std::endl;

    return coalesceCoord;
  } //CoalesceCoordinates

  typedef Kokkos::DefaultNode::DefaultNodeType KDNT;
  typedef Kokkos::DefaultKernels<void,int,Kokkos::DefaultNode::DefaultNodeType>::SparseOps KDKSO;


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Teuchos::FancyOStream> Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MakeFancy(std::ostream & os) {
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
  ArrayRCP<const bool>
  Utils<SC, LO, GO, NO, LMO>::DetectDirichletRows(Matrix const &A, typename Teuchos::ScalarTraits<SC>::magnitudeType const &tol)
  {
    const RCP<const Map> rowMap = A.getRowMap();
    ArrayRCP<bool> boundaryNodes(A.getNodeNumRows(),true);

    for(LO row=0; row < Teuchos::as<LO>(rowMap->getNodeNumElements()); ++row) {

      ArrayView<const LO> indices;
      ArrayView<const SC> vals;
      A.getLocalRowView(row, indices, vals);
      size_t nnz = A.getNumEntriesInLocalRow(row);
      if (nnz > 1) {
        for(size_t col=0; col<nnz; ++col) {
          if ( (indices[col] != row) && Teuchos::ScalarTraits<SC>::magnitude(vals[col]) > tol) {
            boundaryNodes[row] = false;
            break;
          }
        }
      }
    }
    return boundaryNodes;
  } //DetectDirichletRows

#ifdef HAVE_MUELU_EPETRA
//   template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
//   RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap(RCP<Epetra_CrsMatrix> &epAB) {
//     TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
//     return Teuchos::null;
//   }


//   template<>
//   inline RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > Convert_Epetra_CrsMatrix_ToXpetra_CrsMatrixWrap<double,int,int,KDNT,KDKSO > (RCP<Epetra_CrsMatrix> &epAB) {
//     RCP<Xpetra::EpetraCrsMatrix> tmpC1 = rcp(new Xpetra::EpetraCrsMatrix(epAB));
//     RCP<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> > tmpC2 = rcp_implicit_cast<Xpetra::CrsMatrix<double,int,int,KDNT,KDKSO> >(tmpC1);
//     RCP<Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO> > tmpC3 = rcp(new Xpetra::CrsMatrixWrap<double,int,int,KDNT,KDKSO>(tmpC2));
//     return tmpC3;
//   }
#endif

//   template<class T>
//   std::string toString(T const &what) {
//     std::ostringstream buf; buf << what;
//     return buf.str();
//   }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Transpose(RCP<Matrix> const &Op, bool const & optimizeTranspose)
  {
   Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer("YY Entire Transpose"));
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    std::string TorE = "epetra";
#else
    std::string TorE = "tpetra";
#endif

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
    RCP<Epetra_CrsMatrix> epetraOp;
    try {
      epetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2NonConstEpetraCrs(Op);
    }
    catch (...) {
      TorE = "tpetra";
    }
#endif

#ifdef HAVE_MUELU_TPETRA
    RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpetraOp;
    if (TorE=="tpetra") {
      try {
        tpetraOp = Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Op2TpetraCrs(Op);
      }
      catch (...) {
        throw(Exceptions::RuntimeError("Utils::Transpose: Can only transpose Crs matrices"));
      }
    } //if
#endif

    if (TorE == "tpetra") {
#ifdef HAVE_MUELU_TPETRA
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A;
      {
      Teuchos::TimeMonitor tmm(*Teuchos::TimeMonitor::getNewTimer("YY Tpetra Transpose Only"));
      Tpetra::RowMatrixTransposer<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> transposer(*tpetraOp); //more than meets the eye
      A = transposer.createTranspose(optimizeTranspose ? Tpetra::DoOptimizeStorage : Tpetra::DoNotOptimizeStorage); //couldn't have just used a bool...
      }

      //RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > A=Utils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::simple_Transpose(tpetraOp);
      RCP<TpetraCrsMatrix> AA = rcp(new TpetraCrsMatrix(A) );
      RCP<CrsMatrix> AAA = rcp_implicit_cast<CrsMatrix>(AA);
      RCP<Matrix> AAAA = rcp( new CrsMatrixWrap(AAA) );
      if (!AAAA->isFillComplete())
        AAAA->fillComplete(Op->getRangeMap(),Op->getDomainMap());
      return AAAA;
#else
      throw(Exceptions::RuntimeError("Tpetra"));
#endif
    }

    //epetra case
    std::cout << "Utilities::Transpose() not implemented for Epetra" << std::endl;
    return Teuchos::null;

  } //Transpose

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MyOldScaleMatrix_Epetra(RCP<Matrix> &Op, Teuchos::ArrayRCP<SC> const &scalingVector,
                               bool doFillComplete,
                               bool doOptimizeStorage)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MyOldScalematrix and Epetra cannot be used with Scalar != double, LocalOrdinal != int, GlobalOrdinal != int");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(RCP<Matrix> const &A, bool transposeA, SC alpha, RCP<Matrix> &B, SC beta)
  {
    if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
      throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
    }

    if (A->getRowMap()->lib() == Xpetra::UseEpetra) {
      throw(Exceptions::RuntimeError("You cannot use Epetra::MatrixMatrix::Add with Scalar!=double or Ordinal!=int"));
    } else if(A->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpA = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(A);
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpB = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2NonConstTpetraCrs(B);

      Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, beta);
#else
      throw(Exceptions::RuntimeError("MueLu must be compiled with Tpetra."));
#endif
    }

  } //Utils2::TwoMatrixAdd()


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Utils2<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::TwoMatrixAdd(RCP<Matrix> const &A, bool const &transposeA, SC const &alpha,
                           RCP<Matrix> const &B, bool const &transposeB, SC const &beta,
                           RCP<Matrix> &C, bool const &AHasFixedNnzPerRow)
  {
    if ( !(A->getRowMap()->isSameAs(*(B->getRowMap()))) ) {
      throw(Exceptions::Incompatible("TwoMatrixAdd: matrix row maps are not the same."));
    }
    if (C==Teuchos::null) {
      if (!A->isFillComplete() || !B->isFillComplete())
        TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,"Global statistics are not available for estimates.");

      size_t maxNzInA = A->getGlobalMaxNumRowEntries();
      size_t maxNzInB = B->getGlobalMaxNumRowEntries();
      size_t numLocalRows = A->getNodeNumRows();

      RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      fos->setOutputToRootOnly(0);

      if (maxNzInA == 1 || maxNzInB == 1 || AHasFixedNnzPerRow) {
        //first check if either A or B has at most 1 nonzero per row
        //the case of both having at most 1 nz per row is handled by the ``else''
        Teuchos::ArrayRCP<size_t> exactNnzPerRow(numLocalRows);
        if ( (maxNzInA == 1 && maxNzInB > 1) || AHasFixedNnzPerRow) {
          for (size_t i=0; i<numLocalRows; ++i)
            exactNnzPerRow[i] = B->getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInA;
        } else {
          for (size_t i=0; i<numLocalRows; ++i)
            exactNnzPerRow[i] = A->getNumEntriesInLocalRow(Teuchos::as<LO>(i)) + maxNzInB;
        }
        *fos << "Utils::TwoMatrixAdd : special case detected (one matrix has a fixed nnz per row)"
             << ", using static profiling" << std::endl;
        C = rcp( new CrsMatrixWrap(A->getRowMap(), exactNnzPerRow, Xpetra::StaticProfile) );
      }
      else {
        //general case
        double nnzPerRowInA = Teuchos::as<double>(A->getGlobalNumEntries()) / A->getGlobalNumRows();
        double nnzPerRowInB = Teuchos::as<double>(B->getGlobalNumEntries()) / B->getGlobalNumRows();
        LO nnzToAllocate = Teuchos::as<LO>( (nnzPerRowInA + nnzPerRowInB) * 1.5) + Teuchos::as<LO>(1);

        LO maxPossible = A->getGlobalMaxNumRowEntries() + B->getGlobalMaxNumRowEntries();
        //Use static profiling (more efficient) if the estimate is at least as big as the max
        //possible nnz's in any single row of the result.
        Xpetra::ProfileType pft = (maxPossible) > nnzToAllocate ? Xpetra::DynamicProfile : Xpetra::StaticProfile;

        *fos << "nnzPerRowInA = " << nnzPerRowInA << ", nnzPerRowInB = " << nnzPerRowInB << std::endl;
        *fos << "Utils::TwoMatrixAdd : space allocated per row = " << nnzToAllocate
             << ", max possible nnz per row in sum = " << maxPossible
             << ", using " << (pft == Xpetra::DynamicProfile ? "dynamic" : "static" ) << " profiling"
             << std::endl;
        C = rcp( new CrsMatrixWrap(A->getRowMap(), nnzToAllocate, pft) );
      }
      if (transposeB)
        *fos << "Utils::TwoMatrixAdd : ** WARNING ** estimate could be badly wrong because second summand is transposed" << std::endl;
    }

    if (C->getRowMap()->lib() == Xpetra::UseEpetra) {
      throw(Exceptions::RuntimeError("You cannot use Epetra::MatrixMatrix::Add with Scalar!=double or Ordinal!=int"));
    } else if(C->getRowMap()->lib() == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpA = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(A);
      RCP<const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > tpB = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2TpetraCrs(B);
      RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> >       tpC = Utils<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps>::Op2NonConstTpetraCrs(C);

      Tpetra::MatrixMatrix::Add(*tpA, transposeA, alpha, *tpB, transposeB, beta, tpC);
#else
      throw(Exceptions::RuntimeError("MueLu must be compile with Tpetra."));
#endif
    }

    ///////////////////////// EXPERIMENTAL
    if(A->IsView("stridedMaps")) C->CreateView("stridedMaps", A);
    if(B->IsView("stridedMaps")) C->CreateView("stridedMaps", B);
    ///////////////////////// EXPERIMENTAL

  } //Utils2::TwoMatrixAdd()

} //namespace MueLu

#define MUELU_UTILITIES_SHORT
#endif // MUELU_UTILITIES_DEF_HPP

//  LocalWords:  LocalOrdinal
