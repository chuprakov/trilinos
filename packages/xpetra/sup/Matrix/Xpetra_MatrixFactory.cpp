// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
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

#ifndef XPETRA_MATRIXFACTORY_CPP
#define XPETRA_MATRIXFACTORY_CPP

#include "Xpetra_MatrixFactory.hpp"

namespace Xpetra {

  RCP<Xpetra::Matrix<double,int,int> > MatrixFactory2<double,int,int>::BuildCopy(const RCP<const Matrix> A) {
    RCP<const CrsMatrixWrap> oldOp = Teuchos::rcp_dynamic_cast<const CrsMatrixWrap>(A);
    if (oldOp == Teuchos::null)
      throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::CrsMatrixWrap failed");

    RCP<const CrsMatrix> oldCrsOp = oldOp->getCrsMatrix();

#ifdef HAVE_XPETRA_EPETRA
    RCP<const EpetraCrsMatrix> oldECrsOp = Teuchos::rcp_dynamic_cast<const EpetraCrsMatrix>(oldCrsOp);
    if (oldECrsOp != Teuchos::null) {
      // Underlying matrix is Epetra
      RCP<CrsMatrix>     newECrsOp(new EpetraCrsMatrix(*oldECrsOp));
      RCP<CrsMatrixWrap> newOp    (new CrsMatrixWrap  (newECrsOp));

      return newOp;
    }
#endif

#ifdef HAVE_XPETRA_TPETRA
    // Underlying matrix is Tpetra
    RCP<const TpetraCrsMatrix> oldTCrsOp = Teuchos::rcp_dynamic_cast<const TpetraCrsMatrix>(oldCrsOp);
    if (oldTCrsOp != Teuchos::null) {
      RCP<CrsMatrix>     newTCrsOp(new TpetraCrsMatrix(*oldTCrsOp));
      RCP<CrsMatrixWrap> newOp    (new CrsMatrixWrap(newTCrsOp));

      return newOp;
    }
#else
    throw Exceptions::BadCast("Cast from Xpetra::Matrix to Xpetra::EpetraCrsMatrix or Xpetra::TpetraCrsMatrix failed");
#endif

    return Teuchos::null;  // make compiler happy
  }

} // namespace Xpetra

#endif
