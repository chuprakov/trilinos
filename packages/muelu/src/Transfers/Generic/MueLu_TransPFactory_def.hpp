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
#ifndef MUELU_TRANSPFACTORY_DEF_HPP
#define MUELU_TRANSPFACTORY_DEF_HPP

// disable clang warnings
#ifdef __clang__
#pragma clang system_header
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

#include <Xpetra_Matrix.hpp>

#include "MueLu_TransPFactory_decl.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_DisableMultipleCallCheck.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<const ParameterList> TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetValidParameterList(const ParameterList& paramList) const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
    validParamList->set< RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory of the matrix P");
    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level & fineLevel, Level & coarseLevel) const {
    Input(coarseLevel, "P");
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level & fineLevel, Level & coarseLevel) const {
    FactoryMonitor m(*this, "Transpose P", coarseLevel);

    // PFact might be called twice: once for P, once for R. Disabling the debug check.
    RCP<const FactoryBase> PFact = GetFactory("P");
    if (PFact == Teuchos::null) { PFact = coarseLevel.GetFactoryManager()->GetFactory("P"); }
    MueLu::DisableMultipleCallCheck check(rcp_dynamic_cast<const TwoLevelFactoryBase>(PFact));

    //
    //
    //

    RCP<Matrix> P = Get< RCP<Matrix> >(coarseLevel, "P");

    //doesn't work -- bug in EpetraExt?
    // Teuchos::ParameterList matrixList;
    //RCP<Matrix> I = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>("Identity", P->getRangeMap(), matrixList);
    //      RCP<CrsMatrixWrap> I = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsMatrixWrap>("Identity", P->getDomainMap(), matrixList);
    //RCP<Matrix> R = Utils::TwoMatrixMultiply(P, true, I, false); //doesn't work -- bug in EpetraExt?
    //      RCP<Matrix> R = Utils::TwoMatrixMultiply(I, false, P, true);

    RCP<Matrix> R = Utils2::Transpose(*P, true);

    Set(coarseLevel, "R", R);

    ///////////////////////// EXPERIMENTAL
    if(P->IsView("stridedMaps")) R->CreateView("stridedMaps", P, true);
    ///////////////////////// EXPERIMENTAL

  } //Build

} //namespace MueLu

#endif // MUELU_TRANSPFACTORY_DEF_HPP
