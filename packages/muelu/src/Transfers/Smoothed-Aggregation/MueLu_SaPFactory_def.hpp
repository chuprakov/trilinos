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
#ifndef MUELU_SAPFACTORY_DEF_HPP
#define MUELU_SAPFACTORY_DEF_HPP

#include <Xpetra_Matrix.hpp>

#include "MueLu_SaPFactory_decl.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SaPFactory(RCP<const FactoryBase> InitialPFact, RCP<const FactoryBase> AFact)
    : initialPFact_(InitialPFact), AFact_(AFact),
      dampingFactor_(4./3), diagonalView_("current") {
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~SaPFactory() {}
  
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDampingFactor(Scalar dampingFactor) {
    dampingFactor_ = dampingFactor;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDiagonalView(std::string const& diagView) {
    diagonalView_ = diagView;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Scalar SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDampingFactor() {
    return dampingFactor_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDiagonalView() {
    return diagonalView_;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
    fineLevel.DeclareInput("A",AFact_.get(),this);

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = initialPFact_;
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }
    coarseLevel.DeclareInput("P",initialPFact.get(),this);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel,coarseLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void SaPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level &fineLevel, Level &coarseLevel) const {
    FactoryMonitor m(*this, "Prolongator smoothing", coarseLevel);

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    // Get default tentative prolongator factory
    // Getting it that way ensure that the same factory instance will be used for both SaPFactory and NullspaceFactory.
    // -- Warning: Do not use directly initialPFact_. Use initialPFact instead everywhere!
    RCP<const FactoryBase> initialPFact = initialPFact_;
    if (initialPFact == Teuchos::null) { initialPFact = coarseLevel.GetFactoryManager()->GetFactory("Ptent"); }

    // Level Get
    RCP<Matrix> A     = fineLevel.  Get< RCP<Matrix> >("A", AFact_.get());
    RCP<Matrix> Ptent = coarseLevel.Get< RCP<Matrix> >("P", initialPFact.get());

    if(restrictionMode_) {
      SubFactoryMonitor m2(*this, "Transpose A", coarseLevel);
      A = Utils2::Transpose(A,true); // build transpose of A explicitely
    }

    //Build final prolongator
    RCP<Matrix> finalP; // output

    //FIXME Xpetra::Matrix should calculate/stash max eigenvalue
    //FIXME SC lambdaMax = A->GetDinvALambda();

    if (dampingFactor_ != Teuchos::ScalarTraits<Scalar>::zero()) {

      //Teuchos::ParameterList matrixList;
      //RCP<Matrix> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map,CrsMatrixWrap>("Identity",fineLevel.Get< RCP<Matrix> >("A")->getRowMap(),matrixList);
      //RCP<Matrix> newPtent = Utils::TwoMatrixMultiply(I,false,Ptent,false);
      //Ptent = newPtent; //I tried a checkout of the original Ptent, and it seems to be gone now (which is good)

      RCP<Matrix> AP;
      {
        SubFactoryMonitor m2(*this, "MxM: A x Ptentative", coarseLevel);
        //JJH -- If I switch doFillComplete to false, the resulting matrix seems weird when printed with describe.
        //JJH -- The final prolongator is wrong, to boot.  So right now, I fillComplete AP, but avoid fillComplete
        //JJH -- in the scaling.  Long story short, we're doing 2 fillCompletes, where ideally we'd do just one.
        bool doFillComplete=true;
        bool optimizeStorage=false;
        AP = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);
      }

      {
        SubFactoryMonitor m2(*this, "Scaling (A x Ptentative) by D^{-1}", coarseLevel);
        bool doFillComplete=false;
        bool optimizeStorage=false;
        Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(*A);
        Utils::MyOldScaleMatrix(AP,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag
      }

      Scalar lambdaMax;
      {
        SubFactoryMonitor m2(*this, "Eigenvalue estimate", coarseLevel);
        Magnitude stopTol = 1e-4;
        lambdaMax = Utils::PowerMethod(*A, true, (LO) 10, stopTol);
        //Scalar lambdaMax = Utils::PowerMethod(*A, true, (LO) 50,(Scalar)1e-7, true);
        GetOStream(Statistics1, 0) << "Damping factor = " << dampingFactor_/lambdaMax << " (" << dampingFactor_ << " / " << lambdaMax << ")" << std::endl;
      }

      {
        SubFactoryMonitor m2(*this, "M+M: P = (Ptentative) + (D^{-1} x A x Ptentative)", coarseLevel);
        
        bool doTranspose=false; 
        if (AP->isFillComplete())
          Utils2::TwoMatrixAdd(Ptent,doTranspose,Teuchos::ScalarTraits<Scalar>::one(),AP,doTranspose,-dampingFactor_/lambdaMax,finalP);
        else {
          Utils2::TwoMatrixAdd(Ptent,doTranspose,Teuchos::ScalarTraits<Scalar>::one(),AP,-dampingFactor_/lambdaMax);
          finalP = AP;
        }
      }

      {
        SubFactoryMonitor m2(*this, "FillComplete() of P", coarseLevel);
        finalP->fillComplete( Ptent->getDomainMap(), Ptent->getRangeMap() );
      }
      
    } else {
      finalP = Ptent;
    }

    // Level Set
    if(!restrictionMode_)
      {
        // prolongation factory is in prolongation mode
        coarseLevel.Set("P", finalP, this);
	
        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) finalP->CreateView("stridedMaps", Ptent);
        ///////////////////////// EXPERIMENTAL
      }
    else
      {
        // prolongation factory is in restriction mode
        RCP<Matrix> R = Utils2::Transpose(finalP,true); // use Utils2 -> specialization for double
        coarseLevel.Set("R", R, this);
	
        ///////////////////////// EXPERIMENTAL
        if(Ptent->IsView("stridedMaps")) R->CreateView("stridedMaps", Ptent, true);
        ///////////////////////// EXPERIMENTAL
      }

  } //Build()

} //namespace MueLu

#endif // MUELU_SAPFACTORY_DEF_HPP
