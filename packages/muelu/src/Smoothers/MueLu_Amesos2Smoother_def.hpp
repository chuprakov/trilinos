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
#ifndef MUELU_AMESOS2SMOOTHER_DEF_HPP
#define MUELU_AMESOS2SMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"
#if defined (HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_AMESOS2)
#include <Xpetra_Matrix.hpp>

#include <Amesos2_config.h>
#include <Amesos2.hpp>

#include "MueLu_Amesos2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  // Metaprogramming: CheckNodeType<T>::val would be 1 only for the default node type
  template<typename T> struct CheckNodeType                                              { enum { val = 0 }; };
  template <>          struct CheckNodeType<KokkosClassic::DefaultNode::DefaultNodeType> { enum { val = 1 }; };

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Amesos2Smoother(const std::string& type, const Teuchos::ParameterList& paramList)
    : type_(type), paramList_(paramList)
  {
    if (!CheckNodeType<Node>::val)
      throw Exceptions::NotImplemented("Amesos2Smoother works only with serial node type");

    // Set default solver type
    // TODO: It would be great is Amesos2 provides directly this kind of logic for us
    if(type_ == "") {
#if defined(HAVE_AMESOS2_SUPERLUDIST)
      type_ = "Superludist";
#elif defined(HAVE_AMESOS2_KLU2)
      type_ = "Klu";
#elif defined(HAVE_AMESOS2_SUPERLU)
      type_ = "Superlu";
#else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Amesos2Smoother::Amesos2Smoother(): Amesos2 have been compiled without SuperLU_DIST, SuperLU or Klu. "
                                 "By default, MueLu tries to use one of these libraries. Amesos2 must be compiled with one of these solvers or a valid Amesos2 solver have to be specified explicitly.");
#endif
    } // if(type_ == "")

    //TMP: Amesos2 KLU never available but most MueLu tests are using KLU by default
    // (ex: examples driven by ML parameter lists)
    // -> temporarily fallback to SUPERLU
    // Remove this when KLU becomes available.
#if defined(HAVE_AMESOS2_SUPERLU)
    if (type_ == "Klu" && Amesos2::query(type_) == false) {
      type_ = "Superlu";
      this->GetOStream(Warnings0, 0) << "Warning: MueLu::Amesos2Smoother: KLU2 not available. Using SuperLu instead" << std::endl;
    }
#endif // HAVE_AMESOS2_SUPERLU
    // END OF TMP

    // Check the validity of the solver type parameter
    TEUCHOS_TEST_FOR_EXCEPTION(Amesos2::query(type_) == false, Exceptions::RuntimeError, "MueLu::Amesos2Smoother::Amesos2Smoother(): The Amesos2 library reported that the solver '" << type_ << "' is not available. "
                               "Amesos2 have been compiled without the support of this solver or the solver name is misspelled.");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~Amesos2Smoother() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level& currentLevel) const {
    this->Input(currentLevel, "A");
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);
    if (SmootherPrototype::IsSetup() == true)
        this->GetOStream(Warnings0, 0) << "Warning: MueLu::Amesos2Smoother::Setup(): Setup() has already been called" << std::endl;

    RCP<Matrix> A = Factory::Get< RCP<Matrix> >(currentLevel, "A");

    if (CheckNodeType<Node>::val) {
      // template metaprogramming: if node type is inapplicable, it should not instantiate
      RCP<Tpetra_CrsMatrix> tA = Utils::Op2NonConstTpetraCrs(A);

      prec_ = Amesos2::create<Tpetra_CrsMatrix,Tpetra_MultiVector>(type_, tA);
    }

    TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "Could not create Amesos2 solver");

    //TODO      prec_->setParameters(paramList_);
    //TODO
    // int rv = prec_->numericFactorization();
    //       if (rv != 0) {
    //         std::ostringstream buf;
    //         buf << rv;
    //         std::string msg = "Amesos2_BaseSolver::NumericFactorization return value of " + buf.str(); //TODO: BaseSolver or ... ?
    //         throw(Exceptions::RuntimeError(msg));
    //       }

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector& X, const MultiVector& B, const bool& InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::Amesos2Smoother::Apply(): Setup() has not been called");

    if (CheckNodeType<Node>::val) {
      // template metaprogramming: if node type is inapplicable, it should not instantiate
      RCP<Tpetra_MultiVector> tX = Utils::MV2NonConstTpetraMV2(X);
      MultiVector& BNonC = const_cast<MultiVector&>(B);
      RCP<Tpetra_MultiVector> tB = Utils::MV2NonConstTpetraMV2(BNonC);

      prec_->setX(tX);
      prec_->setB(tB);

      prec_->solve();

      prec_->setX(Teuchos::null);
      prec_->setB(Teuchos::null);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp(new Amesos2Smoother(*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    if (CheckNodeType<Node>::val) {
      // template metaprogramming: if node type is inapplicable, it should not instantiate
      if (SmootherPrototype::IsSetup() == true) {
        out << prec_->description();
      } else {
        out << SmootherPrototype::description();
        out << "{type = " << type_ << "}";
      }
    }
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Amesos2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (CheckNodeType<Node>::val) {
      // template metaprogramming: if node type is inapplicable, it should not instantiate
      if (verbLevel & Parameters0)
        out0 << "Prec. type: " << type_ << std::endl;

      if (verbLevel & Parameters1) {
        out0 << "Parameter list: " << std::endl;
        Teuchos::OSTab tab2(out);
        out << paramList_;
      }

      if ((verbLevel & External) && (prec_ != Teuchos::null)) {
        Teuchos::OSTab tab2(out);
        out << *prec_ << std::endl;
      }

      if (verbLevel & Debug)
        out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
            << "-" << std::endl
            << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

} // namespace MueLu

#endif // HAVE_MUELU_TPETRA && HAVE_MUELU_AMESOS2
#endif // MUELU_AMESOS2SMOOTHER_DEF_HPP
