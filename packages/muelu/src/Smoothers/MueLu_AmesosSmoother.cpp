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
#include <algorithm>

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_AMESOS)

#include <Epetra_LinearProblem.h>

#include <Amesos_config.h>
#include <Amesos.h>
#include <Amesos_BaseSolver.h>

#include "MueLu_AmesosSmoother.hpp"

#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  AmesosSmoother::AmesosSmoother(const std::string& type, const Teuchos::ParameterList& paramList)
    : type_(type), paramList_(paramList) {
    if (!type_.empty()) {
      // Transform string to "Abcde" notation
      std::transform(type_.begin(),   type_.end(),   type_.begin(), ::tolower);
      std::transform(type_.begin(), ++type_.begin(), type_.begin(), ::toupper);
    }
    if (type_ == "Amesos_klu")      type_ = "Klu";
    if (type_ == "Klu2")            type_ = "Klu";
    if (type_ == "Amesos_umfpack")  type_ = "Umfpack";
    if (type_ == "Superlu_dist")    type_ = "Superludist";

    // Try to come up with something availble
    // Order corresponds to our preference
    // TODO: It would be great is Amesos provides directly this kind of logic for us
    std::string oldtype = type_;
    if (type_ == "" || Amesos().Query(type_) == false) {
#if defined(HAVE_AMESOS_SUPERLU)
      type_ = "Superlu";
#elif defined(HAVE_AMESOS_KLU)
      type_ = "Klu";
#elif defined(HAVE_AMESOS_SUPERLUDIST)
      type_ = "Superludist";
#elif defined(HAVE_AMESOS_UMFPACK)
      type_ = "Umfpack";
#else
      throw Exceptions::RuntimeError("Amesos have been compiled without SuperLU_DIST, SuperLU, Umfpack or Klu. By default, MueLu tries"
                                     "to use one of these libraries. Amesos2 must be compiled with one of these solvers or"
                                     "a valid Amesos solver have to be specified explicitly.");
#endif
      if (oldtype != "")
        this->GetOStream(Warnings0, 0) << "Warning: MueLu::AmesosSmoother: \"" << oldtype << "\" is not available. Using \"" << type_ << "\" instead" << std::endl;
      else
        this->GetOStream(Warnings0, 0) << "MueLu::AmesosSmoother: using \"" << type_ << "\"" << std::endl;
    }
  }

  AmesosSmoother::~AmesosSmoother() {}

  void AmesosSmoother::DeclareInput(Level &currentLevel) const {
    Input(currentLevel, "A");
  }

  void AmesosSmoother::Setup(Level& currentLevel) {
    FactoryMonitor m(*this, "Setup Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
      GetOStream(Warnings0, 0) << "Warning: MueLu::AmesosSmoother::Setup(): Setup() has already been called" << std::endl;

    A_ = Get< RCP<Matrix> >(currentLevel, "A");

    RCP<Epetra_CrsMatrix> epA = Utils::Op2NonConstEpetraCrs(A_);
    linearProblem_ = rcp( new Epetra_LinearProblem() );
    linearProblem_->SetOperator(epA.get());

    Amesos factory;
    prec_ = rcp(factory.Create(type_, *linearProblem_));
    TEUCHOS_TEST_FOR_EXCEPTION(prec_ == Teuchos::null, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Solver '" + type_ + "' not supported by Amesos");

    // set Reindex flag, if A is distributed with non-contiguous maps
    // unfortunately there is no reindex for Amesos2, yet. So, this only works for Epetra based problems
    if (A_->getRowMap()->isDistributed() == true && A_->getRowMap()->isContiguous() == false)
      paramList_.set("Reindex", true);

    prec_->SetParameters(paramList_);

    int r = prec_->NumericFactorization();
    TEUCHOS_TEST_FOR_EXCEPTION(r != 0, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Setup(): Amesos solver returns value of " +
                               Teuchos::Utils::toString(r) + " during NumericFactorization()");

    SmootherPrototype::IsSetup(true);
  }

  void AmesosSmoother::Apply(MultiVector& X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::AmesosSmoother::Apply(): Setup() has not been called");

    Epetra_MultiVector &epX = Utils::MV2NonConstEpetraMV(X);
    Epetra_MultiVector const &epB = Utils::MV2EpetraMV(B);
    //Epetra_LinearProblem takes the right-hand side as a non-const pointer.
    //I think this const_cast is safe because Amesos won't modify the rhs.
    Epetra_MultiVector &nonconstB = const_cast<Epetra_MultiVector&>(epB);

    linearProblem_->SetLHS(&epX);
    linearProblem_->SetRHS(&nonconstB);

    prec_->Solve();

    // Don't keep pointers to our vectors in the Epetra_LinearProblem.
    linearProblem_->SetLHS(0);
    linearProblem_->SetRHS(0);
  }

  RCP<MueLu::SmootherPrototype<double,int,int> > AmesosSmoother::Copy() const {
    return rcp( new AmesosSmoother(*this) );
  }

  std::string AmesosSmoother::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  //using MueLu::Describable::describe; // overloading, not hiding
  void AmesosSmoother::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }

    if (verbLevel & Parameters1) {
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
    }

    if (verbLevel & External) {
      if (prec_ != Teuchos::null) { prec_->PrintStatus(); prec_->PrintTiming(); } //TODO: redirect output?
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<A_>: " << A_ << std::endl
           << "RCP<linearProblem__>: " << linearProblem_ << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

} // namespace MueLu

#endif // HAVE_MUELU_EPETRA && HAVE_MUELU_AMESOS
