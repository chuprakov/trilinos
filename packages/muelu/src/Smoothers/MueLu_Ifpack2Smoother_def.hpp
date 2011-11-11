#ifndef MUELU_IFPACK2SMOOTHER_DEF_HPP
#define MUELU_IFPACK2SMOOTHER_DEF_HPP

#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_IFPACK2

#include "MueLu_Ifpack2Smoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Ifpack2Smoother(std::string const & type, Teuchos::ParameterList const & paramList, LO const &overlap, RCP<FactoryBase> AFact) //TODO: empty paramList valid for Ifpack??
    : type_(type), paramList_(paramList), overlap_(overlap), AFact_(AFact)
  { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~Ifpack2Smoother() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetParameters(Teuchos::ParameterList const & paramList) {
    paramList_ = paramList;

    if (SmootherPrototype::IsSetup()) {
      // It might be invalid to change parameters after the setup, but it depends entirely on Ifpack implementation.
      // TODO: I don't know if Ifpack returns an error code or exception or ignore parameters modification in this case...

      Teuchos::ParameterList nonConstParamList = paramList; // because Ifpack SetParameters() input argument is not const...
      prec_->setParameters(nonConstParamList);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  Teuchos::ParameterList const & Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetParameters() { return paramList_; }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    currentLevel.DeclareInput("A", AFact_.get());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    Monitor m(*this, "Setup Smoother");
    if (this->IsSetup() == true) this->GetOStream(Warnings0, 0) << "Warning: MueLu::Ifpack2Smoother::Setup(): Setup() has already been called";

    RCP<Operator> A = currentLevel.Get< RCP<Operator> >("A", AFact_.get());

    if (type_ == "CHEBYSHEV") {
      Scalar maxEigenValue = paramList_.get("chebyshev: max eigenvalue",(Scalar)-1.0);
      if (maxEigenValue == -1.0) {
        maxEigenValue = Utils::PowerMethod(*A,true,10,1e-4);
        paramList_.set("chebyshev: max eigenvalue", maxEigenValue);

        this->GetOStream(Statistics1, 0) << "chebyshev: max eigenvalue" << " = " << maxEigenValue << std::endl;
      }
    }

    RCP<const Tpetra::CrsMatrix<SC, LO, GO, NO, LMO> > tpA = Utils::Op2NonConstTpetraCrs(A);
    prec_ = Ifpack2::Factory::create(type_, tpA, overlap_);

    prec_->setParameters(paramList_);
    prec_->initialize();
    prec_->compute();

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector& X, MultiVector const &B, bool const &InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::IfpackSmoother::Apply(): Setup() has not been called");

    // Forward the InitialGuessIsZero option to Ifpack2
    //  TODO: It might be nice to switch back the internal
    //        "zero starting solution" option of the ifpack2 object prec_ to his
    //        initial value at the end but there is no way right now to get
    //        the current value of the "zero starting solution" in ifpack2.
    //        It's not really an issue, as prec_  can only be used by this method.
    Teuchos::ParameterList  paramList;
    if (type_ == "CHEBYSHEV") {
      paramList.set("chebyshev: zero starting solution", InitialGuessIsZero);
    }
    else if (type_ == "RELAXATION") {
      paramList.set("relaxation: zero starting solution", InitialGuessIsZero);
    }
    else if (type_ == "ILUT") {
      if (InitialGuessIsZero == false) {
        if (this->IsPrint(Warnings0, 0)) {
          static int warning_only_once=0;
          if ((warning_only_once++) == 0)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::Ifpack2Smoother::Apply(): ILUT has no provision for a nonzero initial guess." << std::endl;
          // TODO: ILUT using correction equation should be implemented in ifpack2 directly
          //       I think that an option named "zero starting solution"
          //       is also appropriate for ILUT
        }
      }
    } else {
      // TODO: When https://software.sandia.gov/bugzilla/show_bug.cgi?id=5283#c2 is done
      // we should remove the if/else/elseif and just test if this
      // option is supported by current ifpack2 preconditioner
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError,"IfpackSmoother::Apply(): Ifpack preconditioner '"+type_+"' not supported");
    }
    prec_->setParameters(paramList);

    // Apply
    Tpetra::MultiVector<SC,LO,GO,NO> &tpX = Utils::MV2NonConstTpetraMV(X);
    Tpetra::MultiVector<SC,LO,GO,NO> const &tpB = Utils::MV2TpetraMV(B);
    prec_->apply(tpB,tpX);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp(new Ifpack2Smoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }
    
  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << std::endl;
    }
      
    if (verbLevel & Parameters1) { 
      out0 << "Parameter list: " << std::endl; { Teuchos::OSTab tab2(out); out << paramList_; }
      out0 << "Overlap: "        << overlap_ << std::endl;
    }
      
    if (verbLevel & External) {
      if (prec_ != Teuchos::null) { Teuchos::OSTab tab2(out); out << *prec_ << std::endl; }
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl
           << "-" << std::endl
           << "RCP<prec_>: " << prec_ << std::endl;
    }
  }

} // namespace MueLu

#endif //ifdef HAVE_MUELU_IFPACK2
#endif // MUELU_IFPACK2SMOOTHER_DEF_HPP
