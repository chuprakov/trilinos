#ifndef MUELU_TRANSPFACTORY_HPP
#define MUELU_TRANSPFACTORY_HPP

#include <iostream>
#include "Cthulhu_CrsOperator.hpp"
#include "MueLu_RFactory.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

  /*!
    @class TransPFactory class.
    @brief Factory for building restriction operators.

    This factory currently depends on an underlying matrix-matrix multiply with the identity
    matrix to do the transpose.  This should probably be fixed at some point.
  */

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  class TransPFactory : public RFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {

#include "MueLu_UseShortNames.hpp"

    template<class AA, class BB, class CC, class DD, class EE>
    inline friend std::ostream& operator<<(std::ostream& os, TransPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    TransPFactory() {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "TransPFactory: Instantiating a new factory" << std::endl;
    }

    //! Destructor.
    virtual ~TransPFactory() {}
    //@}

    //! @name Build methods.
    //@{
    bool BuildR(Level & fineLevel, Level & coarseLevel) {
      Teuchos::OSTab tab(this->out_);
      Teuchos::ParameterList matrixList;
      RCP<Operator> P = coarseLevel.GetP();
      //doesn't work -- bug in EpetraExt?
      //RCP<CrsOperator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getRangeMap(),matrixList);

      RCP<CrsOperator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>("Identity",P->getDomainMap(),matrixList);

      //RCP<Operator> R = Utils::TwoMatrixMultiply(P,I,true); //doesn't work -- bug in EpetraExt?
      RCP<Operator> R = Utils::TwoMatrixMultiply(I,P,false,true);
      coarseLevel.SetR(R);
      return true;
    }
    //@}

    //! @name Set methods.
    //@{
    void UsePtent(bool ToF) {
      throw(Exceptions::NotImplemented("TransPFactory.UsePtent()")); //TODO
    }
    //@}



  }; //class TransPFactory

  //! Friend print function.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::ostream& operator<<(std::ostream& os, TransPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
    os << "Printing a TransPFactory object" << std::endl;
    return os;
  }

} //namespace MueLu

#define MUELU_TRANSPFACTORY_SHORT

#endif //ifndef MUELU_TRANSPFACTORY_HPP
