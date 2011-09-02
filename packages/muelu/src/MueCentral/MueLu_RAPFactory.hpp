#ifndef MUELU_RAPFACTORY_HPP
#define MUELU_RAPFACTORY_HPP

#include <iostream>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_TwoLevelFactoryBase.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {
/*!
  @class RAPFactory class.
  @brief Factory for building coarse matrices.
*/
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class RAPFactory : public TwoLevelFactoryBase {

#include "MueLu_UseShortNames.hpp"


  //JG to JJH: use Teuchos::Describable instead ?
  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, RAPFactory<AA,BB,CC,DD,EE> &factory);

  public:
    //@{ Constructors/Destructors.
  RAPFactory(RCP<FactoryBase> PFact = Teuchos::null, RCP<FactoryBase> RFact = Teuchos::null) 
    : PFact_(PFact), RFact_(RFact),
      implicitTranspose_(false) {}

    virtual ~RAPFactory() {}
    //@}

    //! Input
    //@{

    void DeclareInput(Level &fineLevel, Level &coarseLevel) const { }

    //@}

    //@{ Build methods.
    bool Build(Level &fineLevel, Level &coarseLevel) const {  //FIXME make fineLevel const!!

      std::ostringstream buf; buf << coarseLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("RAP::Build_"+buf.str()));
      timer->start(true);

      Teuchos::OSTab tab(this->getOStream());
      //MueLu_cout(Teuchos::VERB_LOW) << "call the Epetra matrix-matrix multiply here" << std::endl;
      RCP<Operator> P = coarseLevel.Get< RCP<Operator> >("P", PFact_);
      RCP<Operator> A = fineLevel.Get< RCP<Operator> >("A");
RCP<Teuchos::Time> apTimer = rcp(new Teuchos::Time("RAP::A_times_P_"+buf.str()));
apTimer->start(true);
      RCP<Operator> AP = Utils::TwoMatrixMultiply(A,false,P,false);
apTimer->stop();
MemUtils::ReportTimeAndMemory(*apTimer, *(P->getRowMap()->getComm()));
      //std::string filename="AP.dat";
      //Utils::Write(filename,AP);

      if (implicitTranspose_) {
        //RCP<Operator> RA = Utils::TwoMatrixMultiply(P,true,A,false);
        //filename = "PtA.dat";
        //Utils::Write(filename,AP);
        RCP<Operator> RAP = Utils::TwoMatrixMultiply(P,true,AP,false);
        coarseLevel.Set("A",RAP);
      } else {
        RCP<Operator> R = coarseLevel.Get< RCP<Operator> >("R", RFact_);
RCP<Teuchos::Time> rapTimer = rcp(new Teuchos::Time("RAP::R_times_AP_"+buf.str()));
rapTimer->start(true);
        RCP<Operator> RAP = Utils::TwoMatrixMultiply(R,false,AP,false);
rapTimer->stop();
MemUtils::ReportTimeAndMemory(*rapTimer, *(P->getRowMap()->getComm()));

        coarseLevel.Set("A", RAP, this);
      }

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(P->getRowMap()->getComm()));

      return true;
    }
    //@}

    void SetImplicitTranspose(bool const &implicit) {
      implicitTranspose_ = implicit;
    }


private:
  //! P Factory
  RCP<FactoryBase> PFact_;

  //! R Factory
  RCP<FactoryBase> RFact_;
  
  bool implicitTranspose_;

}; //class RAPFactory

//! Friend print method.
  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, RAPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> &factory) {
  os << "Printing RAPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_RAPFACTORY_SHORT

#endif //ifndef MUELU_RAPFACTORY_HPP
