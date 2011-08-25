#ifndef MUELU_SAPFACTORY_HPP
#define MUELU_SAPFACTORY_HPP

#include <iostream>

#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_CrsOperator.hpp>
#include <Xpetra_Vector.hpp>
#include <Xpetra_VectorFactory.hpp>

#ifdef HAVE_MUELU_EPETRA_AND_EPETRAEXT
#include "EpetraExt_MatrixMatrix.h"
#endif

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_PFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

/*!
  @class SaPFactory class.
  @brief Factory for building Smoothed Aggregation prolongators.

  Right now this factory assumes a 1D problem.  Aggregation is hard-coded to divide
  the # fine dofs by 3.
*/

template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps> //TODO: or BlockSparseOp ?
class SaPFactory : public PFactory {
#include "MueLu_UseShortNames.hpp"

  template<class AA, class BB, class CC, class DD, class EE>
  inline friend std::ostream& operator<<(std::ostream& os, SaPFactory<AA,BB,CC,DD, EE> &factory);

  private:
/*
     TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
     RCP<MueLu::AggregationFactory<LO,GO,NO,LMO> > AggFact_;
     CoalesceFact_
     TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
*/
     RCP<PFactory> initialPFact_;
     std::string diagonalView_;
     bool doQR_;
     Scalar dampingFactor_;
     bool useAFiltered_;
     bool reUseP_;
     bool reUsePtent_;

  public:
    //! @name Constructors/Destructors.
    //@{

    /*! @brief Default constructor.

    */
    SaPFactory() : diagonalView_("current"),
                   //AggFact_(Teuchos::null),
                   doQR_(false), dampingFactor_(4./3), useAFiltered_(false), reUseP_(false),
                   reUsePtent_(false)
                   //, PFactory::reUseGraph_(false), PFactory::reUseAggregates_(false)
    {
      initialPFact_ = rcp(new TentativePFactory()),
      PFactory::reUseGraph_=false;
      PFactory::reUseAggregates_=false;
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "SaPFactory: Instantiating a new factory" << std::endl;
    }

    /*! @brief Constructor.

        User can supply a factory for generating the tentative prolongator.
    */
    SaPFactory(RCP<PFactory> InitialPFact) : initialPFact_(InitialPFact), diagonalView_("current"),
                   //AggFact_(Teuchos::null),
                   doQR_(false), dampingFactor_(4./3), useAFiltered_(false), reUseP_(false),
                   reUsePtent_(false)
                   //, PFactory::reUseGraph_(false), PFactory::reUseAggregates_(false)
    {
      PFactory::reUseGraph_=false;
      PFactory::reUseAggregates_=false;
    }

    //! Destructor.
    virtual ~SaPFactory() {}
    //@}

    //! @name Build methods.
    //@{

    /*!
      @brief Build method.

      Builds smoothed aggregation prolongator and returns it in <tt>coarseLevel</tt>.
      //FIXME what does the return code mean (unclear in MueMat)?
      //FIXME how should nullspace be stored?
    */
  bool Build(Level& fineLevel, Level &coarseLevel) const {
    return BuildP(fineLevel,coarseLevel);
  }

    bool BuildP(Level &fineLevel, Level &coarseLevel) const {
      Teuchos::OSTab tab(this->out_);

      std::ostringstream buf; buf << coarseLevel.GetLevelID();
      RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("SaPFactory::BuildP_"+buf.str()));
      timer->start(true);

      RCP<Operator> finalP;

      if (reUseP_) {
        if (coarseLevel.Get< RCP<Operator> >("P") == Teuchos::null)
          throw(std::runtime_error("SaPFactory: you have asked to reuse P, but it doesn't exist"));
        if (coarseLevel.IsAvailable("Nullspace") == false)
          throw(std::runtime_error("SaPFactory: you have asked to reuse cnull, but it doesn't exist"));
        return true;
      }

      //TODO get or generate fine grid nullspace here
      RCP<MultiVector> fineNullspace;
      if (fineLevel.IsAvailable("Nullspace")) {
        fineLevel.Get("Nullspace",fineNullspace);
      } else {
        //TODO add this functionality
        //throw(Exceptions::NotImplemented("SaPFactory.Build():  nullspace generation not implemented yet"));
        std::cout << "nullspace generation not implemented yet" << std::endl;
      }

      coarseLevel.Request("Ptent");
      coarseLevel.Request("Nullspace");
      initialPFact_->BuildP(fineLevel,coarseLevel);
      RCP<Operator> Ptent;
      coarseLevel.Get("Ptent",Ptent);


      //TMP, to force desallocation of Ptent
      coarseLevel.Set("Ptent", Teuchos::null);

      RCP<MultiVector> coarseNullspace;
      if (reUsePtent_) {
        try {
          coarseLevel.Get("Ptent",Ptent); //FIXME throws an error, replace with recomputation
          coarseLevel.Get("Nullspace",coarseNullspace); //FIXME throws an error, replace with recomputation
        }
        catch(...) {
          throw(Exceptions::NotImplemented("SaPFactory.Build(): regeneration of Ptent/nullspace not implemented yet"));
        }
      }

      if (coarseLevel.IsRequested("Ptent"))
        coarseLevel.Set("Ptent",Ptent);
      

      //Build final prolongator

      //FIXME Xpetra::Operator should calculate/stash max eigenvalue
      //FIXME SC lambdaMax = Op->GetDinvALambda();

      if (dampingFactor_ != 0) {

        RCP<Teuchos::Time> sapTimer;
        //sapTimer = rcp(new Teuchos::Time("SaPFactory:I * Ptent"));
        //sapTimer->start(true);
        //Teuchos::ParameterList matrixList;
        //RCP<Operator> I = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map,CrsOperator>("Identity",fineLevel.Get< RCP<Operator> >("A")->getRowMap(),matrixList);
        //RCP<Operator> newPtent = Utils::TwoMatrixMultiply(I,false,Ptent,false);
        //Ptent = newPtent; //I tried a checkout of the original Ptent, and it seems to be gone now (which is good)
        //sapTimer->stop();
        //MemUtils::ReportTimeAndMemory(*sapTimer, *(Op->getRowMap()->getComm()));

        RCP< Operator > Op = fineLevel.Get< RCP<Operator> >("A");
        sapTimer = rcp(new Teuchos::Time("SaPFactory:APtent_"+buf.str()));
        sapTimer->start(true);

        //JJH -- If I switch doFillComplete to false, the resulting matrix seems weird when printed with describe.
        //JJH -- The final prolongator is wrong, to boot.  So right now, I fillComplete AP, but avoid fillComplete
        //JJH -- in the scaling.  Long story short, we're doing 2 fillCompletes, where ideally we'd do just one.
        bool doFillComplete=true;
        bool optimizeStorage=false;
        RCP<Operator> AP = Utils::TwoMatrixMultiply(Op,false,Ptent,false,doFillComplete,optimizeStorage);
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(Op->getRowMap()->getComm()));

        sapTimer = rcp(new Teuchos::Time("SaPFactory:Dinv_APtent_"+buf.str()));
        sapTimer->start(true);
        doFillComplete=false;
        optimizeStorage=false;
        Teuchos::ArrayRCP<SC> diag = Utils::GetMatrixDiagonal(Op);
        Utils::MyOldScaleMatrix(AP,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(Op->getRowMap()->getComm()));

        sapTimer = rcp(new Teuchos::Time("SaPFactory:eigen_estimate_"+buf.str()));
        sapTimer->start(true);
        Scalar lambdaMax = Utils::PowerMethod(*Op, true, (LO) 10,(Scalar)1e-4);
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(Op->getRowMap()->getComm()));
        RCP<const Teuchos::Comm<int> > comm = Op->getRowMap()->getComm();
        if (comm->getRank() == 0)
          std::cout << "damping factor = " << dampingFactor_/lambdaMax << " ("
                    << dampingFactor_ << " / " << lambdaMax << ")" << std::endl;

        sapTimer = rcp(new Teuchos::Time("SaPFactory:Pt_plus_DinvAPtent_"+buf.str()));
        sapTimer->start(true);

        bool doTranspose=false; 
        if (AP->isFillComplete())
          Utils::TwoMatrixAdd(Ptent,doTranspose,1.0,AP,doTranspose,-dampingFactor_/lambdaMax,finalP);
        else {
          Utils::TwoMatrixAdd(Ptent,doTranspose,1.0,AP,-dampingFactor_/lambdaMax);
          finalP = AP;
        }
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(Op->getRowMap()->getComm()));

        sapTimer = rcp(new Teuchos::Time("SaPFactory:finalP_fillComplete_"+buf.str()));
        sapTimer->start(true);
        finalP->fillComplete( Ptent->getDomainMap(), Ptent->getRangeMap() );
        sapTimer->stop();
        MemUtils::ReportTimeAndMemory(*sapTimer, *(Op->getRowMap()->getComm()));
      }
      else {
        finalP = Ptent;
      }

      coarseLevel.Set("P", finalP);
      //coarseLevel.Set("Nullspace",coarseNullspace);

      //Utils::MatrixPrint(finalP);

      timer->stop();
      MemUtils::ReportTimeAndMemory(*timer, *(finalP->getRowMap()->getComm()));

      return true;
    } //Build()
    //@}

    //! @name Set methods.
    //@{

    //! Set prolongator smoother damping factor.
    void SetDampingFactor(Scalar dampingFactor) {
      dampingFactor_ = dampingFactor;
    }

    void TentativeWithQR(bool value) {
      doQR_ = value;
    }

    //! Change view of diagonal.
    void SetDiagonalView(std::string const& diagView) {
      diagonalView_ = diagView;
    }

    void SetUseAFiltered(bool value) {
      throw(Exceptions::NotImplemented("SetUseAFiltered not fully implemented"));
      useAFiltered_ = value;
      //FIXME add to needs?
    }

    void ReUseP(bool value) {
      reUseP_ = value;
    }

    void ReUsePtent(bool value) {
      reUsePtent_ = value;
    }

    //@}

    //! @name Get methods.
    //@{

    //! Returns prolongator smoother damping factor.
    Scalar GetDampingFactor() {
      return dampingFactor_;
    }

    bool TentativeWithQR() {
      return doQR_;
    }

    bool ReUseP() {
      return reUseP_;
    }

    bool ReUsePtent() {
      return reUsePtent_;
    }

    //! Returns current view of diagonal.
    std::string GetDiagonalView() {
      return diagonalView_;
    }

    //@}
/*
//TODO
function [this] = SaPFactory(CoalesceFact,AggFact, diagonalView) //copy ctor
function SetDiagonalView(this, diagonalView)
*/

}; //class SaPFactory

//! Friend print function.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::ostream& operator<<(std::ostream& os, SaPFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> &factory) {
  os << "Printing an SaPFactory object" << std::endl;
  return os;
}

} //namespace MueLu

#define MUELU_SAPFACTORY_SHORT

#endif //ifndef MUELU_SAPFACTORY_HPP
