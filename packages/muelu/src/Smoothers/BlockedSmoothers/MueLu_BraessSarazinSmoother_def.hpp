/*
 * MueLu_BraessSarazinSmoother_def.hpp
 *
 *  Created on: Apr 16, 2012
 *      Author: wiesner
 */

#ifndef MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_
#define MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_BlockedCrsOperator.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include "MueLu_BraessSarazinSmoother_decl.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_SmootherBase.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BraessSarazinSmoother(LocalOrdinal sweeps, Scalar omega, RCP<FactoryBase> AFact)
    : type_("Braess Sarazin"), nSweeps_(sweeps), omega_(omega), AFact_(AFact), A_(Teuchos::null)
  {
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~BraessSarazinSmoother() {}

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::AddFactoryManager(RCP<const FactoryManagerBase> FactManager) {
    FactManager_ = FactManager;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &currentLevel) const {
    //Global A operator for getting the blocks.
    currentLevel.DeclareInput("A", AFact_.get());
    //Smoother for the Schur complement.
    currentLevel.DeclareInput("PreSmoother",FactManager_->GetFactory("PreSmoother").get());
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Setup(Level &currentLevel) {
    //*********************************************
    // Setup routine can be summarized in 4 steps:
    // - Set the map extractors
    // - Set the blocks
    // - Create and set the inverse of the diagonal of F (note that this means if a different approximation of F^-1 is desired, it must be coded again)
    // - Set the smoother for the Schur Complement

//    typedef Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsMatrixClass; //TODO
//    typedef Xpetra::CrsOperator<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsOperatorClass; //TODO

    FactoryMonitor m(*this, "Setup blocked Braess-Sarazin Smoother", currentLevel);

    if (SmootherPrototype::IsSetup() == true)
            this->GetOStream(Warnings0, 0) << "Warning: MueLu::Amesos2Smoother::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = currentLevel.Get<RCP<Operator> > ("A", AFact_.get());

    RCP<BlockedCrsOperator> bA = Teuchos::rcp_dynamic_cast<BlockedCrsOperator>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA == Teuchos::null, Exceptions::BadCast, "MueLu::BlockedPFactory::Build: input matrix A is not of type BlockedCrsOperator! error.");

    // store map extractors
    rangeMapExtractor_ = bA->getRangeMapExtractor();
    domainMapExtractor_ = bA->getDomainMapExtractor();

    //***************************************************************************
    // Store the blocks in local member variables
    //
    Teuchos::RCP<CrsMatrix> A00 = bA->getMatrix(0, 0);
    Teuchos::RCP<CrsMatrix> A01 = bA->getMatrix(0, 1);
    Teuchos::RCP<CrsMatrix> A10 = bA->getMatrix(1, 0);
    Teuchos::RCP<CrsMatrix> A11 = bA->getMatrix(1, 1);

    Teuchos::RCP<CrsOperator> Op00 = Teuchos::rcp(new CrsOperator(A00));
    Teuchos::RCP<CrsOperator> Op01 = Teuchos::rcp(new CrsOperator(A01));
    Teuchos::RCP<CrsOperator> Op10 = Teuchos::rcp(new CrsOperator(A10));
    Teuchos::RCP<CrsOperator> Op11 = Teuchos::rcp(new CrsOperator(A11));

    F_ = Teuchos::rcp_dynamic_cast<Operator>(Op00);
    G_ = Teuchos::rcp_dynamic_cast<Operator>(Op01);
    D_ = Teuchos::rcp_dynamic_cast<Operator>(Op10);
    Z_ = Teuchos::rcp_dynamic_cast<Operator>(Op11);

    //******************************************************************
    // Create the inverse of the diagonal of F

    // DA: it is better if diagFinv_ is a vector, so as to be used in apply, eq 8.13
    RCP<Vector> diagFVector = VectorFactory::Build(F_->getRowMap());
    F_->getLocalDiagCopy(*diagFVector);       // extract diagonal of F
    diagFVector->reciprocal(*diagFVector);    // build reciprocal
    diagFinv_ = diagFVector;

    //********************************************************************
    // Set the Smoother
    smoo_ = currentLevel.Get<RCP<SmootherBase> > ("PreSmoother", FactManager_->GetFactory("PreSmoother").get());

    SmootherPrototype::IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Apply(MultiVector &X, MultiVector const &B, bool const &InitialGuessIsZero) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(SmootherPrototype::IsSetup() == false, Exceptions::RuntimeError, "MueLu::BraessSarazinSmoother::Apply(): Setup() has not been called");

    // TODO DO NOT use v and p, because it will be used in generic blocked matrices!
    //Range map is ok because in vtemp the result from 1/omega * Fhatinv * rvel is stored.
    RCP<MultiVector> vtemp = MultiVectorFactory::Build(F_->getRowMap(),1);
    RCP<MultiVector> ptemp = MultiVectorFactory::Build(Z_->getRowMap(),1);

    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > residual = MultiVectorFactory::Build(B.getMap(), B.getNumVectors());
    RCP<MultiVector> rcpX = Teuchos::rcpFromRef(X);

#define PRINT_RESIDUAL 0
#if PRINT_RESIDUAL
     //This is necessary to print out the residual
     Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> normsvel(1);
     Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> normspre(1);


#endif

    //run for each sweep
    for (LocalOrdinal run = 0; run < nSweeps_; ++run) {
      //************************
      // Calculate residual
      // note: A_ is the full blocked operator
      // r = B - A*X
      residual->update(1.0,B,0.0); // residual = B
      A_->apply(X, *residual, Teuchos::NO_TRANS, -1.0, 1.0); // residual = B-A*X

      //************************
      // split the data
      // extract corresponding subvectors from X and residual

      Teuchos::RCP<MultiVector> rvel = rangeMapExtractor_->ExtractVector(residual, 0);
      Teuchos::RCP<MultiVector> rpre = rangeMapExtractor_->ExtractVector(residual, 1);


#if PRINT_RESIDUAL
      //print the norm
      rvel->norm2(normsvel);
      rpre->norm2(normspre);
      //This format is in such a way it can be used for plotting.
      this->GetOStream(Runtime0, 0) << run << "\t" << normsvel[0] << "\t" << normspre[0] << std::endl;
#endif

      //DA: What does this do?
      Teuchos::RCP<MultiVector> tXvel = domainMapExtractor_->ExtractVector(rcpX, 0);
      Teuchos::RCP<MultiVector> tXpre = domainMapExtractor_->ExtractVector(rcpX, 1);

      //**********************
      // Calculate qrhs = rpre - D * vtemp (equation 8.14)
      RCP<MultiVector> qrhs = MultiVectorFactory::Build(Z_->getRowMap(),1);
      vtemp->putScalar(0.0);
      vtemp->elementWiseMultiply(1.0/omega_,*diagFinv_,*rvel,0.0); // vtemp = 1/omega*Fhatinv * rvel (equation 8.13)

      //This is used for printing intermediate results (uncomment to use)
      //Teuchos::Array<typename Teuchos::ScalarTraits<Scalar>::magnitudeType> normsall(1);
      //vtemp->norm2(normsall);
      //this->GetOStream(Runtime0, 0) << "BSS::Setup, ||vtemp|| = " << normsall[0];

      D_->apply(*vtemp,*qrhs); // qrhs = D*vtemp (intermediate step)
      qrhs->update(1.0,*rpre,-1.0); // qrhs = rpre - D*vtemp

      //***********************
      //Pressure correction, using the preconditioner
      //RCP<MultiVector> q = MultiVectorFactory::Build(Z_->getDomainMap(),1); // TODO think about this.
      RCP<MultiVector> q = MultiVectorFactory::Build(Z_->getRowMap(),1); // TODO think about this.
      smoo_->Apply(*q,*qrhs);


      //*****************************
      //Update

      vtemp->putScalar(0.0);
      G_->apply(*q,*vtemp);
      vtemp->update(1.0,*rvel,-1.0); //velres - G*q
      RCP<MultiVector> vx = MultiVectorFactory::Build(F_->getRowMap(),1);
      vx->elementWiseMultiply(1.0/omega_,*diagFinv_,*vtemp,1.0); //vx = 1/(omega)*DFhatinv*(velres-G*q)


      //**********************
      // update with values

      tXvel->update(1.0,*vx,1.0);
      tXpre->update(1.0,*q,1.0);

      domainMapExtractor_->InsertVector(tXvel, 0, rcpX);
      domainMapExtractor_->InsertVector(tXpre, 1, rcpX);
    }

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> > BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Copy() const {
    return rcp( new BraessSarazinSmoother(*this) );
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void BraessSarazinSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0) {
      out0 << "Prec. type: " << type_ << /*" Sweeps: " << nSweeps_ << " damping: " << omega_ <<*/ std::endl;
    }

    if (verbLevel & Debug) {
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
    }
  }

} // namespace MueLu

#endif /* MUELU_BRAESSSARAZINSMOOTHER_DEF_HPP_ */
