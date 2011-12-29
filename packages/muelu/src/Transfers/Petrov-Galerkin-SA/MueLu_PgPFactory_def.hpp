#ifndef MUELU_PGPFACTORY_DEF_HPP_
#define MUELU_PGPFACTORY_DEF_HPP_

#include "MueLu_PgPFactory_decl.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SingleLevelFactoryBase.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"

#include "Xpetra_VectorFactory.hpp"
#include "Xpetra_ImportFactory.hpp"
#include "Xpetra_ExportFactory.hpp"

namespace MueLu {

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::PgPFactory(RCP<PFactory> InitialPFact, RCP<SingleLevelFactoryBase> AFact)
: initialPFact_(InitialPFact), AFact_(AFact),
  diagonalView_("current") {
  if(initialPFact_ == Teuchos::null)
  {
    // use tentative P factory as default
    initialPFact_ = rcp(new TentativePFactory());
  }

  min_norm_ = DINVANORM;

  bReUseRowBasedOmegas_ = false;
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::~PgPFactory() {}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetDiagonalView(std::string const& diagView) {
  diagonalView_ = diagView;
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetMinimizationMode(MinimizationNorm minnorm) { min_norm_ = minnorm; }

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
MueLu::MinimizationNorm PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMinimizationMode() { return min_norm_; }

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
std::string PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetDiagonalView() {
  return diagonalView_;
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::DeclareInput(Level &fineLevel, Level &coarseLevel) const {
  fineLevel.DeclareInput("A",AFact_.get(),this);
  coarseLevel.DeclareInput("P",initialPFact_.get(),this);

  /* If PgPFactory is reusing the row based damping parameters omega for
   * restriction, it has to request the data here.
   * we have the following scenarios:
   * 1) Reuse omegas:
   * PgPFactory.DeclareInput for prolongation mode requests A and P0
   * PgPFactory.DeclareInput for restriction mode requests A, P0 and RowBasedOmega (call triggered by GenericRFactory)
   * PgPFactory.Build for prolongation mode calculates RowBasedOmega and stores it (as requested)
   * PgPFactory.Build for restriction mode reuses RowBasedOmega (and Releases the data with the Get call)
   * 2) do not reuse omegas
   * PgPFactory.DeclareInput for prolongation mode requests A and P0
   * PgPFactory.DeclareInput for restriction mode requests A and P0
   * PgPFactory.Build for prolongation mode calculates RowBasedOmega for prolongation operator
   * PgPFactory.Build for restriction mode calculates RowBasedOmega for restriction operator
   */
  if( bReUseRowBasedOmegas_ == true && restrictionMode_ == true ) {
    coarseLevel.DeclareInput("RowBasedOmega", this, this); // RowBasedOmega is calculated by this PgPFactory and requested by this PgPFactory
  }
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::Build(Level& fineLevel, Level &coarseLevel) const {
  Teuchos::RCP<Teuchos::FancyOStream> fos = Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout));
  std::ostringstream buf; buf << coarseLevel.GetLevelID();
  RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("PgPFactory::BuildP_"+buf.str()));
  timer->start(true);

  // Level Get
  RCP<Operator> Ptent = coarseLevel.Get< RCP<Operator> >("P", initialPFact_.get());
  RCP<Operator> A     = fineLevel.  Get< RCP<Operator> >("A", AFact_.get());

  /////////////////// switch from A to A^T in restriction mode (necessary as long as implicit transpose not working for Epetra)
  if(restrictionMode_)
    A = Utils2::Transpose(A,true); // build transpose of A explicitely

  Monitor m(*this, "prolongator smoothing (PG-AMG)");

  /////////////////// calculate D^{-1} A Ptent (needed for smoothing)
  bool doFillComplete=true;
  bool optimizeStorage=false;
  RCP<Operator> DinvAP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

  doFillComplete=true;
  optimizeStorage=false;
  Teuchos::ArrayRCP<Scalar> diag = Utils::GetMatrixDiagonal(A);
  Utils::MyOldScaleMatrix(DinvAP0,diag,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

  /////////////////// calculate local damping factors omega

  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > RowBasedOmega = Teuchos::null;

  if(restrictionMode_ == false || bReUseRowBasedOmegas_ == false) {
    // if in prolongation mode: calculate row based omegas
    // if in restriction mode: calculate omegas only if row based omegas are not used from prolongation mode
    ComputeRowBasedOmega(fineLevel, coarseLevel, A, Ptent, DinvAP0, RowBasedOmega);
  } // if(bReUseRowBasedOmegas == false)
  else  {
    // reuse row based omegas, calculated by this factory in the run before (with restrictionMode_ == false)
    RowBasedOmega = coarseLevel.Get<Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > >("RowBasedOmega", this /*MueLu::NoFactory::get()*/);

    // RowBasedOmega is now based on row map of A (not transposed)
    // for restriction we use A^T instead of A
    // -> recommunicate row based omega

    // exporter: overlapping row map to nonoverlapping domain map (target map is unique)
    // since A is already transposed we use the RangeMap of A
    Teuchos::RCP<const Export> exporter =
        Xpetra::ExportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(RowBasedOmega->getMap(), A->getRangeMap());

    Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > noRowBasedOmega =
        Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(A->getRangeMap());

    noRowBasedOmega->doExport(*RowBasedOmega,*exporter,Xpetra::INSERT);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer =
        Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(A->getRangeMap(),A->getRowMap());

    // doImport target->doImport(*source, importer, action)
    RowBasedOmega->doImport(*noRowBasedOmega,*importer,Xpetra::INSERT);
  }

  Teuchos::ArrayRCP< Scalar > RowBasedOmega_local = RowBasedOmega->getDataNonConst(0);

  /////////////////// prolongator smoothing using local damping parameters omega
  RCP<Operator> P_smoothed = Teuchos::null;
  Utils::MyOldScaleMatrix(DinvAP0,RowBasedOmega_local,false,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

  Utils::TwoMatrixAdd(Ptent, false, Teuchos::ScalarTraits<Scalar>::one(),
      DinvAP0, false, -Teuchos::ScalarTraits<Scalar>::one(),
      P_smoothed);
  P_smoothed->fillComplete(Ptent->getDomainMap(), Ptent->getRangeMap());

  //////////////////// store results in Level

  // Level Set
  if(!restrictionMode_)
  {
    // prolongation factory is in prolongation mode
    coarseLevel.Set("P", P_smoothed, this);
  }
  else
  {
    // prolongation factory is in restriction mode
    RCP<Operator> R = Utils2::Transpose(P_smoothed,true); // use Utils2 -> specialization for double
    coarseLevel.Set("R", R, this);
  }

  timer->stop();
  MemUtils::ReportTimeAndMemory(*timer, *(P_smoothed->getRowMap()->getComm()));

}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ComputeRowBasedOmega(Level& fineLevel, Level &coarseLevel, const RCP<Operator>& A, const RCP<Operator>& P0, const RCP<Operator>& DinvAP0, RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > & RowBasedOmega) const {
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Numerator = Teuchos::null;
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > Denominator = Teuchos::null;

  switch (min_norm_)
  {
  case ANORM: {
    // MUEMAT mode (=paper)
    // Minimize with respect to the (A)' A norm.
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' (A' A) D^{-1} A P0)
    //   omega =   ------------------------------------------
    //             diag( P0' A' D^{-1}' ( A'  A) D^{-1} A P0)
    //
    // expensive, since we have to recalculate AP0 due to the lack of an explicit scaling routine for DinvAP0

    // calculate A * P0
    bool doFillComplete=true;
    bool optimizeStorage=false;
    RCP<Operator> AP0 = Utils::TwoMatrixMultiply(A,false,P0,false,doFillComplete,optimizeStorage);

    // compute A * D^{-1} * A * P0
    RCP<Operator> ADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAP0,false,doFillComplete,optimizeStorage);

    Numerator =   Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(ADinvAP0->getColMap(),true);
    Denominator = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(ADinvAP0->getColMap(),true);
    MultiplyAll(AP0, ADinvAP0, Numerator);
    MultiplySelfAll(ADinvAP0, Denominator);
  }
  break;
  case L2NORM: {

    // ML mode 1 (cheapest)
    // Minimize with respect to L2 norm
    //                  diag( P0' D^{-1} A P0)
    //   omega =   -----------------------------
    //             diag( P0' A' D^{-1}' D^{-1} A P0)
    //
    Numerator =   Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(DinvAP0->getColMap(),true);
    Denominator = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(DinvAP0->getColMap(),true);
    MultiplyAll(P0, DinvAP0, Numerator);
    MultiplySelfAll(DinvAP0, Denominator);
  }
  break;
  case DINVANORM: {
    // ML mode 2
    // Minimize with respect to the (D^{-1} A)' D^{-1} A norm.
    // Need to be smart here to avoid the construction of A' A
    //
    //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //   omega =   --------------------------------------------------------
    //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
    //

    // compute D^{-1} * A * D^{-1} * A * P0
    bool doFillComplete=true;
    bool optimizeStorage=false;
    Teuchos::ArrayRCP<Scalar> diagA = Utils::GetMatrixDiagonal(A);
    RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAP0,false,doFillComplete,optimizeStorage);
    Utils::MyOldScaleMatrix(DinvADinvAP0,diagA,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

    Numerator =   Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(DinvADinvAP0->getColMap(),true);
    Denominator = Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(DinvADinvAP0->getColMap(),true);

    MultiplyAll(DinvAP0, DinvADinvAP0, Numerator);
    MultiplySelfAll(DinvADinvAP0, Denominator);
  }
  break;
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::Build: minimization mode not supported. error");
  }


  //////////// build Column based omegas /////////////
  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > ColBasedOmega =
      Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Numerator->getMap()/*DinvAP0->getColMap()*/,true);

  ColBasedOmega->putScalar(-666*Teuchos::ScalarTraits<Scalar>::one());

  Teuchos::ArrayRCP< const Scalar > Numerator_local = Numerator->getData(0);
  Teuchos::ArrayRCP< const Scalar > Denominator_local = Denominator->getData(0);
  Teuchos::ArrayRCP< Scalar >       ColBasedOmega_local = ColBasedOmega->getDataNonConst(0);
  LocalOrdinal zero_local = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar min_local = Teuchos::ScalarTraits<Scalar>::one() * 1000000;
  Scalar max_local = Teuchos::ScalarTraits<Scalar>::zero();
  for(LocalOrdinal i = 0; i < Teuchos::as<LocalOrdinal>(Numerator->getLocalLength()); i++) {
    ColBasedOmega_local[i] = Numerator_local[i] / Denominator_local[i];
    if(ColBasedOmega_local[i] < Teuchos::ScalarTraits<Scalar>::zero()) { // negative omegas are not valid. set them to zero
      ColBasedOmega_local[i] = Teuchos::ScalarTraits<Scalar>::zero();
      zero_local++; // count zero omegas
    }
    if(ColBasedOmega_local[i] < min_local) { min_local = ColBasedOmega_local[i]; }
    if(ColBasedOmega_local[i] > max_local) { max_local = ColBasedOmega_local[i]; }
  }

  { // be verbose
    GlobalOrdinal zero_all;
    Scalar min_all;
    Scalar max_all;
    sumAll(A->getRowMap()->getComm(),zero_local,zero_all);
    minAll(A->getRowMap()->getComm(),min_local, min_all);
    maxAll(A->getRowMap()->getComm(),max_local, max_all);

    GetOStream(MueLu::Statistics1,0) << "PgPFactory: smoothed aggregation (scheme: ";
    switch (min_norm_)
    {
    case ANORM:     { GetOStream(MueLu::Statistics1,0) << "Anorm)"     << std::endl;   }   break;
    case L2NORM:    { GetOStream(MueLu::Statistics1,0) << "L2norm)"    << std::endl;   }   break;
    case DINVANORM: { GetOStream(MueLu::Statistics1,0) << "DinvAnorm)" << std::endl;   }    break;
    default:          GetOStream(MueLu::Statistics1,0) << "unknown)" << std::endl;
    }
    GetOStream(MueLu::Statistics1,0) << "Damping parameter: min = " << min_all << ", max = " << max_all << ", (" << zero_all << " zeros out of " << ColBasedOmega->getGlobalLength() << " column-based omegas)" << std::endl;
  }

  if(coarseLevel.IsRequested("ColBasedOmega", this)) {
    coarseLevel.Set("ColBasedOmega", ColBasedOmega, this);
  }

  //////////// build Row based omegas /////////////
  // transform column based omegas to row based omegas
  RowBasedOmega =
      Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(DinvAP0->getRowMap(),true);

  RowBasedOmega->putScalar(-666*Teuchos::ScalarTraits<Scalar>::one());

  bool bAtLeastOneDefined = false;
  Teuchos::ArrayRCP< Scalar > RowBasedOmega_local = RowBasedOmega->getDataNonConst(0);
  for(LocalOrdinal row = 0; row<Teuchos::as<LocalOrdinal>(A->getNodeNumRows()); row++) {
    Teuchos::ArrayView<const LocalOrdinal> lindices;
    Teuchos::ArrayView<const Scalar> lvals;
    DinvAP0->getLocalRowView(row, lindices, lvals);
    bAtLeastOneDefined = false;
    for(size_t j=0; j<Teuchos::as<size_t>(lindices.size()); j++) {
      Scalar omega = ColBasedOmega_local[lindices[j]];
      if (omega != -666*Teuchos::ScalarTraits<Scalar>::one()) {
        bAtLeastOneDefined = true;
        if(RowBasedOmega_local[row] == -666*Teuchos::ScalarTraits<Scalar>::one())    RowBasedOmega_local[row] = omega;
        else if(omega < RowBasedOmega_local[row]) RowBasedOmega_local[row] = omega;
      }
    }
    if(bAtLeastOneDefined == true) {
      if(RowBasedOmega_local[row] < Teuchos::ScalarTraits<Scalar>::zero()) RowBasedOmega_local[row] = Teuchos::ScalarTraits<Scalar>::zero();
    }
  }

  if(coarseLevel.IsRequested("RowBasedOmega", this)) {
    coarseLevel.Set("RowBasedOmega", RowBasedOmega, this);
  }
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MultiplySelfAll(const RCP<Operator>& Op, Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& InnerProdVec) const {

  // note: InnerProdVec is based on column map of Op
  TEUCHOS_TEST_FOR_EXCEPTION(!InnerProdVec->getMap()->isSameAs(*Op->getColMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplySelfAll: map of InnerProdVec must be same as column map of operator. error");

  Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);

  for(size_t n=0; n<Op->getNodeNumRows(); n++) {
    Teuchos::ArrayView<const LocalOrdinal> lindices;
    Teuchos::ArrayView<const Scalar> lvals;
    Op->getLocalRowView(n, lindices, lvals);

    for(size_t i=0; i<Teuchos::as<size_t>(lindices.size()); i++) {
      InnerProd_local[lindices[i]] += lvals[i]*lvals[i];
    }
  }

  // exporter: overlapping map to nonoverlapping map (target map is unique)
  Teuchos::RCP<const Export> exporter =
      Xpetra::ExportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Op->getColMap(), Op->getDomainMap());

  Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > nonoverlap =
      Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(Op->getDomainMap());

  nonoverlap->doExport(*InnerProdVec,*exporter,Xpetra::ADD);

  // importer: nonoverlapping map to overlapping map

  // importer: source -> target maps
  Teuchos::RCP<const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer =
      Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(Op->getDomainMap(),Op->getColMap());

  // doImport target->doImport(*source, importer, action)
  InnerProdVec->doImport(*nonoverlap,*importer,Xpetra::INSERT);

}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::MultiplyAll(const RCP<Operator>& left, const RCP<Operator>& right, Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& InnerProdVec) const {

  TEUCHOS_TEST_FOR_EXCEPTION(!left->getDomainMap()->isSameAs(*right->getDomainMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: domain maps of left and right do not match. Error.");
  TEUCHOS_TEST_FOR_EXCEPTION(!left->getRowMap()->isSameAs(*right->getRowMap()), Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: row maps of left and right do not match. Error.");

  if(InnerProdVec->getMap()->isSameAs(*left->getColMap())) {

    Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);

    for(size_t n=0; n<left->getNodeNumRows(); n++)
    {
      Teuchos::ArrayView<const LocalOrdinal> lindices_left;
      Teuchos::ArrayView<const Scalar> lvals_left;
      left->getLocalRowView(n, lindices_left, lvals_left);

      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;
      right->getLocalRowView(n, lindices_right, lvals_right);

      for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
      {
        for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
        {
          GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
          GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
          if(left_gid == right_gid)
          {
            InnerProd_local[lindices_left[i]] += lvals_left[i]*lvals_right[j];
            break; // skip remaining gids of right operator
          }
        }
      }
    }

    // exporter: overlapping map to nonoverlapping map (target map is unique)
    Teuchos::RCP<const Export> exporter =
        Xpetra::ExportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(left->getColMap(), left->getDomainMap()); // TODO: change left to right?

    Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > nonoverlap =
        Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(left->getDomainMap()); // TODO: change left to right?

    nonoverlap->doExport(*InnerProdVec,*exporter,Xpetra::ADD);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer =
        Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(left->getDomainMap(),left->getColMap()); // TODO: change left to right?

    // doImport target->doImport(*source, importer, action)
    InnerProdVec->doImport(*nonoverlap,*importer,Xpetra::INSERT);
  }
  else if(InnerProdVec->getMap()->isSameAs(*right->getColMap())) {
    Teuchos::ArrayRCP< Scalar > InnerProd_local = InnerProdVec->getDataNonConst(0);

    for(size_t n=0; n<left->getNodeNumRows(); n++)
    {
      Teuchos::ArrayView<const LocalOrdinal> lindices_left;
      Teuchos::ArrayView<const Scalar> lvals_left;
      left->getLocalRowView(n, lindices_left, lvals_left);

      Teuchos::ArrayView<const LocalOrdinal> lindices_right;
      Teuchos::ArrayView<const Scalar> lvals_right;
      right->getLocalRowView(n, lindices_right, lvals_right);

      for(size_t i=0; i<Teuchos::as<size_t>(lindices_left.size()); i++)
      {
        for(size_t j=0; j<Teuchos::as<size_t>(lindices_right.size()); j++)
        {
          GlobalOrdinal left_gid = left->getColMap()->getGlobalElement(lindices_left[i]);
          GlobalOrdinal right_gid= right->getColMap()->getGlobalElement(lindices_right[j]);
          if(left_gid == right_gid)
          {
            InnerProd_local[lindices_right[j]] += lvals_left[i]*lvals_right[j];
            break; // skip remaining gids of right operator
          }
        }
      }
    }

    // exporter: overlapping map to nonoverlapping map (target map is unique)
    Teuchos::RCP<const Export> exporter =
        Xpetra::ExportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(right->getColMap(), right->getDomainMap()); // TODO: change left to right?

    Teuchos::RCP<Xpetra::Vector<Scalar,LocalOrdinal,GlobalOrdinal,Node> > nonoverlap =
        Xpetra::VectorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(right->getDomainMap()); // TODO: change left to right?

    nonoverlap->doExport(*InnerProdVec,*exporter,Xpetra::ADD);

    // importer: nonoverlapping map to overlapping map

    // importer: source -> target maps
    Teuchos::RCP<const Xpetra::Import<LocalOrdinal,GlobalOrdinal,Node> > importer =
        Xpetra::ImportFactory<LocalOrdinal,GlobalOrdinal,Node>::Build(right->getDomainMap(),right->getColMap()); // TODO: change left to right?

    // doImport target->doImport(*source, importer, action)
    InnerProdVec->doImport(*nonoverlap,*importer,Xpetra::INSERT);
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::PgPFactory::MultiplyAll: map of InnerProdVec must be same as column map of left or right operator? Error.");
  }
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::BuildP(Level &fineLevel, Level &coarseLevel) const {
  std::cout << "TODO: remove me" << std::endl;
}

template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
void PgPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::ReUseDampingParameters(bool bReuse) {
  bReUseRowBasedOmegas_ = bReuse;
}

#if 0
case ANORM: {
  // MUEMAT mode (=paper)
  // Minimize with respect to the (A)' A norm.
  // Need to be smart here to avoid the construction of A' A
  //
  //                   diag( P0' (A' A) D^{-1} A P0)
  //   omega =   ------------------------------------------
  //             diag( P0' A' D^{-1}' ( A'  A) D^{-1} A P0)
  //
  // expensive, since we have to recalculate AP0 due to the lack of an explicit scaling routine for DinvAP0

  // calculate A * Ptent
  bool doFillComplete=true;
  bool optimizeStorage=false;
  RCP<Operator> AP0 = Utils::TwoMatrixMultiply(A,false,Ptent,false,doFillComplete,optimizeStorage);

  // compute A * D^{-1} * A * P0
  RCP<Operator> ADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);

  Numerator = MultiplyAll(AP0, ADinvAP0, GID2localgid);
  Denominator = MultiplySelfAll(ADinvAP0, GID2localgid);
}
break;
case L2NORM: {
  // ML mode 1 (cheapest)
  // Minimize with respect to L2 norm
  //                  diag( P0' D^{-1} A P0)
  //   omega =   -----------------------------
  //             diag( P0' A' D^{-1}' D^{-1} A P0)
  //
  Numerator = MultiplyAll(Ptent, DinvAPtent, GID2localgid);
  Denominator = MultiplySelfAll(DinvAPtent, GID2localgid);
}
break;
case DINVANORM: {
  // ML mode 2
  // Minimize with respect to the (D^{-1} A)' D^{-1} A norm.
  // Need to be smart here to avoid the construction of A' A
  //
  //                   diag( P0' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
  //   omega =   --------------------------------------------------------
  //             diag( P0' A' D^{-1}' ( A' D^{-1}' D^{-1} A) D^{-1} A P0)
  //

  // compute D^{-1} * A * D^{-1} * A * P0
  bool doFillComplete=true;
  bool optimizeStorage=false;
  RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);
  Utils::MyOldScaleMatrix(DinvADinvAP0,diagA,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

  Numerator = MultiplyAll(DinvAPtent, DinvADinvAP0, GID2localgid);
  Denominator = MultiplySelfAll(DinvADinvAP0, GID2localgid);
}
break;
case ATDINVTPLUSDINVANORM: {
  // ML mode 3 (most expensive)
  //             diag( P0' ( A'D' + DA) D A P0)
  //   omega =   -----------------------------
  //             diag( P0'A'D' ( A'D' + DA) D A P0)
  //
  //             diag( DinvAP0'DinvAP0 + P0'DinvADinvAP0)
  //         =   -----------------------------
  //                2*diag( DinvADinvAP0'DinvAP0)
  //
  //

  // compute D^{-1} * A * D^{-1} * A * P0
  bool doFillComplete=true;
  bool optimizeStorage=false;
  RCP<Operator> DinvADinvAP0 = Utils::TwoMatrixMultiply(A,false,DinvAPtent,false,doFillComplete,optimizeStorage);
  Utils::MyOldScaleMatrix(DinvADinvAP0,diagA,true,doFillComplete,optimizeStorage); //scale matrix with reciprocal of diag

  Numerator = MultiplyAll(Ptent, DinvADinvAP0, GID2localgid);
  RCP<Teuchos::Array<Scalar> > Numerator2= MultiplySelfAll(DinvAPtent, GID2localgid);
  TEUCHOS_TEST_FOR_EXCEPTION(Numerator->size() != Numerator2->size(), Exceptions::RuntimeError, "PgPFactory::ComputeRowBasedOmegas: size of Numerator and Numerator2 different. Error");
  for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
    (*Numerator)[i] += (*Numerator2)[i];
  Denominator = MultiplyAll(DinvAPtent,DinvADinvAP0, GID2localgid);
  for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
    (*Denominator)[i] *= 2.;

}
break;
}

/////////////////// DEBUG: check for zeros in denominator
size_t zeros_in_denominator = 0;
for(size_t i=0; i<Teuchos::as<size_t>(Denominator->size()); i++)
{
  if((*Denominator)[i] == Teuchos::ScalarTraits<Scalar>::zero()) zeros_in_denominator ++;
}
if(zeros_in_denominator>Teuchos::ScalarTraits<Scalar>::zero())
  GetOStream(Warnings0, 0) << "Found " << zeros_in_denominator<< " zeros in Denominator. very suspicious!" << std::endl;

/////////////////// build ColBasedOmegas
RCP<Teuchos::ArrayRCP<Scalar> > ColBasedOmegas = Teuchos::rcp(new Teuchos::ArrayRCP<Scalar>(Numerator->size(),Teuchos::ScalarTraits<Scalar>::zero()));
for(size_t i=0; i<Teuchos::as<size_t>(Numerator->size()); i++)
{
  (*ColBasedOmegas)[i] = (*Numerator)[i]/(*Denominator)[i];
  if((*ColBasedOmegas)[i] < Teuchos::ScalarTraits<Scalar>::zero())
    (*ColBasedOmegas)[i] = Teuchos::ScalarTraits<Scalar>::zero();
}
#endif

} //namespace MueLu

#endif /* MUELU_PGPFACTORY_DEF_HPP_ */
