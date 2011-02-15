#ifndef MUELU_TENTATIVEPFACTORY_HPP
#define MUELU_TENTATIVEPFACTORY_HPP

#include <iostream>
#include "Cthulhu_Map.hpp"
#include "Cthulhu.hpp"
#include "MueLu_OperatorFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_PFactory.hpp"
#include "Cthulhu_MultiVectorFactory.hpp"
#include "Teuchos_ScalarTraits.hpp"

//FIXME this is a complete hack ... I am not proud
#include "./ml_lapack.h"

//FIXME these should be added to the Teuchos BLAS wrappers
#if defined(INTEL_CXML)
#  define BLAS_PREFIX __stdcall
#else
#  define BLAS_PREFIX
#endif
#define DGEQRF_F77  F77_BLAS_MANGLE(dgeqrf,DGEQRF)
#define DORGQR_F77  F77_BLAS_MANGLE(dorgqr,DORGQR)

void BLAS_PREFIX DGEQRF_F77(int *, int *, double *, int *,
                 double *, double *, int *, int *);

void BLAS_PREFIX DORGQR_F77(int *m, int *n, int *k, double * a,
                 int *lda, double *tau, double *work, int *lwork,
                 int *info);

namespace MueLu {

  /*!
    @class TentativePFactory class.
    @brief Factory for building tentative prolongator.
  */

template<class ScalarType, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
class TentativePFactory : public PFactory<ScalarType,LocalOrdinal,GlobalOrdinal,Node, LocalMatOps> {
#include "MueLu_UseShortNames.hpp"

  private:

    /*
    TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
    coalesceFact_;
    aggregationFact_;
    TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO 
    */
    //! use QR decomposition for improving nullspace information per default
    bool QR_;

  public:
    //! @name Constructors/Destructors.
    //@{

    //! Constructor.
    TentativePFactory(/*TODO*/ /*CoalesceFact, AggregationFact*/ /*TODO*/) : QR_(false) {
      //Teuchos::OSTab tab(this->out_);
      //MueLu_cout(Teuchos::VERB_HIGH) << "TentativePFactory: Instantiating a new factory" << std::endl;
      /*
        if (CoalesceFact != Teuchos::null) coalesceFact_ = CoalesceFact;
        else                               coalesceFact_ = rcp(new CoalesceDropFactory());
        if (AggregationFact != Teuchos::null) aggregationFact_ = AggregationFact;
        else                                  aggregationFact_ = rcp(new AggregationFactory());
      */
    }

    //! Destructor.
    virtual ~TentativePFactory() {}
    //@}

    //! @name Set/Get methods.
    //@{

    void TentativeWithQR(bool value) {
      QR_ = value;
    }

    bool TentativeWithQR() {
      return QR_;
    }

    //@}

    //! @name Build methods.
    //@{
    bool BuildP(Level & fineLevel, Level & coarseLevel) const {

      //TODO get or generate fine grid nullspace here
      RCP<MultiVector> fineNullspace;
      if (fineLevel.IsSaved("Nullspace"))
        fineLevel.CheckOut("Nullspace",fineNullspace);
      else {
        //TODO add this functionality
        //throw(Exceptions::NotImplemented("TenativePFactory.BuildP():  nullspace generation not implemented
        //yet"));
        std::cout << "nullspace generation not implemented yet" << std::endl;
      }
/*
       //TODO build aggregates
            % 1) build aggregates
            AggInfo = TentativePFactory.BuildAggregates(FineLevel, CoarseLevel,...
                this.CoalesceFact_, this.AggFact_, ...
                this.GetOutputLevel(), Specs, ...
                this.ReUseAggregates(), this.ReUseGraph());
*/
      RCP<Operator> Ptent = MakeTentative(fineLevel);
      if (coarseLevel.IsRequested("Ptent"))
        coarseLevel.Save("Ptent",Ptent);
      coarseLevel.SetP(Ptent);
      //coarseLevel.Save("nullspace",cnull);

      return true;
    }
    //@}

    //! @name Static methods.
    //@{
    /*! @brief Make tenative prolongator.
        TODO this signature does not match MueMat
    */
    static RCP<Operator> MakeTentative(Level const &currentLevel)
    {
      Teuchos::RCP< Operator > Op = currentLevel.GetA();
      GO nFineDofs = Op->getGlobalNumRows();
      GO nCoarseDofs = nFineDofs/3;
      if (nCoarseDofs*3 != nFineDofs)
        throw(Exceptions::NotImplemented("MakeTentative: currently #fine DOFS must be a multiple of 3"));
      Teuchos::RCP< Operator > Ptent = Teuchos::rcp( new CrsOperator(Op->getRowMap(), 2) );
      std::vector<SC> Values(1);
      Values[0] = 1.0/sqrt(3.0);
      std::vector<LO> Indices(1);
      Teuchos::ArrayView<SC> av(&Values[0],1);
      Teuchos::ArrayView<GO> iv(&Indices[0],1);
      for (int i=0; i<nFineDofs; i++) {
        Indices[0] = i / 3;
        Ptent->insertGlobalValues(i,iv,av);
      }
      //TODO replace this with a factory or something else
#ifdef HAVE_CTHULHU_EPETRA
      RCP<Map> domainMap = rcp( new Cthulhu::EpetraMap(nCoarseDofs,Op->getRowMap()->getIndexBase(),Op->getRowMap()->getComm()) );

      Ptent->fillComplete(domainMap, Op->getRowMap());
#else
#warning
#endif

      //MatrixPrint(Op);
      return Ptent;
    } //MakeTentative()

    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType Magnitude;

    /*! @brief Make tentative prolongator with QR.
        FIXME once completed, this should replace MakeTentative
        FIXME I make no attempt to detect if the aggregate is too small to support the NS
        FIXME This needs to use Jeremie's Aggregate class.
        FIXME Aggregates of size 3 are hard-wired in. (AFAIK, no assumptions of size 3, however.)
    */
    //FIXME return Operator instead of void
    //static RCP<Operator> MakeTentativeWithQR(Level &fineLevel, Level &coarseLevel)
    static void MakeTentativeWithQR(Level &fineLevel, Level &coarseLevel)
    {
      using Teuchos::ArrayRCP;

      //TODO checkout aggregates from Level hash table

      Teuchos::RCP< Operator > fineA = fineLevel.GetA();
      GO nFineDofs = fineA->getGlobalNumRows();
      if ((nFineDofs/3)*3 != nFineDofs)
        throw(Exceptions::NotImplemented("MakeTentative: currently #fine DOFS must be a multiple of 3"));

      GO numAggs = nFineDofs / 3;  //FIXME  should come from aggregate class:
                                   //FIXME  numAggs = aggregates->GetNumAggregates();

      //get the fine grid nullspace
      RCP<MultiVector> fineNullspace;
      if (fineLevel.IsSaved("Nullspace"))
        fineLevel.CheckOut("Nullspace",fineNullspace);
      else {
        throw(Exceptions::NotImplemented("MakeTentativeWithQR:  nullspace generation not implemented yet"));
      }
      const size_t NSDim = fineNullspace->getNumVectors();
      GO nCoarseDofs = numAggs*NSDim;
      // Create array of aggregate sizes.
      // TODO Should this come from aggregate class?
      ArrayRCP<GO> aggSizes(numAggs);
      // FIXME Right now, we hardwire in size of 3 for every aggregate.
      for (typename Teuchos::ArrayRCP<GO>::iterator i=aggSizes.begin(); i!=aggSizes.end(); ++i)
        *i = 3;

      // Find the largest aggregate.
      // TODO Should this come from aggregate class?
      LO maxAggSize=0;
      for (typename Teuchos::ArrayRCP<GO>::iterator i=aggSizes.begin(); i!=aggSizes.end(); ++i)
        if (*i > maxAggSize) maxAggSize = *i;

      // Create a lookup table to determine the rows (fine DOFs) that belong to a given aggregate.
      // aggToRowMap[i][j] is the jth DOF in aggregate i
      // TODO Question:  should the aggregate class provide this?
      ArrayRCP< ArrayRCP<GO> > aggToRowMap(numAggs);
      ArrayRCP< GO > numDofs(numAggs,0);  //Temporary counter that tracks how many DOFS have been recorded
                                          //for each each aggregate in aggToRowMap.  This is only used
                                          //while aggToRowMap is being populated.
      GO t=0;
      for (typename ArrayRCP<ArrayRCP<GO> >::iterator a2r=aggToRowMap.begin(); a2r!=aggToRowMap.end(); ++a2r)
        *a2r = ArrayRCP<GO>(aggSizes[t++]);
      for (GO i=0; i<nFineDofs; ++i) {
        GO myAgg = i/3; //TODO use aggregate class
        aggToRowMap[myAgg][ numDofs[myAgg] ] = i;
        ++(numDofs[myAgg]);
      }
      // Allocate workspace for LAPACK QR routines.
      ArrayRCP<SC> localQR(maxAggSize*NSDim); // The submatrix of the nullspace to be orthogonalized.
      LO           workSize = NSDim;          // Length of work. Must be at least dimension of nullspace.
                                              // LAPACK may calculate better value, returned in work[0].
      ArrayRCP<SC> work(workSize);            // (in/out) work vector
                                              // FIXME DGEQRF documentation says this should be min(M,N) where B=MxN
      ArrayRCP<SC> tau(NSDim);                // (out) scalar factors of elementary reflectors, input to DORGQR
      LO           info=0;                      // (out) =0: success; =i, i<0: i-th argument has illegal value

      // Pull out the nullspace vectors so that we can have random access.
      // (Question -- do we have to do this?)
      ArrayRCP< ArrayRCP<const SC> > fineNS(NSDim);
      for (size_t i=0; i<NSDim; ++i)
        fineNS[i] = fineNullspace->getData(i);

      //Allocate storage for the coarse nullspace.
      LO indexBase=fineA->getRowMap()->getIndexBase();
      //FIXME We shouldn't rely explicitly on EpetraMap
      RCP<const Map > coarseMap = rcp( new Cthulhu::EpetraMap(Teuchos::OrdinalTraits<Cthulhu::global_size_t>::invalid(),
                                                              nCoarseDofs,indexBase,fineA->getRowMap()->getComm())       );

      RCP<MultiVector> coarseNullspace = MultiVectorFactory::Build(coarseMap,NSDim);
      ArrayRCP< ArrayRCP<SC> > coarseNS(NSDim);
      for (size_t i=0; i<NSDim; ++i)
        coarseNS[i] = coarseNullspace->getDataNonConst(i);

      //Allocate temporary storage for the tentative prolongator.
      // TODO Right now, this is stored as a point matrix.
      // TODO we should be able to allocate something other than a CRS matrix
      ArrayRCP<GO> rowPtr(nFineDofs+1);
      for (GO i=0; i<=nFineDofs; ++i)
        rowPtr[i] = i*NSDim;
      ArrayRCP<GO> colPtr(nFineDofs*NSDim,0);
      ArrayRCP<SC> valPtr(nFineDofs*NSDim,0.);

      RCP<Operator> Ptentative = rcp(new CrsOperator(fineA->getDomainMap(), NSDim));

      //used in the case of just one nullspace vector
      Teuchos::Array<Magnitude> norms(1);

      //TODO Is it more efficient to extract the nullspace in aggregate order all at once,
      //TODO instead of one aggregate at a time?
      //*****************************************************************
      //Loop over all aggregates and calculate local QR decompositions.
      //*****************************************************************
      for (LO agg=0; agg<numAggs; ++agg)
      {
        LO myAggSize = aggSizes[agg];

        // For each aggregate, extract the corresponding piece of the nullspace and put it in the flat array,
        // "localQR" (in column major format) for the QR routine.
           for (size_t j=0; j<NSDim; ++j) {
             for (LO k=0; k<myAggSize; ++k) {
                //aggToRowMap[agg][k] is the kth DOF in the ith aggregate
                //fineNS[j][n] is the nth entry in the jth NS vector
                localQR[j* myAggSize + k] = fineNS[j][ aggToRowMap[agg][k] ];
             } //for (LO k=0 ...
           } //for (LO j=0 ...

        int intFineNSDim = NSDim;

        if (NSDim == 1) {
          //only one nullspace vector, so normalize by hand
          fineNullspace->norm2(norms);
          tau[0] = localQR[0];
          localQR[0] = norms[0];
        } else {
          //Perform the QR.  Upon return, R is stored explicitly, Q is stored implicitly as product
          //of reflection matrices.
          DGEQRF_F77( &myAggSize, &intFineNSDim, localQR.getRawPtr(), &myAggSize,
                      tau.getRawPtr(), work.getRawPtr(), &workSize, &info );
        }

        if (info != 0) {
          std::ostringstream buf;
          buf << info;
          std::string msg = "MakeTentativeWithQR: dgeqrf (LAPACK QR routine) returned error code " + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }

        // LAPACK may have determined a better length for the work array.  Returns it in work[0],
        // so we cast to avoid compiler warnings.
        if ( work[0] > workSize) {
          workSize = (int) work[0];
          work = ArrayRCP<SC>(workSize);
        } else
          workSize = (int) work[0]; //TODO huh, think about this -- should it ever shrink?

        // Extract R, the coarse nullspace.  This is stored in upper triangular part of localQR.
        // Note:  coarseNS[i][.] is the ith coarse nullspace vector, which may be counter to your intuition.
        Cthulhu::global_size_t offset=agg*NSDim;
        for (size_t j=0; j<NSDim; ++j) {
          for (size_t k=0; k<=j; ++k) {
            coarseNS[j][offset+k] = localQR[ myAggSize*j + k ];
          }
        }

        // Calculate Q, the tentative prolongator, explicitly.  This requires calling a second LAPACK routine.

        if (NSDim == 1) {
          //again, only one nullspace vector, so calculate Q by hand
          norms[0] = localQR[0];
          localQR[0] = tau[0];
          norms[0] = 1. / norms[0];
          for (LO i=0; i<myAggSize; ++i)
            localQR[i] *= norms[0];
        } else {
          DORGQR_F77( &myAggSize, &intFineNSDim, &intFineNSDim, localQR.getRawPtr(),
                      &myAggSize, tau.getRawPtr(), work.getRawPtr(), &workSize, &info );
        }

        if (info != 0) {
          std::ostringstream buf;
          buf << info;
          std::string msg = "MakeTentativeWithQR: dorgqr (LAPACK auxiliary QR routine) returned error code "
                             + buf.str();
          throw(Exceptions::RuntimeError(msg));
        }

        // LAPACK may have determined a better length for the work array.  Returns it in work[0],
        // so we cast to avoid compiler warnings.
        if ( work[0] > workSize) {
          workSize = (int) work[0];
          work = ArrayRCP<SC>(workSize);
        } else
          workSize = (int) work[0]; //TODO huh, think about this -- should it ever shrink?

        //Save the part of tentative P (the Q factor) corresponding to the current aggregate
        //in a CSR format.  This saves the entire Q factor as if it were dense.
        GO index;
        for (GO j=0; j<myAggSize; ++j) {
          for (size_t k=0; k<NSDim; ++k) {
            index = rowPtr[aggToRowMap[agg][j]]+k;
            colPtr[index] = agg * NSDim + k;
            valPtr[index] = localQR[k*myAggSize+j];
          }
        }

      } // for (LO agg=0; agg<numAggs; ++agg)

      // ***********************************************************
      // ************* end of aggregate-wise QR ********************
      // ***********************************************************

      //Remove all zeros in the tentative prolongator data arrays.
      GO k = rowPtr[0];
      GO nNonzeros=0;
      GO index = k;
      for (GO i=0; i<nFineDofs; ++i)
      {
        for (GO j=k; j< rowPtr[i+1]; ++j) {
          if (valPtr[j] != 0.0) {
            valPtr[index] = valPtr[j];
            colPtr[index] = colPtr[j];
            ++index;
            ++nNonzeros;
          }
        }
        k = rowPtr[i+1];
        rowPtr[i+1] = index;
      } //for (GO i=0...

      //Insert into the final tentative prolongator matrix.
      //TODO I don't know how to efficiently insert data into a matrix.
      //TODO In Epetra, you can use Views of data pointers.  In Tpetra, I don't know whether this can be done.
      //TODO So for right now, I just insert row by row.
      RCP<const Map> fineRowMap = fineA->getRowMap();
      for (GO i=0; i<nFineDofs; ++i)
      {
        //j=rowPtr[i], j1 = rowPtr[i+1]
        //columns are in colPtr[j] ... colPtr[j1]-1
        //values are in valPtr[j] ... valPtr[j1]-1
        GO start = rowPtr[i];
        GO nnz = rowPtr[i+1] - rowPtr[i];
        Ptentative->insertGlobalValues(fineRowMap->getGlobalElement(i),
                                      colPtr.view(start,nnz),
                                      valPtr.view(start,nnz));
      }

      Ptentative->fillComplete(coarseMap,fineA->getRangeMap());

      coarseLevel.Save("Nullspace",coarseNullspace);
      coarseLevel.Save("Ptent",Ptentative);

    } //MakeTentativeWithQR()

    static void BuildAggregates() {
      throw(Exceptions::NotImplemented("TentativePFactory: BuildAggregates not implemented"));
    }
    static void MakeNoQRTentative() {
      throw(Exceptions::NotImplemented("TentativePFactory: MakeNoQRTentative not implemented"));
    }
    //@}

}; //class TentativePFactory

} //namespace MueLu

#define MUELU_TENTATIVEPFACTORY_SHORT

#endif //ifndef MUELU_TENTATIVEPFACTORY_HPP
