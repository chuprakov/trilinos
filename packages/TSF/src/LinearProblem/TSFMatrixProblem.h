#ifndef TSFMATRIXPROBLEM_H
#define TSFMATRIXPROBLEM_H

#include "TSFDefs.h"
#include "TSFLinearProblemBase.h"
#include "TSFMatrixOperator.h"
#include "TSFVector.h"


namespace TSF
{


  /** \ingroup CoreSubtypes
   * TSFMatrixProblem
   *
   */

  class TSFMatrixProblem : public TSFLinearProblemBase
    {
    public:
      /** empty ctor only */
      TSFMatrixProblem(TSFMatrixOperator* matrix);

      /** virtual dtor */
      virtual ~TSFMatrixProblem(){;}

      /** returns the right-hand side vector b */
      virtual TSFVector getRHS() const ;

      /** returns a RHS vector of a type specified by the input
       * vector space */
      virtual TSFVector getRHS(const TSFVectorSpace& space) const ;

      /** By default, the solution is not known, so the base class
       * implementation throws an exception. */
      virtual TSFVector getKnownSolution() const ;

      /** By default, the solution is not known, so the base class
       * implementation throws an exception. */
      virtual TSFVector getKnownSolution(const TSFVectorSpace& space) const ;

      /** returns the linear operator A */
      virtual TSFLinearOperator getOperator() const ;

    protected:
      /** */
      virtual int nGlobalRows() const = 0 ;

      /** */
      virtual int nGlobalColumns() const {return nGlobalRows();}

      /** */
      virtual int nLocalRows() const = 0 ;

      /** */
      virtual int lowestLocalRow() const = 0 ;

      /** fill the matrix with values */
      virtual void fillMatrix() const = 0 ;

      /** fill the RHS vector with values */
      virtual void fillRHS(TSFVector& rhs) const = 0 ;

      /** fill the known solution vector with values. The base class method
       * throws an exception, since in general we will not know the soln */
      virtual void fillKnownSolution(TSFVector& soln) const ;


      /** create the update list */
      virtual TSFSmartPtr<TSFArray<int> > formUpdateList() const = 0 ;


      /** create the update list */
      virtual TSFSmartPtr<TSFArray<int> > formBandwidth() const = 0 ;

      void buildMatrixStructure() const ;

      /* the operator */
      mutable TSFLinearOperator op_;
      /* the matrix */
      mutable TSFMatrixOperator* matrix_;

      /* */
      mutable TSFSmartPtr<TSFArray<int> > localRowIndices_;

      /* */
      mutable TSFSmartPtr<TSFArray<TSFNonDupArray<int> > > graph_;

      /* */
      mutable bool matrixReady_;


    };
}

#endif
