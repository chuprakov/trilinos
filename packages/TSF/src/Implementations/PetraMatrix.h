#ifndef PETRAMATRIX_H
#define PETRAMATRIX_H



#include "TSFDefs.h"
#include "TSFNonDupArray.h"
#include "TSFSmartPtr.h"
#include "TSFMatrixOperator.h"
#include "TSFPreconditioner.h"
#include "TSFArray.h"



#define PETRA_BOOL_SUPPORTED

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"

namespace TSF
{
  using std::string;

  /** \ingroup Petra
   * Wrapper for Petra matrices
   */

  class PetraMatrix : public TSFMatrixOperator
    {
    public:
      /** construct with domain and range space. These must be Petra vector
       * spaces, or else an error will be thrown */
      PetraMatrix(const TSFVectorSpace& domain,
                  const TSFVectorSpace& range);

      /** the usual virtual dtor */
      virtual ~PetraMatrix(){;}

      /** apply operator to a vector in the domain space, returning a vector
       * in the range space */
      virtual void apply(const TSFVector& in,
                         TSFVector& out) const ;

      /** apply adjoint operator to a vector in the range space, returning
       * a vector in the domain space. The default implementation throws an
       * exception */
      virtual void applyAdjoint(const TSFVector& in,
                                TSFVector& out) const ;

      /** \name matrix configuration interface */
      //@{

      /** inform caller if a graph is needed to configure this matrix */
      virtual bool requiresGraph() const {return false;}

      /** inform caller if row bandwidth is needed to configure this matrix */
      virtual bool requiresBandwidth() const {return true;}

      /** inform caller if this matrix type can handle non-square matrices */
      virtual bool supportsNonSquare() const {return true;}

      /** set the columns to be used in a given row initialize to zero*/
      virtual void setRowStructure(int globalRowIndex, int bandwidth,
                                   const int* columnIndices);
      /** set the columns to be used in a given row initialized
          to values given*/
      virtual void setRowStructure(int globalRowIndex, int bandwidth,
                                   const int* columnIndices,
                                   const double* values);

      /** set the bandwith of all rows */
      virtual void setBandwidth(int nLocalRows, const int* bandwidth) ;

      /** finalize values. This makes a call to Petra's FillComplete(). */
      virtual void freezeValues();

      /** finalize structure. The Petra matrix object is actually constructed
       * at this point */
      virtual void freezeStructure();
      //@}

      /** \name matrix loading interface */
      //@{

      /** set an element  */
      virtual void setElement(int i, int j, const TSFReal& aij);


      /** add to selected elements of a row in the matrix */
      virtual void addToRow(int globalRowIndex,
                            int nCols,
                            const int* globalColumnIndices,
                            const TSFReal* a);
      virtual void getRow(int row, TSFArray<int>& indices,
                          TSFArray<TSFReal>& values) const;

      virtual TSFLinearOperator* getTranspose();

      /** set all elements to zero */
      virtual void zero();

      //@}

      virtual void print(ostream& os) const ;

      virtual TSFLinearOperatorBase* deepCopy() const ;

      /** \name incomplete factorization preconditioning interface */
      //@{

      /** create a k-level incomplete factorization. Default is to throw
       * an error. */
      virtual void getILUKPreconditioner(int fillLevels,
                                         int overlapFill,
                                         TSFPreconditioner& rtn) const ;

      //@}

      // VEH
      /** create a k-level incomplete factorization for a right preconditioner.
       * Default is to throw
       * an error. */
      virtual void getILUKRightPreconditioner(int fillLevels,
                                              int overlapFill,
                                              TSFPreconditioner& rtn) const ;

      //@}


      /** append my timings to a list of timings */
      static void collectTimings(TSFArray<TSFTimer>& timers) ;

      static Epetra_CrsMatrix* getConcrete(const TSFLinearOperator& A);

      /* VEH/RST */
      /* Insert an epetra matrix into a TSF PetraMatrix. */
      void setPetraMatrix(Epetra_CrsMatrix* A, bool ownership)
        {
          matrix_ = TSFSmartPtr<Epetra_CrsMatrix>(A,ownership);
        }


      /** timer for matrix-vector multiplies */
      static TSFTimer& mvMultTimer();
      /** timer for forming ILU preconditioner */
      static TSFTimer& iluTimer();
    private:
      void petraCheck(int ierr, const string& methodName) const ;

      void ptrCheck(const string& methodName) const ;
      TSFSmartPtr<Epetra_Map> rowMap_;
      TSFSmartPtr<Epetra_Map> columnMap_;
      TSFSmartPtr<Epetra_CrsMatrix> matrix_;
      TSFArray<int> nCols_;
      TSFSmartPtr<Epetra_Comm> petraComm_;
      bool transposed_;
      TSFLinearOperator opTrp_;
      mutable bool hasCachedPreconditioner_;
      mutable TSFPreconditioner cachedPreconditioner_;
      /*      static TSFTimer mvMultTimer_; */
      /*      static TSFTimer ILUTimer_; */
      bool hasTranspose_;

    };


}

#endif
