#ifndef TSFMATRIXVIEW_H
#define TSFMATRIXVIEW_H

#include "TSFDefs.h"
#include "TSFSmartPtr.h"
#include "TSFMatrixOperator.h"
#include "TSFLinearOperator.h"

namespace TSF
{
  class TSFPreconditioner;


  /**
   * TSFMatrixView provides high-level access to the matrix configuration, loading,
   * and factoring methods of TSFMatrixOperator. A TSFMatrixView is constructed from
   * a TSFLinearOperator. In the constructor call, the operator's implem pointer is
   * dynamically cast to a TSFMatrixOperator pointer; is the operator is not a
   * TSFMatrixOperator this cast will fail and an error will result.
   */

  class TSFMatrixView
    {
    public:
      /** \name Constructors */
      //@{
      /** construct from an existing operator */
      TSFMatrixView(const TSFLinearOperator& op);
      //@}

      /** \name matrix configuration interface */
      //@{
      /** inform caller if a full graph is needed to configure this matrix */
      bool requiresGraph() const {return matrix_->requiresGraph();}

      /** inform caller if a bandwidth array
          is needed to configure this matrix */
      bool requiresBandwidth() const {return matrix_->requiresBandwidth();}

      /** inform caller if this matrix type can handle non-square matrices */
      bool supportsNonSquare() const {return matrix_->supportsNonSquare();}

      /** Set the full sparsity graph, giving the column indices for each
       * row owned by the current processor. */
      void setGraph(int nLocalRows, const int* bandwidth,
                    const int** columnIndices)
        {matrix_->setGraph(nLocalRows, bandwidth, columnIndices);}

      /** set the columns to be used in a given row */
      void setRowStructure(int globalRowIndex, int bandwidth,
                           const int* columnIndices)
        {matrix_->setRowStructure(globalRowIndex, bandwidth, columnIndices);}

      /** set the bandwith of each row */
      void setBandwidth(int nLocalRows, const int* bandwidth)
        {matrix_->setBandwidth(nLocalRows, bandwidth);}

      /** hook for implementation-dependent structure finalization call */
      void freezeStructure() {matrix_->freezeStructure();}

      /** hook for implementation-dependent structure finalization call */
      void freezeValues() {matrix_->freezeValues();}
      //@}

      /** \name matrix loading interface */
      //@{

      /** add to selected elements of a row in the matrix */
      void addToRow(int globalRowIndex,
                    int nCols,
                    const int* globalColumnIndices,
                    const TSFReal* a)
        {matrix_->addToRow(globalRowIndex, nCols, globalColumnIndices, a);}

      /** set a single element */
      void setElement(int i, int j, const TSFReal& aij) {matrix_->setElement(i,j,aij);}

      /** set all elements to zero */
      void zero() {matrix_->zero();}
      //@}

      /** \name incomplete factorization interface */
      //@{
      /** create a k-level incomplete factorization. */
      void getILUKPreconditioner(int fillLevels,
                                 int overlapFill,
                                 TSFPreconditioner& rtn) const
        {matrix_->getILUKPreconditioner(fillLevels, overlapFill, rtn);}
      //@}

      /** \name factoring interface */
      //@{
      /** indicate whether the matrix is now stored in factored form */
      bool isFactored() const {return matrix_->isFactored();}

      /** factor the matrix */
      void factor() {matrix_->factor();}
      //@}

    protected:
      TSFMatrixOperator* matrix_;
      TSFLinearOperator op_;
    };

}

#endif
