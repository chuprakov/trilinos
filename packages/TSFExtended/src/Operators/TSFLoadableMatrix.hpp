/* @HEADER@ */
/* @HEADER@ */

#ifndef TSFLOADABLEMATRIX_HPP
#define TSFLOADABLEMATRIX_HPP

#include "TSFConfigDefs.hpp"

namespace TSFExtended
{
  /** 
   * Class LoadableMatrix provides an abstract interface for configuration
   * and loading of matrices. 
   */
  template <class Scalar>
  class LoadableMatrix 
  {
  public:
    /** Virtual dtor */
    virtual ~LoadableMatrix(){;}

    /** 
     * Set the locations of all my nonzero elements. 
     * @param nLocalRows number of locally-owned rows
     * @param globalRowIndex array of global indices of the local rows
     * @param numNonzeros array of number of nonzeros for each row
     * @param array of arrays of column indices for each row
     */
    virtual void setGraph(int nLocalRows,
                          const int* globalRowIndex,
                          const int* numNonzeros,
                          const int** columnIndices) = 0 ;
    
    

    /** Insert a set of elements in a row, overwriting any previously
     * existing values. 
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     * @param elements array of element values. Must be nElemsToInsert in
     * length
     */
    virtual void setRowValues(int globalRowIndex,
                              int nElemsToInsert,
                              const int* globalColumnIndices,
                              const Scalar* elementValues) = 0 ;

    /** Insert a set of elements in a row, adding to any previously
     * existing values. 
     * @param globalRowIndex the global index of the row to which these
     * elements belong.
     * @param nElemsToInsert the number of elements being inserted in this
     * step
     * @param globalColumnIndices array of column indices. Must 
     * be nElemsToInsert in length. 
     * @param elements array of element values. Must be nElemsToInsert in
     * length
     */
    virtual void addToRow(int globalRowIndex,
                          int nElemsToInsert,
                          const int* globalColumnIndices,
                          const Scalar* elementValues) = 0 ;

    /** Set all elements to zero, preserving the existing structure */
    virtual void zero() = 0 ;
    

    /** Finalize values of the matrix. This is a hook for any
     * implementation-dependent steps that must be done after
     * loading of elements. */
    virtual void freezeValues() = 0 ;

  private:
    
    
  };
}

#endif
