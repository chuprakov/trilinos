/* @HEADER@ */
/* ***********************************************************************
// 
//           TSFExtended: Trilinos Solver Framework Extended
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// **********************************************************************/
/* @HEADER@ */

#ifndef TSFLOADABLEMATRIX_HPP
#define TSFLOADABLEMATRIX_HPP

#include "TSFConfigDefs.hpp"
#include <set>

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

    /** */
    virtual void configure(int lowestRow,
                           const std::vector<std::set<int> >& nonzeros);

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


    /** 
     *
     */
    virtual void addElementBatch(int numRows, 
                                 int rowBlockSize,
                                 const int* globalRowIndices,
                                 int numColumnsPerRow,
                                 const int* globalColumnIndices,
                                 const Scalar* values,
                                 const int* skipRow);
                                 

    /** Set all elements to zero, preserving the existing structure */
    virtual void zero() = 0 ;
    

    /** Finalize values of the matrix. This is a hook for any
     * implementation-dependent steps that must be done after
     * loading of elements. */
    virtual void freezeValues() = 0 ;

  private:
    
    
  };


  /* Default implementation of addElementBatch */
  template <class Scalar>
  void LoadableMatrix<Scalar>::addElementBatch(int numRows, 
                                               int rowBlockSize,
                                               const int* globalRowIndices,
                                               int numColumnsPerRow,
                                               const int* globalColumnIndices,
                                               const Scalar* values,
                                               const int* skipRow)
  {
    int numRowBlocks = numRows/rowBlockSize;
    int row = 0;

    for (int rb=0; rb<numRowBlocks; rb++)
      {
        const int* cols = globalColumnIndices + rb*numColumnsPerRow;
        for (int r=0; r<rowBlockSize; r++, row++)
          {
            if (skipRow[row]) continue;
            const double* rowVals = values + row*numColumnsPerRow;
            addToRow(globalRowIndices[row], numColumnsPerRow,
                     cols, rowVals);
          }
      }
  }

  /* Default implementation of configure */
  template <class Scalar>
  void LoadableMatrix<Scalar>::configure(int lowestRow,
                                         const std::vector<std::set<int> >& nonzeros)
  {
    std::vector<double> zeros;
    std::vector<int> colIndices;
    int maxSize = 0;

    for (int i=0; i<nonzeros.size(); i++)
      {
        std::set<int>::const_iterator iter;
        const std::set<int>& s = nonzeros[i];
        colIndices.resize(0);
        for (iter=s.begin(); iter != s.end(); iter++) 
          {
            colIndices.push_back(*iter);
          }
        if (colIndices.size() > maxSize) 
          {
            zeros.resize(colIndices.size());
            for (int j=maxSize; j<zeros.size(); j++) zeros[j] = 0.0;
            maxSize = zeros.size();
          }
        setRowValues(lowestRow + i, colIndices.size(),
                     &(colIndices[0]), &(zeros[0]));
      }
  }
}

#endif
