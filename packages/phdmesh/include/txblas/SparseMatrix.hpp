/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2007) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/*  This library is free software; you can redistribute it and/or modify  */
/*  it under the terms of the GNU Lesser General Public License as        */
/*  published by the Free Software Foundation; either version 2.1 of the  */
/*  License, or (at your option) any later version.                       */
/*                                                                        */
/*  This library is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     */
/*  Lesser General Public License for more details.                       */
/*                                                                        */
/*  You should have received a copy of the GNU Lesser General Public      */
/*  License along with this library; if not, write to the Free Software   */
/*  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307   */
/*  USA                                                                   */
/*------------------------------------------------------------------------*/
/**
 * @author H. Carter Edwards
 * @date   August 2007
 */

#ifndef txblas_SparseMatrix_hpp
#define txblas_SparseMatrix_hpp

#include <vector>
#include <txblas/sparse_mxv.h>
#include <util/Parallel.hpp>

namespace phdmesh {

/** Simple partitioning of an array of 'nglobal' values.
 *    partition[ p_rank ]     = begin ordinal
 *    partition[ p_rank + 1 ] = end   ordinal
 */
void simple_partition(
  ParallelMachine         comm ,
  const unsigned          nglobal ,
  std::vector<unsigned> & partition );

class SparseMatrix {
private:
  SparseMatrix( const SparseMatrix & );
  SparseMatrix & operator = ( const SparseMatrix & );

  // Global:
  ParallelMachine       m_comm ;
  unsigned              m_comm_size ;
  unsigned              m_comm_rank ;
  bool                  m_sparse ;
  std::vector<unsigned> m_partition ;
  std::vector<int>      m_work_disp ;
  std::vector<int>      m_send_disp ;
  std::vector<int>      m_send_map ;

  // Local:
  unsigned              m_row_first ;
  unsigned              m_row_size ;
  unsigned              m_block_size ;
  std::vector<unsigned> m_blocking ;
  std::vector<txblas_SparseMatrixEntry> m_matrix ;

public:

  ~SparseMatrix();
  SparseMatrix();

  void allocate( ParallelMachine , 
                 const std::vector<unsigned> & partition ,
                 const unsigned local_non_zeros ,
                 const unsigned local_block_size );

  const unsigned * partition() const { return & m_partition[0] ; }

  unsigned non_zeros() const  { return m_matrix.size(); }
  unsigned row_size() const   { return m_row_size ; }
  unsigned block_size() const { return m_block_size ; }
  txblas_SparseMatrixEntry * matrix() { return & m_matrix[0] ; }

  void commit();

  void multiply( const double * const x , double * const y ) const ;
};


}

#endif

