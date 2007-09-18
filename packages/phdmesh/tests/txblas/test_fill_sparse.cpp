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

#include <iostream>
#include <algorithm>

#include <util/ParallelComm.hpp>

#include <txblas/SparseMatrix.hpp>

using namespace phdmesh ;

void test_fill_sparse_band(
  ParallelMachine comm ,
  const std::vector<unsigned> & partition ,
  const unsigned iband ,
  const unsigned nband ,
  const unsigned stride ,
  const double evalue ,
  const unsigned block_size ,
  SparseMatrix & mat )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  const unsigned local_irow = partition[ p_rank ];
  const unsigned local_nrow = partition[ p_rank + 1 ] - local_irow ;

  const unsigned nglobal = partition[ p_size ];
  const unsigned nzrow = 1 + 2 * nband ;
  const unsigned nztotal = nzrow * local_nrow ;

  mat.allocate( comm , partition , nztotal , block_size );

  txblas_SparseMatrixEntry * const coeff = mat.matrix();

  for ( unsigned i = 0 ; i < nztotal ; ++i ) { coeff[i].val = -1 ; }

  for ( unsigned i = 0 ; i < local_nrow ; ++i ) {
    const unsigned irow = local_irow + i ;

    unsigned k = i * nzrow ;

    coeff[k].row = i ;
    coeff[k].col = irow ;
    coeff[k].val = evalue + ( nzrow - 1 );
    ++k ;

    for ( unsigned j = 0 ; j < nband ; ++j ) {
      const unsigned b = iband + j * stride ;
      coeff[k].row = i ;
      coeff[k].col = ( irow + b ) % nglobal ;
      ++k ;

      coeff[k].row = i ;
      coeff[k].col = ( irow + nglobal - b ) % nglobal ;
      ++k ;
    }
  }

  mat.commit();
}


