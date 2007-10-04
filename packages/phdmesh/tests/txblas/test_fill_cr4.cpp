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

#include <txblas/CR4Matrix.hpp>

using namespace phdmesh ;

void test_fill_cr4_band(
  ParallelMachine comm ,
  const std::vector<unsigned> & partition ,
  const unsigned iband ,
  const unsigned nband ,
  const unsigned stride ,
  const double evalue ,
  std::vector<unsigned> & prefix ,
  std::vector<txblas_cr4> & matrix )
{
  const unsigned p_size = parallel_machine_size( comm );
  const unsigned p_rank = parallel_machine_rank( comm );

  const unsigned local_irow = partition[ p_rank ];
  const unsigned local_nrow = partition[ p_rank + 1 ] - local_irow ;

  const unsigned nglobal = partition[ p_size ];
  const unsigned nzrow = 1 + 2 * nband ;
  const unsigned nzrow_alloc = ( nzrow / 4 ) + ( nzrow % 4 ? 1 : 0 );
  const unsigned local_alloc = local_nrow * nzrow_alloc ;

  prefix.resize( local_nrow + 1 );
  matrix.resize( local_alloc );

  for ( unsigned i = 0 ; i <= local_nrow ; ++i ) {
    prefix[i] = i * nzrow_alloc ;
  }

  std::vector< std::pair<unsigned,double> > row( nzrow );

  for ( unsigned i = 0 ; i < local_nrow ; ++i ) {
    const unsigned irow = local_irow + i ;

    row[0].first  = irow ;
    row[0].second = evalue + ( nzrow - 1 );

    for ( unsigned j = 0 ; j < nband ; ++j ) {
      const unsigned b = iband + j * stride ;
      row[j+1].first  = ( irow + b ) % nglobal ;
      row[j+1].second = -1 ;

      row[j+1+nband].first = ( irow + nglobal - b ) % nglobal ;
      row[j+1+nband].second = -1 ;
    }

    std::sort( row.begin() , row.end() );

    const unsigned k = i * nzrow_alloc ;

    unsigned j = 0 ;
    for ( ; j < nzrow ; ++j ) {
      const unsigned jd4 = j / 4 ;
      const unsigned jm4 = j % 4 ;
      matrix[k+jd4].col[jm4] = row[j].first ;
      matrix[k+jd4].val[jm4] = row[j].second ;
    }
    for ( ; j < nzrow_alloc * 4 ; ++j ) {
      const unsigned jd4 = j / 4 ;
      const unsigned jm4 = j % 4 ;
      matrix[k+jd4].col[jm4] = matrix[k+jd4].col[0] ;
      matrix[k+jd4].val[jm4] = 0 ;
    }
  }
}


