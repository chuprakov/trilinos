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
 */

#include <pthread.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <util/taskpool.h>
#include <txblas/sparse_mxv.h>

/*--------------------------------------------------------------------*/
/* sort: ( row / BN , col , row % BN ) */

static int txblas_rbcr_compare( const void * lhs , const void * rhs )
{
  static unsigned block_size = 0 ;

  int diff = 0 ;

  if ( NULL != lhs ) {
    if ( NULL != rhs ) {
      const txblas_SparseMatrixEntry * const l =
        (const txblas_SparseMatrixEntry *) lhs ;

      const txblas_SparseMatrixEntry * const r =
        (const txblas_SparseMatrixEntry *) rhs ;

      diff = (int) ( l->row / block_size ) - (int) ( r->row / block_size );
      if ( diff ) {
        fprintf(stderr,"txblas_rbcr_compare argument error");
        exit(-1);
      }
      if ( ! diff ) { diff = (int) l->col - (int) r->col ; }
      if ( ! diff ) { diff = (int) l->row - (int) r->row ; }
    }
    else {
      diff = block_size = *((const unsigned *) lhs );
    }
  }
  return diff ;
}

static txblas_SparseMatrixEntry * txblas_row_split(
  txblas_SparseMatrixEntry * i ,
  txblas_SparseMatrixEntry * j ,
  const unsigned row_split )
{
  txblas_SparseMatrixEntry tmp ;
  while ( i < j ) {
    while ( i < --j && row_split < j->row );

    while ( i < j && i->row < row_split ) { ++i ; }

    tmp.row = i->row ;  tmp.col = i->col ;  tmp.val = i->val ;
    i->row  = j->row ;  i->col  = j->col ;  i->val  = j->val ;
    i->row  = tmp.row ; i->col  = tmp.col ; i->val  = tmp.val ;
  }
  return i ;
}

/*--------------------------------------------------------------------*/

typedef struct txblasTask_rbcr_MatrixStruct {
  unsigned                         number_row ;
  unsigned                         block_size ;
  const unsigned                 * p_begin ;
  const txblas_SparseMatrixEntry * a_begin ;
  const double                   * x_begin ;
        double                   * y_begin ;
} txblasTask_rbcr_Matrix ;

static int txblas_task_rbcr_mxv(
  void * data , unsigned p_size , unsigned p_rank )
{
  if ( p_size <= p_rank ) return -1 ;

  /* Local 'const' copies of shared data */

  txblasTask_rbcr_Matrix * const t = (txblasTask_rbcr_Matrix*) data ;

  const unsigned p_next     = p_rank + 1 ;
  const unsigned number_row = t->number_row ;
  const unsigned block_size = t->block_size ;
  const unsigned number_block = ( number_row / block_size ) +
                                ( number_row % block_size ? 1 : 0 );
  const unsigned beg_block = ( number_block * p_rank ) / p_size ;
  const unsigned end_block = ( number_block * p_next ) / p_size ;

  const txblas_SparseMatrixEntry * const a_beg = t->a_begin ;
  const unsigned                 * const p_beg = t->p_begin ;
  const double                   * const x_beg = t->x_begin ;
        double                   * const y_beg = t->y_begin ;

  double ytmp[ block_size ];

  const txblas_SparseMatrixEntry * a = a_beg + p_beg[ beg_block ];

  for ( unsigned iblock = beg_block ; iblock < end_block ; ) {
    const unsigned row_beg = block_size * iblock ;
    const unsigned row_num = ( block_size < number_row - row_beg )
                             ? block_size : number_row - row_beg ;

    double * const y_end = ytmp + row_num ;
    double * y ;

    for ( y = ytmp ; y < y_end ; ++y ) { *y = 0 ; }

    const txblas_SparseMatrixEntry * const a_end = a_beg + p_beg[ ++iblock ];

    while ( a != a_end ) {
      ytmp[ a->row - row_beg ] += a->val * x_beg[ a->col ];
      ++a ;
    }

    double * y_out = y_beg + row_beg ;

    for ( y = ytmp ; y < y_end ; ++y , ++y_out ) { *y_out = *y ; }
  }

  return 0 ;
}

void txblas_rbcr_mxv(
  const unsigned                 number_row  /* Number rows */ ,
  const unsigned                 block_size  /* block size */ ,
  const txblas_SparseMatrixEntry a[] ,
  const unsigned                 blocking[] ,
  const double                   x[] ,  /* Input vector */
        double                   y[] )  /* Output vector */
{
  const unsigned nlock = 0 ;

  txblasTask_rbcr_Matrix data ;
  data.number_row = number_row ;
  data.block_size = block_size ;
  data.p_begin    = blocking ;
  data.a_begin    = a ;
  data.x_begin    = x ;
  data.y_begin    = y ;

  phdmesh_taskpool_run( & txblas_task_rbcr_mxv , & data , nlock );
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/

typedef struct txblasTask_MatrixPrepStruct {
  unsigned                   number_row ;
  unsigned                   block_size ;
  unsigned                   iter_level ;
  unsigned                 * p_begin ;
  txblas_SparseMatrixEntry * a_begin ;
} txblasTask_MatrixPrep ;

static int txblas_task_rbcr_prep(
  void * data , unsigned p_size , unsigned p_rank )
{
  if ( p_size <= p_rank ) return -1 ;

  txblasTask_MatrixPrep * const t = (txblasTask_MatrixPrep *) data ;

  const unsigned                   number_row = t->number_row ;
  const unsigned                   block_size = t->block_size ;
  const unsigned                   level      = t->iter_level ;
  unsigned                 * const p_begin    = t->p_begin ;
  txblas_SparseMatrixEntry * const a_begin    = t->a_begin ;

  const unsigned part_num     = 1u << ( level - 1 );
  const unsigned level_denom  = 1u << level ;
  const unsigned number_block = number_row / block_size +
                              ( number_row % block_size ? 1 : 0 );

  const size_t size = ( (unsigned char *)(a_begin+1)) -
                      ( (unsigned char *)(a_begin  ));

  const unsigned part_beg = ( part_num * p_rank ) / p_size ;
  const unsigned part_end = ( part_num * ( p_rank + 1 ) ) / p_size ;

  for ( unsigned block = part_beg ; block < part_end ; ++block ) {

    unsigned k = 2 * block ;

    const double f_begin = (double) k        / (double) level_denom ;
    const double f_split = (double)( k + 1 ) / (double) level_denom ;
    const double f_end   = (double)( k + 2 ) / (double) level_denom ;

    const unsigned iblock_begin = (unsigned)( number_block * f_begin );
    const unsigned iblock_split = (unsigned)( number_block * f_split );
    const unsigned iblock_end   = (unsigned)( number_block * f_end );

    txblas_SparseMatrixEntry * i = a_begin + p_begin[iblock_begin];
    txblas_SparseMatrixEntry * j = a_begin + p_begin[iblock_end];

    if ( iblock_begin < iblock_split && iblock_split < iblock_end ) {

      /* Splitting the span among blocks */

      const unsigned row_split = iblock_split * block_size ;

      p_begin[ iblock_split ] = txblas_row_split( i, j, row_split ) - a_begin ;
    }
    else if ( iblock_begin < iblock_end ) {

     /* Final sort of a single block */

     const size_t num = j - i ;

     qsort( i , num , size , & txblas_rbcr_compare );
    }
  }
  return 0 ;
}

int txblas_rbcr_prep(
  const unsigned number_row ,
  const unsigned block_size ,
  const unsigned number_coef ,
        txblas_SparseMatrixEntry a[] ,
        unsigned blocking[] )
{
  txblas_rbcr_compare( & block_size , NULL );

  txblasTask_MatrixPrep data ;

  data.number_row = number_row ;
  data.block_size = block_size ;
  data.p_begin    = blocking ;
  data.a_begin    = a ;

  const unsigned nlock = 0 ;

  const unsigned nblock = number_row / block_size +
                        ( number_row % block_size ? 1 : 0 );

  blocking[0] = 0 ;
  blocking[ nblock ] = number_coef ;

  for ( data.iter_level = 0 ; ( 1u << data.iter_level ) < nblock ; ) {
    ++(data.iter_level);
    phdmesh_taskpool_run( & txblas_task_rbcr_prep , & data , nlock );
  }
  ++(data.iter_level);
  phdmesh_taskpool_run( & txblas_task_rbcr_prep , & data , nlock );

  for ( unsigned i = 1 ; i < number_coef ; ++i ) {
    const unsigned i_block = a[i-1].row / block_size ;
    const unsigned j_block = a[i].row / block_size ;
    if ( j_block < i_block ) {
      abort();
    }
    else if ( j_block == i_block ) {
      if ( a[i].col < a[i-1].col ) {
        abort();
      }
      else if ( a[i].col == a[i-1].col ) {
        if ( a[i].row <= a[i-1].row ) {
          abort();
        }
      }
    }
  }

  return 0 ;
}

/*--------------------------------------------------------------------*/
/*--------------------------------------------------------------------*/
/* sort: ( row , col ) */

static int txblas_cr_compare( const void * lhs , const void * rhs )
{
  const txblas_SparseMatrixEntry * const l =
    (const txblas_SparseMatrixEntry *) lhs ;

  const txblas_SparseMatrixEntry * const r =
    (const txblas_SparseMatrixEntry *) rhs ;

  return (int) l->col - (int) r->col ;
}

static int txblas_task_cr_mxv(
  void * data , unsigned p_size , unsigned p_rank )
{
  if ( p_size <= p_rank ) return -1 ;

  /* Local 'const' copies of shared data */

  txblasTask_rbcr_Matrix * const t = (txblasTask_rbcr_Matrix*) data ;

  const unsigned p_next  = p_rank + 1 ;
  const unsigned all_row = t->number_row ;
  const unsigned beg_row = ( all_row * p_rank ) / p_size ;
  const unsigned end_row = ( all_row * p_next ) / p_size ;

  const double * const x_beg = t->x_begin ;

  const unsigned *       p     = t->p_begin + beg_row ;
  const unsigned * const p_end = t->p_begin + end_row ;
  const txblas_SparseMatrixEntry * const a_begin = t->a_begin ;

        double                   * y = t->y_begin + beg_row ;
  const txblas_SparseMatrixEntry * a = a_begin + *p ;

  while ( p < p_end ) {
    double ytmp = 0 ;

    const txblas_SparseMatrixEntry * const a_end = a_begin + *++p ;

    while ( a < a_end ) {
      ytmp += a->val * x_beg[ a->col ];
      ++a ;
    }

    *y++ = ytmp ;
  }

  return 0 ;
}


void txblas_cr_mxv(
  const unsigned                 number_row  /* Number rows */ ,
  const txblas_SparseMatrixEntry a[] ,
  const unsigned                 blocking[] ,
  const double                   x[] ,  /* Input vector */
        double                   y[] )  /* Output vector */
{
  const unsigned nlock = 0 ;

  txblasTask_rbcr_Matrix data ;
  data.number_row = number_row ;
  data.block_size = 0 ;
  data.p_begin = blocking ;
  data.a_begin = a ;
  data.x_begin = x ;
  data.y_begin = y ;

  phdmesh_taskpool_run( & txblas_task_cr_mxv , & data , nlock );
}

static int txblas_task_cr_prep(
  void * data , unsigned p_size , unsigned p_rank )
{
  if ( p_size <= p_rank ) return -1 ;

  txblasTask_MatrixPrep * const t = (txblasTask_MatrixPrep *) data ;

  const unsigned                   number_row = t->number_row ;
  const unsigned                   level      = t->iter_level ;
  unsigned                 * const p_begin    = t->p_begin ;
  txblas_SparseMatrixEntry * const a_begin    = t->a_begin ;

  const unsigned p_next   = p_rank + 1 ;
  const unsigned num_part = 1u << level ;

  const unsigned part_beg = ( num_part * p_rank ) / p_size ;
  const unsigned part_end = ( num_part * p_next ) / p_size ;

  for ( unsigned ipart = part_beg ; ipart < part_end ; ++ipart ) {
    const double ifraction = ((double) ipart ) / ((double) num_part );
    const double jfraction = ((double) ipart + 1 ) / ((double) num_part );
    const unsigned row_beg = (unsigned)( number_row * ifraction );
    const unsigned row_end = (unsigned)( number_row * jfraction );

    txblas_SparseMatrixEntry * i = a_begin + p_begin[row_beg];
    txblas_SparseMatrixEntry * j = a_begin + p_begin[row_end];

    if ( 1 < row_end - row_beg ) {
      const double kfraction = ( ifraction + jfraction ) * 0.5 ;
      const unsigned row_split = (unsigned)( number_row * kfraction );
      p_begin[ row_split ] = txblas_row_split( i , j , row_split ) - a_begin ;
    }
    else {
      const unsigned number = j - i ;

      const size_t size = ( (unsigned char *)(i+1)) -
                          ( (unsigned char *)(i  ));

      qsort( i, number , size, & txblas_cr_compare );
    }
  }

  return 0 ;
}

int txblas_cr_prep(
  const unsigned number_row ,
  const unsigned number_coef ,
        txblas_SparseMatrixEntry a[] ,
        unsigned blocking[] )
{
  txblasTask_MatrixPrep data ;

  data.number_row = number_row ;
  data.block_size = 0 ;
  data.p_begin    = blocking ;
  data.a_begin    = a ;

  const unsigned nlock = 0 ;

  blocking[0] = 0 ;
  blocking[ number_row ] = number_coef ;

  data.iter_level = 0 ;
  phdmesh_taskpool_run( & txblas_task_cr_prep , & data , nlock );

  while ( ( 1u << data.iter_level ) < number_row ) {
    ++( data.iter_level );
    phdmesh_taskpool_run( & txblas_task_cr_prep , & data , nlock );
  }

#if 0
  {
    unsigned i = 0 ;
    unsigned j = 0 ;
    while ( j < number_row ) {
      if ( blocking[j] != i ) {
        abort();
      }
      while ( i < number_coef && a[i].row == j ) {
        ++i ;
      }
      ++j ;
    }
    if ( blocking[j] != number_coef ) {
      abort();
    }
  }
  {
    for ( unsigned i = 1 ; i < number_coef ; ++i ) {
      if ( a[i].row < a[i-1].row ) {
        abort();
      }
      else if ( a[i].row == a[i-1].row ) {
        if ( a[i].col <= a[i-1].col ) {
          abort();
        }
      }
    }
  }
#endif

  return 0 ;
}

