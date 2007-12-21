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

#include <stddef.h>
#include <math.h>
#include <txblas/reduction.h>

struct TaskXY {
  const double   alpha ;
  const double * x_beg ;
        double * y_beg ;
        unsigned number ;
};

static void task_axpby_work( void * , TPI_ThreadPool , int );

void txdaxpy( TPI_ThreadPool pool ,
              unsigned n , double a , const double * x , double * y )
{
  struct TaskXY data = { a , x , y , n };
  TPI_Run( pool , & task_axpby_work , & data , 0 );
}

/*--------------------------------------------------------------------*/

static void task_axpby_work( void * arg , TPI_ThreadPool pool , int p_rank )
{
  int p_size ;

  if ( ! TPI_Pool_size( pool , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const unsigned p_next = p_rank + 1 ;
    const unsigned n = t->number ;
    const unsigned i = ( n * p_rank ) / p_size ;
    const unsigned j = ( n * p_next ) / p_size ;

    const double a = t->alpha ;

    const double * const xe = t->x_beg + j ;
    const double *       x  = t->x_beg + i ;
          double *       y  = t->y_beg + i ;

    while ( xe != x ) { *y += a * *x ; ++x ; ++y ; }
  }
}


