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
  const double   beta ;
  const double * x_beg ;
        double * y_beg ;
        unsigned number ;
};

static void task_axpby_work( void * , TPI_ThreadPool );

void tdaxpby( TPI_ThreadPool pool ,
              unsigned n , double a , const double * x , double b , double * y )
{
  struct TaskXY data = { a , b , x , y , n };
  TPI_Run( pool , & task_axpby_work , & data );
}

/*--------------------------------------------------------------------*/

static void task_axpby_work( void * arg , TPI_ThreadPool pool )
{
  enum { STRIDE = 8 };
  int p_size ;
  int p_rank ;

  if ( ! TPI_Pool_rank( pool , & p_rank , & p_size ) ) {

    struct TaskXY * const t = (struct TaskXY *) arg ;

    const unsigned n_beg = ( t->number * ( p_rank     ) ) / p_size ;
    const unsigned n_end = ( t->number * ( p_rank + 1 ) ) / p_size ;

    const double a = t->alpha ;
    const double b = t->beta ;

    const double * const xe = t->x_beg + n_end ;
    const double * const xb = xe - ( n_end - n_beg ) % STRIDE ;
    const double *       x  = t->x_beg + n_beg ;
          double *       y  = t->y_beg + n_beg ;

    for ( ; x < xb ; x += STRIDE , y += STRIDE ) {
      y[0] = a * x[0] + b * y[0] ;
      y[1] = a * x[1] + b * y[1] ;
      y[2] = a * x[2] + b * y[2] ;
      y[3] = a * x[3] + b * y[3] ;
      y[4] = a * x[4] + b * y[4] ;
      y[5] = a * x[5] + b * y[5] ;
      y[6] = a * x[6] + b * y[6] ;
      y[7] = a * x[7] + b * y[7] ;
    }

    for ( ; x < xe ; ++x , ++y ) {
      *y = a * *x + b * *y ;
    }
  }
}


