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
#include <util/taskpool.h>
#include <txblas/reduction.h>

/*--------------------------------------------------------------------*/

#define SUM_ADD( v , a ) \
  { \
    const double S_SUM_ADD = v[0] + a ; \
    if ( a < v[0] ) { v[1] += a    - ( S_SUM_ADD - v[0] ); } \
    else            { v[1] += v[0] - ( S_SUM_ADD - a    ); } \
    v[0]  = S_SUM_ADD + v[1] ; \
    v[1] += S_SUM_ADD - v[0] ; \
  }

#define SUM_SUB( v , a ) \
  { \
    const double S_SUM_SUB = v[0] - a ; \
    if ( a < fabs(v[0]) ) { v[1] +=  -a  - ( S_SUM_SUB - v[0] ); } \
    else                  { v[1] += v[0] - ( S_SUM_SUB + a    ); } \
    v[0]  = S_SUM_SUB + v[1] ; \
    v[1] += S_SUM_SUB - v[0] ; \
  }

/*--------------------------------------------------------------------*/

void xdsum_add_value( double * v , double a )
{
  if ( a < 0 ) { a = -a ; v += 2 ; }
  SUM_ADD( v , a );
}

void xdsum_add_dsum( double * v , const double * const a )
{
  SUM_ADD( v , a[0] );
  SUM_ADD( v , a[1] );
  v += 2 ;
  SUM_ADD( v , a[2] );
  SUM_ADD( v , a[3] );
}

void xdsum_get_value( double * const y , const double * const v )
{
  y[0] = v[0] ;
  y[1] = v[1] ;
  SUM_SUB( y , v[2] );
  SUM_SUB( y , v[3] );
}

/*--------------------------------------------------------------------*/

struct TaskX {
        double * x_sum ;
  const double * x_beg ;
        unsigned number ;
};

struct TaskXY {
        double * xy_sum ;
  const double * x_beg ;
  const double * y_beg ;
        unsigned number ;
};

/*--------------------------------------------------------------------*/

static
void add_array( double * const s , const double * x , const double * const xe )
{
  while ( xe != x ) {
    double   a = *x ; ++x ;
    double * p = s ;
    if ( a < 0 ) { a = -a ; p += 2 ; }
    SUM_ADD( p , a );
  }
}

void xdsum_add_array( double * s , unsigned n , const double * x )
{ add_array( s , x , x + n ); }

static int task_sum_work( void * arg , unsigned p_size , unsigned p_rank )
{
  struct TaskX * const t  = (struct TaskX *) arg ;

  const unsigned p_next = p_rank + 1 ;
  const unsigned n = t->number ;
  const double * const xb = t->x_beg + ( n * p_rank ) / p_size ;
  const double * const xe = t->x_beg + ( n * p_next ) / p_size ;
        double * const v  = t->x_sum ;

  double partial[4] = { 0 , 0 , 0 , 0 };

  add_array( partial , xb , xe );

  phdmesh_taskpool_lock(0,NULL);
  xdsum_add_dsum( v , partial );
  phdmesh_taskpool_unlock(0,NULL);

  return 0 ;
}

void txdsum_add_array( double * s , unsigned n , const double * x )
{
  struct TaskX data ;
  data.x_sum  = s ;
  data.x_beg  = x ;
  data.number = n ;
  phdmesh_taskpool_run( & task_sum_work , & data , 1 );
}

/*--------------------------------------------------------------------*/

static void norm1( double * s , const double * x , const double * const xe )
{
  while ( xe != x ) { const double a = fabs(*x); ++x ; SUM_ADD( s , a ); }
}

void xdnorm1( double * s2 , unsigned n , const double * x )
{ norm1( s2 , x , x + n ); }

static int task_norm1_work( void * arg , unsigned p_size , unsigned p_rank )
{
  struct TaskX * const t  = (struct TaskX *) arg ;

  const unsigned p_next = p_rank + 1 ;
  const unsigned n = t->number ;
  const double * const xb = t->x_beg + ( n * p_rank ) / p_size ;
  const double * const xe = t->x_beg + ( n * p_next ) / p_size ;
        double * const v  = t->x_sum ;

  double partial[2] = { 0 , 0 };

  norm1( partial , xb , xe );

  phdmesh_taskpool_lock(0,NULL);
  SUM_ADD( v , partial[0] );
  SUM_ADD( v , partial[1] );
  phdmesh_taskpool_unlock(0,NULL);

  return 0 ;
}

void txdnorm1( double * s , unsigned n , const double * x )
{
  struct TaskX data ;
  data.x_sum  = s ;
  data.x_beg  = x ;
  data.number = n ;
  phdmesh_taskpool_run( & task_norm1_work , & data , 1 );
}

/*--------------------------------------------------------------------*/

static void dot1( double * s , const double * x , const double * const xe )
{
  while ( xe != x ) { const double a = *x * *x ; ++x ; SUM_ADD( s , a ); }
}

void xddot1( double * s2 , unsigned n , const double * x )
{ dot1( s2 , x , x + n ); }

static int task_dot_x_work( void * arg , unsigned p_size , unsigned p_rank )
{
  struct TaskX * const t  = (struct TaskX *) arg ;

  const unsigned p_next = p_rank + 1 ;
  const unsigned n = t->number ;
  const double * const xb = t->x_beg + ( n * p_rank ) / p_size ;
  const double * const xe = t->x_beg + ( n * p_next ) / p_size ;
        double * const v  = t->x_sum ;

  double partial[2] = { 0 , 0 };

  dot1( partial , xb , xe );

  phdmesh_taskpool_lock(0,NULL);
  SUM_ADD( v , partial[0] );
  SUM_ADD( v , partial[1] );
  phdmesh_taskpool_unlock(0,NULL);

  return 0 ;
}

void txddot1( double * s , unsigned n , const double * x )
{
  struct TaskX data ;
  data.x_sum  = s ;
  data.x_beg  = x ;
  data.number = n ;
  phdmesh_taskpool_run( & task_dot_x_work , & data , 1 );
}

/*--------------------------------------------------------------------*/

static void dot( double * const s ,
                 const double * x ,
                 const double * y , const double * const ye )
{
  while ( ye != y ) {
    double   a = *x * *y ; ++x ; ++y ;
    double * p = s ;
    if ( a < 0 ) { a = -a ; p += 2 ; }
    SUM_ADD( p , a );
  }
}

void xddot( double * s4 , unsigned n , const double * x , const double * y )
{ dot( s4 , x , y , y + n ); }

static int task_dot_xy_work( void * arg , unsigned p_size , unsigned p_rank )
{
  struct TaskXY * const t  = (struct TaskXY *) arg ;

  const unsigned p_next = p_rank + 1 ;
  const unsigned n = t->number ;
  const unsigned i = ( n * p_rank ) / p_size ;
  const unsigned j = ( n * p_next ) / p_size ;

  const double * const xb = t->x_beg + i ;
  const double * const yb = t->y_beg + i ;
  const double * const ye = t->y_beg + j ;
        double * const s  = t->xy_sum ;

  double partial[4] = { 0 , 0 , 0 , 0 };

  dot( partial , xb , yb , ye );

  phdmesh_taskpool_lock(0,NULL);
  xdsum_add_dsum( s , partial );
  phdmesh_taskpool_unlock(0,NULL);

  return 0 ;
}

void txddot( double * s ,
             unsigned n , const double * x , const double * y )
{
  struct TaskXY data ;
  data.xy_sum = s ;
  data.x_beg  = x ;
  data.x_beg  = y ;
  data.number = n ;
  phdmesh_taskpool_run( & task_dot_xy_work , & data , 1 );
}

/*--------------------------------------------------------------------*/

