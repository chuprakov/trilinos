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
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   November 2006
 */

#ifndef util_ParallelReduce_hpp
#define util_ParallelReduce_hpp

#include <cstddef>
#include <iosfwd>
#include <util/FixedArray.hpp>
#include <util/Parallel.hpp>

//------------------------------------------------------------------------

namespace phdmesh {

/** Write string from any or all processors
 *  to the ostream on the root processor.
 */
void all_write( ParallelMachine     arg_comm ,
                std::ostream      & arg_root_os ,
                const std::string & arg_msg );

void all_reduce_sum( ParallelMachine ,
                     const double * local , double * global , unsigned count );

void all_reduce_sum( ParallelMachine ,
                     const float * local , float * global , unsigned count );

void all_reduce_sum( ParallelMachine ,
                     const int * local , int * global , unsigned count );

void all_reduce_bor( ParallelMachine ,
                     const unsigned * local ,
                     unsigned * global , unsigned count );

/** Aggregated parallel in-place reduce-to-all-processors operations.
 *
 *  example:
 *    ParallelMachine comm = ... ;
 *    double a[5] ;
 *    int    b[3] ;
 *    all_reduce( comm , ReduceSum<5>( a ) & ReduceMax<3>( b ) );
 *
 *  Reduction options include:
 *    ReduceSum   <N>( T * )   // Summation
 *    ReduceProd  <N>( T * )   // Product
 *    ReduceMax   <N>( T * )   // Maximum
 *    ReduceMin   <N>( T * )   // Minimum
 *    ReduceBitOr <N>( T * )   // Bit-wise OR
 *    ReduceBitAnd<N>( T * )   // Bit-wise AND
 */
template < class ReduceOp >
void all_reduce( ParallelMachine , const ReduceOp & );

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {
namespace {
// Blank namespace so that this class produces local symbols,
// avoiding complaints from a linker of multiple-define symbols.

struct ReduceEnd {
  struct WorkType {};
  void copyin(  WorkType & ) const {}
  void copyout( WorkType & ) const {}
  static void op( WorkType & , WorkType & ) {}
};

// Workhorse class for aggregating reduction operations.

template <class Op, typename T, class Next>
struct Reduce {

  typedef T Type ;
  enum { N = Op::N };

  struct WorkType {
    typename Next::WorkType m_next ;
    Type                    m_value[N];
  };

  Next   m_next ;
  Type * m_value ;

  // Copy values into buffer:
  void copyin( WorkType & w ) const
    { Copy<N>( w.m_value , m_value ); m_next.copyin( w.m_next ); }
      
  // Copy value out from buffer:
  void copyout( WorkType & w ) const
    { Copy<N>( m_value , w.m_value ); m_next.copyout( w.m_next ); }

  // Reduction function
  static void op( WorkType & out , WorkType & in )
    { Op( out.m_value , in.m_value ); Next::op( out.m_next , in.m_next ); }

  // Aggregate reduction operations, use '&' for left-to-right evaluation
  template<class OpB, typename TB>
  Reduce<OpB, TB, Reduce<Op,T,Next> >
    operator & ( const Reduce<OpB,TB,ReduceEnd> & rhs )
      { return Reduce<OpB, TB, Reduce<Op,T,Next> >( rhs , *this ); }

  // Constructor for aggregation:
  Reduce( const Reduce<Op,T, ReduceEnd> & arg_val , const Next & arg_next )
    : m_next( arg_next ), m_value( arg_val.m_value ) {}

  // Constructor for aggregate member:
  explicit Reduce( Type * arg_value )
   : m_next(), m_value( arg_value ) {}

  static void void_op( void*inv, void*inoutv, int*, ParallelDatatype*);
};

template <class Op, typename T, class Next>
void Reduce<Op,T,Next>::void_op( void*inv, void*inoutv,int*,ParallelDatatype*)
{
  op( * reinterpret_cast<WorkType*>( inoutv ) ,
      * reinterpret_cast<WorkType*>( inv ) );
}

}
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace phdmesh {

template<unsigned N, typename T>
inline
Reduce< Sum<N> , T, ReduceEnd> ReduceSum( T * value )
{ return Reduce< Sum<N>, T, ReduceEnd >( value ); }

template<unsigned N, typename T>
inline
Reduce< Prod<N>, T, ReduceEnd > ReduceProd( T * value )
{ return Reduce< Prod<N>, T, ReduceEnd >( value ); }

template<unsigned N, typename T>
inline
Reduce< Max<N>, T, ReduceEnd> ReduceMax( T * value )
{ return Reduce< Max<N>, T, ReduceEnd>( value ); }

template<unsigned N, typename T>
inline
Reduce< Min<N>, T, ReduceEnd> ReduceMin( T * value )
{ return Reduce<Min<N>, T, ReduceEnd>( value ); }

template<unsigned N, typename T>
inline
Reduce< BitOr<N>, T, ReduceEnd> ReduceBitOr( T * value )
{ return Reduce< BitOr<N>, T, ReduceEnd>( value ); }

template<unsigned N, typename T>
inline
Reduce< BitAnd<N>, T, ReduceEnd> ReduceBitAnd( T * value )
{ return Reduce< BitAnd<N>, T, ReduceEnd>( value ); }

//----------------------------------------------------------------------
// all_reduce( comm , ReduceSum<5>( A ) & ReduceMax<3>( B ) );

extern "C" {
typedef void (*ParallelReduceOp)
  ( void * inv , void * outv , int * , ParallelDatatype * );
}

void all_reduce( ParallelMachine  arg_comm ,
                 ParallelReduceOp arg_op ,
                 void           * arg_in ,
                 void           * arg_out ,
                 unsigned         arg_len );

namespace {

template < class ReduceOp >
void all_reduce_driver( ParallelMachine comm , const ReduceOp & op )
{
  typedef typename ReduceOp::WorkType WorkType ;

  WorkType inbuf , outbuf ;

  ParallelReduceOp f =
    reinterpret_cast<ParallelReduceOp>( & ReduceOp::void_op );
  op.copyin( inbuf );
  all_reduce( comm , f , & inbuf, & outbuf, sizeof(WorkType) );
  op.copyout( outbuf );
}

}

template < class ReduceOp >
inline
void all_reduce( ParallelMachine comm , const ReduceOp & op )
{ all_reduce_driver<ReduceOp>( comm , op ); }

}

//----------------------------------------------------------------------

#endif

