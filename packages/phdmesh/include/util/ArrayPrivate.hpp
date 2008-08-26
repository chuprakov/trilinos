/*------------------------------------------------------------------------*/
/*      phdMesh : Parallel Heterogneous Dynamic unstructured Mesh         */
/*                Copyright (2008) Sandia Corporation                     */
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
 * @date   June 2008
 */

#ifndef util_ArrayPrivate_hpp
#define util_ArrayPrivate_hpp

namespace phdmesh {

//----------------------------------------------------------------------

size_t array_stride_size( unsigned rank , const size_t * const stride );

void array_stride_to_natural_dimensions(
  unsigned rank , const size_t * const stride , unsigned * const dim );

void array_stride_to_natural_indices(
  unsigned rank , const size_t * const stride ,
  size_t offset , unsigned * const indices );

//----------------------------------------------------------------------

void array_check_rank(    unsigned rank , unsigned test_rank );
void array_check_ordinal( unsigned rank , unsigned test_ordinal );
void array_check_index(   size_t size , unsigned test_ordinal );

void array_check_indices( const bool ,
                          const unsigned ,
                          const size_t * const ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 ,
                          const unsigned = 0 );

//----------------------------------------------------------------------
/** \function array_check_order_known
 *  \brief Template function defined for known array ordering.
 */
template< ArrayOrder > void array_check_order_known();

template<> inline void array_check_order_known< NaturalOrder >() {}
template<> inline void array_check_order_known< FortranOrder >() {}

//----------------------------------------------------------------------

template< unsigned , unsigned > struct array_check_ordinal_is_less ;

template<> struct array_check_ordinal_is_less<0,8> {};
template<> struct array_check_ordinal_is_less<1,8> {};
template<> struct array_check_ordinal_is_less<2,8> {};
template<> struct array_check_ordinal_is_less<3,8> {};
template<> struct array_check_ordinal_is_less<4,8> {};
template<> struct array_check_ordinal_is_less<5,8> {};
template<> struct array_check_ordinal_is_less<6,8> {};
template<> struct array_check_ordinal_is_less<7,8> {};

template<> struct array_check_ordinal_is_less<0,7> {};
template<> struct array_check_ordinal_is_less<1,7> {};
template<> struct array_check_ordinal_is_less<2,7> {};
template<> struct array_check_ordinal_is_less<3,7> {};
template<> struct array_check_ordinal_is_less<4,7> {};
template<> struct array_check_ordinal_is_less<5,7> {};
template<> struct array_check_ordinal_is_less<6,7> {};

template<> struct array_check_ordinal_is_less<0,6> {};
template<> struct array_check_ordinal_is_less<1,6> {};
template<> struct array_check_ordinal_is_less<2,6> {};
template<> struct array_check_ordinal_is_less<3,6> {};
template<> struct array_check_ordinal_is_less<4,6> {};
template<> struct array_check_ordinal_is_less<5,6> {};

template<> struct array_check_ordinal_is_less<0,5> {};
template<> struct array_check_ordinal_is_less<1,5> {};
template<> struct array_check_ordinal_is_less<2,5> {};
template<> struct array_check_ordinal_is_less<3,5> {};
template<> struct array_check_ordinal_is_less<4,5> {};

template<> struct array_check_ordinal_is_less<0,4> {};
template<> struct array_check_ordinal_is_less<1,4> {};
template<> struct array_check_ordinal_is_less<2,4> {};
template<> struct array_check_ordinal_is_less<3,4> {};

template<> struct array_check_ordinal_is_less<0,3> {};
template<> struct array_check_ordinal_is_less<1,3> {};
template<> struct array_check_ordinal_is_less<2,3> {};

template<> struct array_check_ordinal_is_less<0,2> {};
template<> struct array_check_ordinal_is_less<1,2> {};

template<> struct array_check_ordinal_is_less<0,1> {};

//----------------------------------------------------------------------

template< class , unsigned > struct ArrayTagAt ;

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,0>
{ typedef Tag1 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,1>
{ typedef Tag2 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,2>
{ typedef Tag3 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,3>
{ typedef Tag4 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,4>
{ typedef Tag5 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,5>
{ typedef Tag6 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,6>
{ typedef Tag7 type ; };

template< typename Scalar , ArrayOrder order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTagAt<Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>,7>
{ typedef Tag8 type ; };

//----------------------------------------------------------------------

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,TApp> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,TApp,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,Tag5,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,TApp,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,Tag4,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,TApp,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,Tag3,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,TApp,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,Tag2,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,Tag2,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,TApp,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,Tag1,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,TApp,void,void,void,void,void,void> type ;
};

template< typename Scalar , class TApp >
struct ArrayAppend<
  Array<Scalar,RankZero,void,void,void,void,void,void,void,void> , TApp >
{
  typedef
    Array<Scalar,NaturalOrder,TApp,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class TApp >
struct ArrayAppend<
  Array<Scalar,NaturalOrder,void,void,void,void,void,void,void,void> , TApp >
{ /* Cannot append to runtime-ranked array */ };

template< typename Scalar , class TApp >
struct ArrayAppend<
  Array<Scalar,FortranOrder,void,void,void,void,void,void,void,void> , TApp >
{ /* Cannot append to runtime-ranked array */ };

//----------------------------------------------------------------------

template< class A > struct ArrayTruncate ;

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTruncate<
  Array<Scalar,NaturalOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> >
{
  typedef
    Array<Scalar,NaturalOrder,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8,void> type ;
};

template< typename Scalar , class Tag1 >
struct ArrayTruncate<
  Array<Scalar,NaturalOrder,Tag1,void,void,void,void,void,void,void> >
{
  typedef Array<Scalar,RankZero,void,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> >
{
  typedef Array<Scalar,RankZero,void,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,void,void,void,void,void,void,void> type ;
};

template< typename Scalar , class Tag1 , class Tag2 , class Tag3 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,void,void,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,void,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,void,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void> type ;
};

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct ArrayTruncate<
  Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> >
{
  typedef
    Array<Scalar,FortranOrder,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void> type ;
};

//----------------------------------------------------------------------

template< ArrayOrder , unsigned Rank , unsigned Ordinal = 0 >
struct ArrayStrideDim ;

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<RankZero,Rank,Ordinal> {
  static size_t dimension( const size_t * ) { return 0 ; }
  static size_t dimension( const size_t * , unsigned ) { return 0 ; }
};

template< unsigned Rank >
struct ArrayStrideDim<FortranOrder,Rank,0> {
  static size_t dimension( const size_t * stride ) { return stride[0]; }

  static size_t dimension( const size_t * stride , const unsigned & ordinal )
    { return ordinal ? stride[ordinal] / stride[ordinal-1] : stride[0] ; }
};

template< unsigned Rank >
struct ArrayStrideDim<NaturalOrder,Rank,0> {
  static size_t dimension( const size_t * stride ) { return stride[0]; }

  static size_t dimension( const size_t * stride , const unsigned & ordinal )
    {
      const unsigned i = ( Rank - 1 ) - ordinal ;
      return i ? stride[i] / stride[i-1] : stride[0] ;
    }
};

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<FortranOrder,Rank,Ordinal> {
  static size_t dimension( const size_t * stride )
    { return stride[Ordinal] / stride[Ordinal-1]; }
};

template< unsigned Rank , unsigned Ordinal >
struct ArrayStrideDim<NaturalOrder,Rank,Ordinal> {
  static size_t dimension( const size_t * stride )
    {
      enum { I = ( Rank - 1 ) - Ordinal };
      return stride[I] / stride[I-1];
    }
};

//----------------------------------------------------------------------

namespace {

template< class Tag > const ArrayDimTag * array_dim_tag();

template<>
inline
const ArrayDimTag * array_dim_tag<void>() { return NULL ; }

template< class Tag >
inline
const ArrayDimTag * array_dim_tag() { return & Tag::tag(); }

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
const ArrayDimTag * const * array_dim_tags()
{
  static const ArrayDimTag * t[8] =
    {
      array_dim_tag< Tag1 >() ,
      array_dim_tag< Tag2 >() ,
      array_dim_tag< Tag3 >() ,
      array_dim_tag< Tag4 >() ,
      array_dim_tag< Tag5 >() ,
      array_dim_tag< Tag6 >() ,
      array_dim_tag< Tag7 >() ,
      array_dim_tag< Tag8 >()
    };

  return t ;
}

}

//----------------------------------------------------------------------
/** \cond */
//----------------------------------------------------------------------

}

#endif

