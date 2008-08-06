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

#ifndef util_ArrayVector_hpp
#define util_ArrayVector_hpp

#include <util/Array.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class ArrayNatural< std::vector<Scalar> ,
                    Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8 >
{
public:
  // Required by all arrays:
  typedef Scalar           value_type ;
  typedef int              index_type ;
  typedef size_t           size_type ;
  typedef const ArrayDimTag * tag_type ;

  // Required by compile-time knowledge arrays:

  typedef ArrayDimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> TagList ;

  enum { Rank       = TagList::Rank };
  enum { Natural    = true };
  enum { Reverse    = false };
  enum { Contiguous = true };

  template < unsigned Ord >
  struct Tag {
    typedef typename ArrayDimTagListAt<TagList,Ord>::type type ;
  };

private:

  typedef ArrayHelper< TagList::Rank > Helper ;

  std::vector<Scalar> m_vec ;
  size_t m_stride[ TagList::Rank ];

  inline Scalar & value( index_type i ) const
    { return const_cast< Scalar & >( m_vec[i] ); }

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayNatural ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return Rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( Rank , ord );
      return array_tags<TagList>()[ord] ;
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  template < unsigned Ord > size_type dimension() const
    {
      array_bounds_check_ordinal_is_less<Rank,Ord>();
      enum { I = ( Rank - 1 ) - Ord };
      return I ? m_stride[I] / m_stride[I-1] : m_stride[I] ;
    }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( Rank , ord );
      const int i = ( Rank - 1 ) - ord ;
      return i ? m_stride[i] / m_stride[i-1] : m_stride[i] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( Rank );
      n[Rank-1] = m_stride[0] ;
      for ( int i = 1 ; i < Rank ; ++i ) {
        n[ ( Rank - 1 ) - i ] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return & value(0); }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return value(i);
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return value( ARRAY_INDEX_NATURAL_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) );
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return value( ARRAY_INDEX_NATURAL_7(m_stride,i1,i2,i3,i4,i5,i6,i7) );
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return value( ARRAY_INDEX_NATURAL_6(m_stride,i1,i2,i3,i4,i5,i6) );
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return value( ARRAY_INDEX_NATURAL_5(m_stride,i1,i2,i3,i4,i5) );
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return value( ARRAY_INDEX_NATURAL_4(m_stride,i1,i2,i3,i4) );
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return value( ARRAY_INDEX_NATURAL_3(m_stride,i1,i2,i3) );
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return value( ARRAY_INDEX_NATURAL_2(m_stride,i1,i2) );
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      ARRAY_BOUNDS_CHECKING_1(i1);
      return value( i1 );
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayNatural() : m_vec() { Helper::zero( m_stride ); }

  /** \brief Deep copy */
  ArrayNatural( const ArrayNatural & rhs )
    : m_vec( rhs.m_vec ) { Helper::copy( m_stride , rhs.m_stride ); }

  /** \brief Deep copy */
  ArrayNatural & operator = ( const ArrayNatural & rhs )
    {
      m_vec = rhs.m_vec ;
      Helper::copy( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------

  typedef typename Helper::template
    Truncate<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::natural_type
      TruncatedViewType ;

  TruncatedViewType truncated_view( index_type i ) const
    {
      return TruncatedViewType( & value( m_stride[ Rank - 2 ] * i ),
                                ArrayStride(), m_stride );
    }

  //----------------------------------
  // Class specific constructors:

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 , size_type n6 , size_type n7 , size_type n8 )
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_STRIDE_NATURAL_8( m_stride, n1, n2, n3, n4, n5, n6, n7, n8 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 , size_type n6 , size_type n7 )
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_STRIDE_NATURAL_7( m_stride, n1, n2, n3, n4, n5, n6, n7 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 , size_type n6 )
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_STRIDE_NATURAL_6( m_stride, n1, n2, n3, n4, n5, n6 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 )
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_STRIDE_NATURAL_5( m_stride, n1, n2, n3, n4, n5 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 )
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_STRIDE_NATURAL_4( m_stride, n1, n2, n3, n4 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 )
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_STRIDE_NATURAL_3( m_stride, n1, n2, n3 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 )
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_STRIDE_NATURAL_2( m_stride, n1, n2 );
      m_vec.resize( size() );
    }

  void resize(size_type n1 )
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      m_stride[0] = n1 ;
      m_vec.resize( size() );
    }

  ArrayNatural( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 , size_type n6 , size_type n7 , size_type n8 )
    : m_vec() { resize(n1,n2,n3,n4,n5,n6,n7,n8); }

  ArrayNatural( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 , size_type n6 , size_type n7 )
    : m_vec() { resize(n1,n2,n3,n4,n5,n6,n7); }

  ArrayNatural( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 , size_type n6 )
    : m_vec() { resize(n1,n2,n3,n4,n5,n6); }

  ArrayNatural( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 )
    : m_vec() { resize(n1,n2,n3,n4,n5); }

  ArrayNatural( size_type n1 , size_type n2 , size_type n3 , size_type n4 )
    : m_vec() { resize(n1,n2,n3,n4); }

  ArrayNatural( size_type n1 , size_type n2 , size_type n3 )
    : m_vec() { resize(n1,n2,n3); }

  ArrayNatural( size_type n1 , size_type n2 )
    : m_vec() { resize(n1,n2); }

  explicit ArrayNatural( size_type n1 )
    : m_vec() { resize(n1); }

  //----------------------------------
  // Other:

  template< class TA >
  struct Append {
    typedef typename Helper::template
      Append<std::vector<Scalar>,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8,TA>::natural_type
        type ;
  };
};

//----------------------------------------------------------------------

template< typename Scalar ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class ArrayFortran< std::vector<Scalar> ,
                    Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7, Tag8 >
{
public:
  // Required by all arrays:
  typedef Scalar           value_type ;
  typedef int              index_type ;
  typedef size_t           size_type ;
  typedef const ArrayDimTag * tag_type ;

  // Required by compile-time knowledge arrays:

  typedef ArrayDimTagList<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> TagList ;

  enum { Rank       = TagList::Rank };
  enum { Natural    = false };
  enum { Reverse    = true };
  enum { Contiguous = true };

  template < unsigned Ord >
  struct Tag {
    typedef typename ArrayDimTagListAt<TagList,Ord>::type type ;
  };

private:

  typedef ArrayHelper< TagList::Rank > Helper ;

  std::vector<Scalar> m_vec ;
  size_t m_stride[ TagList::Rank ];

  inline
  Scalar & value( index_type i ) const
    { return const_cast< Scalar & >( m_vec[i] ); }

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayNatural ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return Rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( Rank , ord );
      return array_tags<TagList>()[ord] ;
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  template < unsigned Ord > size_type dimension() const
    {
      array_bounds_check_ordinal_is_less<Rank,Ord>();
      return Ord ? m_stride[Ord] / m_stride[Ord-1] : m_stride[Ord] ;
    }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( Rank , ord );
      return ord ? m_stride[ord] / m_stride[ord-1] : m_stride[ord] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( Rank );
      n[0] = m_stride[0] ;
      for ( int i = 1 ; i < Rank ; ++i ) {
        n[i] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return & value(0); }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return value( i );
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return value( ARRAY_INDEX_FORTRAN_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) );
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return value( ARRAY_INDEX_FORTRAN_7(m_stride,i1,i2,i3,i4,i5,i6,i7) );
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return value( ARRAY_INDEX_FORTRAN_6(m_stride,i1,i2,i3,i4,i5,i6) );
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return value( ARRAY_INDEX_FORTRAN_5(m_stride,i1,i2,i3,i4,i5) );
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return value( ARRAY_INDEX_FORTRAN_4(m_stride,i1,i2,i3,i4) );
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return value( ARRAY_INDEX_FORTRAN_3(m_stride,i1,i2,i3) );
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return value( ARRAY_INDEX_FORTRAN_2(m_stride,i1,i2) );
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      ARRAY_BOUNDS_CHECKING_1(i1);
      return value( i1 );
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayFortran() : m_vec() { Helper::zero( m_stride ); }

  ArrayFortran( const ArrayFortran & rhs )
    : m_vec( rhs.m_vec ) { Helper::copy( m_stride , rhs.m_stride ); }

  ArrayFortran & operator = ( const ArrayFortran & rhs )
    {
      m_vec = rhs.m_vec ;
      Helper::copy( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename Helper::template
    Truncate<Scalar,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::fortran_type
      TruncatedViewType ;

  TruncatedViewType truncated_view( index_type i ) const
    {
      return TruncatedViewType( & value( m_stride[ Rank - 2 ] * i ),
                                ArrayStride() , m_stride );
    }

  //----------------------------------
  // Class specific constructors:

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 , size_type n6 , size_type n7 , size_type n8 )
    {
      array_bounds_check_rank_is_equal<Rank,8>();
      ARRAY_STRIDE_FORTRAN_8( m_stride, n1, n2, n3, n4, n5, n6, n7, n8 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 , size_type n6 , size_type n7 )
    {
      array_bounds_check_rank_is_equal<Rank,7>();
      ARRAY_STRIDE_FORTRAN_7( m_stride, n1, n2, n3, n4, n5, n6, n7 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 , size_type n6 )
    {
      array_bounds_check_rank_is_equal<Rank,6>();
      ARRAY_STRIDE_FORTRAN_6( m_stride, n1, n2, n3, n4, n5, n6 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
               size_type n5 )
    {
      array_bounds_check_rank_is_equal<Rank,5>();
      ARRAY_STRIDE_FORTRAN_5( m_stride, n1, n2, n3, n4, n5 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 , size_type n4 )
    {
      array_bounds_check_rank_is_equal<Rank,4>();
      ARRAY_STRIDE_FORTRAN_4( m_stride, n1, n2, n3, n4 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 , size_type n3 )
    {
      array_bounds_check_rank_is_equal<Rank,3>();
      ARRAY_STRIDE_FORTRAN_3( m_stride, n1, n2, n3 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 , size_type n2 )
    {
      array_bounds_check_rank_is_equal<Rank,2>();
      ARRAY_STRIDE_FORTRAN_2( m_stride, n1, n2 );
      m_vec.resize( size() );
    }

  void resize( size_type n1 )
    {
      array_bounds_check_rank_is_equal<Rank,1>();
      m_stride[0] = n1 ;
      m_vec.resize( size() );
    }

  ArrayFortran( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 , size_type n6 , size_type n7 , size_type n8 )
    : m_vec() { resize(n1,n2,n3,n4,n5,n6,n7,n8); }

  ArrayFortran( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 , size_type n6 , size_type n7 )
    : m_vec() { resize(n1,n2,n3,n4,n5,n6,n7); }

  ArrayFortran( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 , size_type n6 )
    : m_vec() { resize(n1,n2,n3,n4,n5,n6); }

  ArrayFortran( size_type n1 , size_type n2 , size_type n3 , size_type n4 ,
                size_type n5 )
    : m_vec() { resize(n1,n2,n3,n4,n5); }

  ArrayFortran( size_type n1 , size_type n2 , size_type n3 , size_type n4 )
    : m_vec() { resize(n1,n2,n3,n4); }

  ArrayFortran( size_type n1 , size_type n2 , size_type n3 )
    : m_vec() { resize(n1,n2,n3); }

  ArrayFortran( size_type n1 , size_type n2 )
    : m_vec() { resize(n1,n2); }

  explicit ArrayFortran( size_type n1 )
    : m_vec() { resize(n1); }

  //----------------------------------
  // Other:

  template< class TA >
  struct Append {
    typedef typename Helper::template
      Append<std::vector<Scalar>,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8,TA>::fortran_type
        type ;
  };
};

//----------------------------------------------------------------------

template< typename Scalar >
class ArrayNatural< std::vector<Scalar> ,
                    void,void,void,void,void,void,void,void>
{
public:

  typedef Scalar           value_type ;
  typedef int              index_type ;
  typedef size_t           size_type ;
  typedef const ArrayDimTag * tag_type ;

  enum { Natural    = true };
  enum { Reverse    = false };
  enum { Contiguous = true };

private:

  std::vector<Scalar> m_vec ;
  size_type m_rank ;
  size_type m_stride[8];
  tag_type  m_tag[8] ;

  inline Scalar & value( index_type i ) const
    { return const_cast< Scalar & >( m_vec[i] ); }

  typedef ArrayHelper<8> Helper ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayNatural ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return m_rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( m_rank , ord );
      return m_tag[ ( m_rank - 1 ) - ord ];
    }

  //----------------------------------

  size_type size() const { return m_stride[ m_rank - 1 ]; }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( m_rank , ord );
      const int i = ( m_rank - 1 ) - ord ;
      return i ? m_stride[i] / m_stride[i-1] : m_stride[i] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( m_rank );
      n[m_rank-1] = m_stride[0] ;
      for ( int i = 1 ; i < m_rank ; ++i ) {
        n[ ( m_rank - 1 ) - i ] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return & value(0); }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return value( i );
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return value( ARRAY_INDEX_NATURAL_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) );
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return value( ARRAY_INDEX_NATURAL_7(m_stride,i1,i2,i3,i4,i5,i6,i7) );
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return value( ARRAY_INDEX_NATURAL_6(m_stride,i1,i2,i3,i4,i5,i6) );
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return value( ARRAY_INDEX_NATURAL_5(m_stride,i1,i2,i3,i4,i5) );
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return value( ARRAY_INDEX_NATURAL_4(m_stride,i1,i2,i3,i4) );
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return value( ARRAY_INDEX_NATURAL_3(m_stride,i1,i2,i3) );
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return value( ARRAY_INDEX_NATURAL_2(m_stride,i1,i2) );
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      ARRAY_BOUNDS_CHECKING_1(i1);
      return value( i1 );
    }

  //----------------------------------
  /** \brief Truncate view of the array */

  typedef ArrayNatural<Scalar> TruncateViewType ;

  TruncateViewType truncated_view( index_type i ) const
    {
      return TruncateViewType( & value( m_stride[ m_rank - 2 ] * i ),
                               m_rank - 1 ,
                               ArrayStride() ,
                               m_stride ,
                               m_tag );
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayNatural()
    : m_vec(), m_rank(0)
    {
      Helper::zero( m_stride );
      Helper::zero( m_tag );
    }

  ArrayNatural( const ArrayNatural & rhs )
    : m_vec( rhs.m_vec ), m_rank( rhs.m_rank )
    {
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
    }

  ArrayNatural & operator = ( const ArrayNatural & rhs )
    {
      m_vec  = rhs.m_vec ;
      m_rank = rhs.m_rank ;
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
      return *this ;
    }

  //----------------------------------
  // Class specific constructors:

  ArrayNatural( const std::vector<size_type> & arg_size ,
                const std::vector<tag_type>  & arg_tag )
    : m_vec(),
      m_rank( arg_size.size() == arg_tag.size() ? arg_size.size() : 0 )
    {
      Helper::stride_natural( m_rank , m_stride , & arg_size[0] );
      Helper::tag_natural(    m_rank , m_tag ,    & arg_tag[0] );
      m_vec.resize( size() );
    }
};

//----------------------------------------------------------------------

template< typename Scalar >
class ArrayFortran< std::vector<Scalar>,
                    void,void,void,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef int                 index_type ;
  typedef size_t              size_type ;
  typedef const ArrayDimTag * tag_type ;

  enum { Natural    = false };
  enum { Reverse    = true };
  enum { Contiguous = true };

private:

  std::vector<Scalar> m_vec ;
  size_type m_rank ;
  size_type m_stride[8];
  tag_type  m_tag[8] ;

  inline Scalar & value( index_type i ) const
    { return const_cast< Scalar & >( m_vec[i] ); }

  typedef ArrayHelper<8> Helper ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayNatural ;

  template< typename ,
            class , class , class , class ,
            class , class , class , class >
  friend class ArrayFortran ;

public:

  //----------------------------------

  unsigned rank() const { return m_rank ; }
  bool natural() const { return Natural ; }
  bool reverse() const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  tag_type tag( unsigned ord ) const
    {
      array_bounds_checking( m_rank , ord );
      return m_tag[ ord ];
    }

  //----------------------------------

  size_type size() const { return m_stride[ m_rank - 1 ]; }

  size_type dimension( index_type ord ) const
    {
      array_bounds_checking( m_rank , ord );
      return ord ? m_stride[ord] / m_stride[ord-1] : m_stride[ord] ;
    }

  void dimensions( std::vector<size_type> & n )
    {
      n.resize( m_rank );
      n[0] = m_stride[0] ;
      for ( int i = 1 ; i < m_rank ; ++i ) {
        n[i] = m_stride[i] / m_stride[i-1] ;
      }
    }
  
  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return & value(0); }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_BOUNDS_CHECKING_ORDINAL( size() , i );
      return value( i );
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 , index_type i8 ) const
    {
      ARRAY_BOUNDS_CHECKING_8(i1,i2,i3,i4,i5,i6,i7,i8);
      return value( ARRAY_INDEX_FORTRAN_8(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) );
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ,
                           index_type i7 ) const
    {
      ARRAY_BOUNDS_CHECKING_7(i1,i2,i3,i4,i5,i6,i7);
      return value( ARRAY_INDEX_FORTRAN_7(m_stride,i1,i2,i3,i4,i5,i6,i7) );
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 , index_type i6 ) const
    {
      ARRAY_BOUNDS_CHECKING_6(i1,i2,i3,i4,i5,i6);
      return value( ARRAY_INDEX_FORTRAN_6(m_stride,i1,i2,i3,i4,i5,i6) );
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ,
                           index_type i5 ) const
    {
      ARRAY_BOUNDS_CHECKING_5(i1,i2,i3,i4,i5);
      return value( ARRAY_INDEX_FORTRAN_5(m_stride,i1,i2,i3,i4,i5) );
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 , index_type i4 ) const
    {
      ARRAY_BOUNDS_CHECKING_4(i1,i2,i3,i4);
      return value( ARRAY_INDEX_FORTRAN_4(m_stride,i1,i2,i3,i4) );
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ,
                           index_type i3 ) const
    {
      ARRAY_BOUNDS_CHECKING_3(i1,i2,i3);
      return value( ARRAY_INDEX_FORTRAN_3(m_stride,i1,i2,i3) );
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( index_type i1 , index_type i2 ) const
    {
      ARRAY_BOUNDS_CHECKING_2(i1,i2);
      return value( ARRAY_INDEX_FORTRAN_2(m_stride,i1,i2) );
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( index_type i1 ) const
    {
      ARRAY_BOUNDS_CHECKING_1(i1);
      return value( i1 );
    }

  //----------------------------------
  /** \brief Truncate view of the array */

  typedef ArrayFortran<Scalar> TruncateViewType ;

  TruncateViewType truncated_view( index_type i ) const
    {
      return TruncateViewType( & value( m_stride[ m_rank - 2 ] * i ),
                               m_rank - 1 ,
                               ArrayStride() ,
                               m_stride ,
                               m_tag );
    }

  //----------------------------------
  // Required constructors and assignment operators:

  ArrayFortran()
    : m_vec(), m_rank(0)
    {
      Helper::zero( m_stride );
      Helper::zero( m_tag );
    }

  ArrayFortran( const ArrayFortran & rhs )
    : m_vec( rhs.m_vec ), m_rank( rhs.m_rank )
    {
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
    }

  ArrayFortran & operator = ( const ArrayFortran & rhs )
    {
      m_vec  = rhs.m_vec ;
      m_rank = rhs.m_rank ;
      Helper::copy( m_stride , rhs.m_stride );
      Helper::copy( m_tag ,    rhs.m_tag );
      return *this ;
    }

  //----------------------------------
  // Class specific constructors:

  ArrayFortran( const std::vector<size_type> & arg_size ,
                const std::vector<tag_type>  & arg_tag )
    : m_vec(),
      m_rank( arg_size.size() == arg_tag.size() ? arg_size.size() : 0 )
    {
      Helper::stride_fortran( m_rank , m_stride , & arg_size[0] );
      Helper::tag_fortran(    m_rank , m_tag ,    & arg_tag[0] );
      m_vec.resize();
    }
};

//----------------------------------------------------------------------

/** \endcond */

}

#endif

