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

#ifndef util_Array_hpp
#define util_Array_hpp

#ifdef ARRAY_BOUNDS_CHECKING
#define ARRAY_CHECK( X ) X
#else
#define ARRAY_CHECK( X )
#endif

#include <vector>
#include <string>
#include <util/SimpleArrayOps.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

/** \enum ArrayOrder
 *  \brief Define natural (C-language) or Fortran ordering of array dimensions.
 *  A RankZero array does not have an ordering.
 */
enum ArrayOrder { NaturalOrder , FortranOrder , RankZero };

//----------------------------------------------------------------------
/** \class Array
 *  \brief Multidimensional array index mapping into a contiguous storage.
 *
 *  \tparam Scalar Type of the array's members.
 *  \tparam Order  How to order the array's multidimensions and multi-indices,
 *                 either NaturalOrder or FortranOrder.
 *  \tparam Tag1   Array dimension tag for the first dimension.
 *  \tparam Tag2   Array dimension tag for the second dimension.
 *  \tparam Tag8   Array dimension tag for the eighth dimension.
 *
 *  Members of an Array type include the following.
 *  - Types:
 *    -# value_type ;   Type for the array's members.
 *    -# Tag<Ordinal>::type ; TagN where N = Ordinal + 1
 *
 *  - Enumerations for compile-time traits:
 *    -# Rank ;       Rank of the array, the number of multi-indices.
 *    -# Natural ;    If natural (a.k.a. C) multi-index ordering is used.
 *    -# Reverse ;    If reverse (a.k.a. Fortran) multi-index ordering is used.
 *    -# Contiguous ; If members are stored in contiguous memory.
 *
 *  - Functions for run-time traits:
 *    -# unsigned rank() const ;
 *    -# bool natural() const ; 
 *    -# bool reverse() const ; 
 *    -# bool contiguous() const ; 
 *    -# unsigned dimension<Ordinal>() const ;
 *    -# unsigned dimension( unsigned ordinal ) const ;
 *    -# void dimensions( std::vector<unsigned> & ) const ;
 *
 *  - Member data access functions:
 *    -# value_type * contiguous_data() const ; If data is contiguous
 *    -# value_type & operator[]( size_type ) const ;\n
 *       Member access by fully-ordered offset.
 *    -# value_type & operator()( const unsigned i1 , ... ) const ;\n
 *       Member access by multi-index of the proper rank.
 *
 *  - Constructors and assignment operators:
 *    -# Default constructor for empty array.
 *    -# Standard copy constructor and assignment operator
 *       generate a new shared view into existing storage.
 *    -# Copy constructor and assignment operator for compatible arrays
 *       of the reversed type.  A compatible reversed array type
 *       exchanges the NaturalOrder for FortranOrder and reverses
 *       the dimension tags to generate a shared compatible view
 *       of the existing storage.
 *    -# Contructor accepting pointer to contiguous storage and
 *       multidimension argument list to generate an array view of the storage.
 *    -# Contructor accepting pointer to contiguous storage and
 *       multidimension array argument to generate an array view of the storage.
 *
 *  - Truncation function:
 *    -# Array::Truncate truncate( const unsigned i ) const ;\n
 *       Generates an array view into the existing array which
 *       is offset by "i" in the slowest changing dimension
 *       (first Natural or last Fortran dimnension) and has the
 *       slowest changing dimension truncated from the returned array.
 */
template< typename Scalar , ArrayOrder Order , 
          class Tag1 = void , class Tag2 = void ,
          class Tag3 = void , class Tag4 = void ,
          class Tag5 = void , class Tag6 = void ,
          class Tag7 = void , class Tag8 = void >
class Array ;

//----------------------------------------------------------------------

template< class ArrayType , class Tag > struct ArrayAppend ;

//----------------------------------------------------------------------
/** \class ArrayDimTag
 *  \brief Virtual base class for array dimension tags.
 *  A derived array dimension tag class, for example
 *  'class MyTag : public phdmesh::ArrayDimTag', must provide
 *  a static method 'const MyTag & tag();' that returns a
 *  singleton for the derived class.  For example,
 *
 *  const MyTag & MyTag::tag() { static const MyTag t ; return t ; }
 */
struct ArrayDimTag {

  /** Name of the tag, typically the name of the derived class */
  virtual const char * name() const = 0 ;

  /** Given a dimension and index produce a string for output. */
  virtual std::string to_string( unsigned dimension ,
                                 unsigned index ) const ;

  /** Given a dimension and input strige produce an index */
  virtual unsigned to_index( unsigned , const std::string & ) const ; 
protected:
  virtual ~ArrayDimTag();
  ArrayDimTag() {}
private:
  ArrayDimTag( const ArrayDimTag & );
  ArrayDimTag & operator = ( const ArrayDimTag & );
};

/** \class  ArrayDimension
 *  \brief  An anonymous array dimension tag,
 *          which is NOT the recommended usage.
 */
struct ArrayDimension : public ArrayDimTag {

  const char * name() const ;

  static const ArrayDimension & tag();

private:
  ~ArrayDimension();
  ArrayDimension() {}
  ArrayDimension( const ArrayDimension & );
  ArrayDimension & operator = ( const ArrayDimension & );
};

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
/** \cond */
//----------------------------------------------------------------------
//----------------------------------------------------------------------

#include <util/ArrayPrivate.hpp>

namespace phdmesh {

// Rank 8 array:

template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class Array
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 8 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 , const unsigned i8 ) const
    { return m_ptr[
        array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 ,
         const unsigned n7 , const unsigned n8 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::
        assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 , n7 , n8 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 7:

template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
class Array<Scalar,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 7 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 7 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 ) const
    { return m_ptr[
        array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 ,
         const unsigned n7 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>::
        assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 , n7 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,void>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 6:

template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
class Array<Scalar,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 6 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 6 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ) const
    { return m_ptr[
        array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 , const unsigned n6 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>::
        assign( m_stride , n1 , n2 , n3 , n4 , n5 , n6 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,void,void>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 5:

template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
class Array<Scalar,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 5 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 5 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4,i5) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 ,
         const unsigned n5 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>::
        assign( m_stride , n1 , n2 , n3 , n4 , n5 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,Tag5,void,void,void>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 4:

template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 , class Tag4 >
class Array<Scalar,array_order,Tag1,Tag2,Tag3,Tag4,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 4 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,Tag2,Tag3,Tag4,void,void,void,void>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 4 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2,i3,i4) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 , const unsigned n4 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,void,void,void,void>::
        assign( m_stride , n1 , n2 , n3 , n4 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,Tag4,void,void,void,void>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 3:

template< typename Scalar , ArrayOrder array_order ,
          class Tag1 , class Tag2 , class Tag3 >
class Array<Scalar,array_order,Tag1,Tag2,Tag3,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 3 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,Tag2,Tag3,void,void,void,void,void>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 3 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2,i3) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ,
         const unsigned n1 , const unsigned n2 ,
         const unsigned n3 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,void,void,void,void,void>::
        assign( m_stride , n1 , n2 , n3 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,Tag3,void,void,void,void,void>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 2:

template< typename Scalar , ArrayOrder array_order , class Tag1 , class Tag2 >
class Array<Scalar,array_order,Tag1,Tag2,void,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 2 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,Tag2,void,void,void,void,void,void>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 2 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1,i2) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr , const unsigned n1 , const unsigned n2 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,void,void,void,void,void,void>::
        assign( m_stride , n1 , n2 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,Tag2,void,void,void,void,void,void>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 1:

template< typename Scalar , ArrayOrder array_order , class Tag1 >
class Array<Scalar,array_order,Tag1,void,void,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 1 };
  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------
  // ArrayType::Tag<K>::type 
  //
  template < unsigned ordinal >
  struct Tag { typedef typename ArrayTagAt<Array,ordinal>::type type ; };

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return array_dim_tags<Tag1,void,void,void,void,void,void,void>()[ordinal];
    }

  //----------------------------------

  size_type size() const { return m_stride[ Rank - 1 ]; }

  // ArrayType::dimension<K>();
  template < unsigned ordinal > unsigned dimension() const
    {
      array_check_ordinal_is_less<ordinal,Rank>();
      return ArrayStrideDim<array_order,Rank,ordinal>::dimension(m_stride);
    }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( Rank , ordinal );
      return ArrayStrideDim<array_order,Rank>::dimension(m_stride,ordinal);
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( Rank );
      for ( unsigned i = 0 ; i < Rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  /** \brief Access member via Rank 1 multi-index */
  value_type & operator()( const unsigned i1 ) const
    { return m_ptr[ array_offset<array_order,Rank>(m_stride,i1) ]; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) { Copy<Rank>( m_stride , (size_type) 0 ); }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ) { Copy<Rank>( m_stride , rhs.m_stride ); }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      Copy<Rank>( m_stride , rhs.m_stride );
      return *this ;
    }

  //----------------------------------
  // Truncated view

  typedef typename ArrayTruncate<Array>::type TruncateType ;

  TruncateType truncate( const unsigned i ) const
    {
      TruncateType tmp ;
      tmp.m_ptr = m_ptr + m_stride[ Rank - 2 ] * i ;
      Copy<Rank-1>( tmp.m_stride , m_stride );
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr , const unsigned n1 )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,void,void,void,void,void,void,void>::
        assign( m_stride , n1 );
    }

  Array( value_type * arg_ptr , const unsigned * const dims )
    : m_ptr( arg_ptr )
    {
      Array<void,array_order,Tag1,void,void,void,void,void,void,void>::
        assign( m_stride , dims );
    }

protected:

  Scalar  * m_ptr ;
  size_type m_stride[ Rank ];

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
// Rank 0:

template< typename Scalar >
class Array<Scalar,RankZero,void,void,void,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Rank       = 0 };
  enum { Natural    = false };
  enum { Reverse    = false };
  enum { Contiguous = true };

  unsigned rank()   const { return Rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  size_type size() const { return m_ptr ? 1 : 0 ; }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via Rank 0 multi-index */
  value_type & operator()() const { return *m_ptr ; }

  //----------------------------------
  // Required constructors and assignment operators:

  Array() : m_ptr(NULL) {}

  Array( const Array & rhs ) : m_ptr( rhs.m_ptr ) {}

  Array & operator = ( const Array & rhs )
    { m_ptr = rhs.m_ptr ; return *this ; }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * arg_ptr ) : m_ptr( arg_ptr ) {}

protected:

  Scalar * m_ptr ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Runtime rank and tag information:

template< typename Scalar , ArrayOrder array_order >
class Array<Scalar,array_order,void,void,void,void,void,void,void,void>
{
public:

  typedef Scalar              value_type ;
  typedef unsigned            size_type ;
  typedef const ArrayDimTag * tag_type ;

  //----------------------------------

  enum { Natural    = NaturalOrder == array_order };
  enum { Reverse    = FortranOrder == array_order };
  enum { Contiguous = true };

  unsigned rank()   const { return m_rank ; }
  bool natural()    const { return Natural ; }
  bool reverse()    const { return Reverse ; }
  bool contiguous() const { return Contiguous ; }

  //----------------------------------

  typedef typename ArrayReverse< Array >::type ReverseType ;

  //----------------------------------

  tag_type tag( const unsigned ordinal ) const
    {
      array_check_ordinal( m_rank , ordinal );
      const int i = Natural ? ( m_rank - 1 ) - ordinal : ordinal ;
      return m_tag[i];
    }

  //----------------------------------

  size_type size() const { return m_stride[ m_rank - 1 ]; }

  // ArrayType::dimension(K);
  unsigned dimension( const unsigned ordinal ) const
    {
      array_check_ordinal( m_rank , ordinal );
      const int i = Natural ? ( m_rank - 1 ) - ordinal : ordinal ;
      return i ? m_stride[i] / m_stride[i-1] : m_stride[i] ;
    }

  void dimensions( std::vector<unsigned> & n )
    {
      n.resize( m_rank );
      for ( unsigned i = 0 ; i < m_rank ; ++i ) { n[i] = dimension(i); }
    }

  //----------------------------------
  /** \brief Access member data */
  value_type * contiguous_data() const { return m_ptr ; }

  /** \brief Access member via full ordering of members. */
  value_type & operator[]( size_type i ) const
    {
      ARRAY_CHECK( array_check_offset(size(),i) );
      return m_ptr[ i ];
    }

  //----------------------------------
  /** \brief Access member via Rank 8 multi-index */
  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 , const unsigned i8 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 8 ) );
      return m_ptr[
        array_offset<array_order,8>(m_stride,i1,i2,i3,i4,i5,i6,i7,i8) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ,
                           const unsigned i7 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 7 ) );
      return m_ptr[
        array_offset<array_order,7>(m_stride,i1,i2,i3,i4,i5,i6,i7) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 , const unsigned i6 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 6 ) );
      return m_ptr[ array_offset<array_order,6>(m_stride,i1,i2,i3,i4,i5,i6) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ,
                           const unsigned i5 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 5 ) );
      return m_ptr[ array_offset<array_order,5>(m_stride,i1,i2,i3,i4,i5) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 , const unsigned i4 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 4 ) );
      return m_ptr[ array_offset<array_order,4>(m_stride,i1,i2,i3,i4) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ,
                           const unsigned i3 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 3 ) );
      return m_ptr[ array_offset<array_order,3>(m_stride,i1,i2,i3) ];
    }

  value_type & operator()( const unsigned i1 , const unsigned i2 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 2 ) );
      return m_ptr[ array_offset<array_order,2>(m_stride,i1,i2) ];
    }

  value_type & operator()( const unsigned i1 ) const
    {
      ARRAY_CHECK( array_check_rank( m_rank , 1 ) );
      return m_ptr[ array_offset<array_order,1>(m_stride,i1) ];
    }

  //----------------------------------
  // Required constructors and assignment operators:

  Array()
    : m_ptr(NULL), m_rank(0)
    {
      Copy<8>( m_stride , (size_type) 0 );
      Copy<8>( m_tag , (tag_type) NULL );
    }

  Array( const Array & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
    }

  Array & operator = ( const Array & rhs )
    {
      m_ptr = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
      return *this ;
    }

  Array( const ReverseType & rhs )
    : m_ptr( rhs.m_ptr ), m_rank( rhs.m_rank )
    {
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
    }

  Array & operator = ( const ReverseType & rhs )
    {
      m_ptr = rhs.m_ptr ;
      m_rank = rhs.m_rank ;
      Copy<8>( m_stride , rhs.m_stride );
      Copy<8>( m_tag , rhs.m_tag );
      return *this ;
    }

  //----------------------------------

  template< ArrayOrder order ,
            class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
            class Tag5 , class Tag6 , class Tag7 , class Tag8 >
  Array(
    const Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & rhs )
  : m_ptr( rhs.m_ptr ), m_rank( 0 )
  {
    typedef Array<Scalar,order,Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> a_t ;
    enum { inRank    = a_t::Rank };
    enum { inNatural = a_t::Natural };
    m_rank = inRank ;
    Copy< inRank >(     m_stride , rhs.m_stride );
    Copy< 8 - inRank >( m_stride + inRank , (size_type) 0 );
    unsigned i = 0 ;
    if ( inNatural ) {
      for ( ; i < inRank ; ++i ) { m_tag[i] = rhs.tag((inRank-1)-i); }
    }
    else {
      for ( ; i < inRank ; ++i ) { m_tag[i] = rhs.tag(i); }
    }
    for ( ; i < 8 ; ++i ) { m_tag[i] = NULL ; }
  }

  //----------------------------------
  // Truncated view

  Array truncate( const unsigned i ) const
    {
      Array tmp ;
      if ( 1 < m_rank ) {
        tmp.m_ptr  = m_ptr + m_stride[ m_rank - 2 ] * i ;
        tmp.m_rank = m_rank - 1 ;
        unsigned k ;
        for ( k = 0 ; k < m_rank - 1 ; ++k ) { tmp.m_stride[i] = m_stride[i] ; }
        for (       ; k < 8          ; ++k ) { tmp.m_stride[i] = 0 ; }
        for ( k = 0 ; k < m_rank - 1 ; ++k ) { tmp.m_tag[i] = m_tag[i] ; }
        for (       ; k < 8          ; ++k ) { tmp.m_tag[i] = NULL ; }
      }
      return tmp ;
    }

  //----------------------------------
  // Class specific constructors:

  Array( value_type * ptr ,
         const unsigned rank ,
         const unsigned * const dims ,
         const tag_type  * const tags )
    : m_ptr( ptr ), m_rank( rank )
    {
      if ( Natural ) {
        size_type n = 1 ;
        unsigned i ;
        for ( i = 0 ; i < rank ; ++i ) { m_stride[i] = n *= dims[(rank-1)-i]; }
        for (       ; i < 8    ; ++i ) { m_stride[i] = 0 ; }
        for ( i = 0 ; i < rank ; ++i ) { m_tag[i] = tags[(rank-1)-i]; }
        for (       ; i < 8    ; ++i ) { m_tag[i] = NULL ; }
      }
      else {
        size_type n = 1 ;
        unsigned i ;
        for ( i = 0 ; i < rank ; ++i ) { m_stride[i] = n *= dims[i] ; }
        for (       ; i < 8    ; ++i ) { m_stride[i] = 0 ; }
        for ( i = 0 ; i < rank ; ++i ) { m_tag[i] = tags[i]; }
        for (       ; i < 8    ; ++i ) { m_tag[i] = NULL ; }
      }
    }

protected:

  Scalar  * m_ptr ;
  unsigned  m_rank ;
  size_type m_stride[8];
  tag_type  m_tag[8] ;

  template< typename , ArrayOrder ,
            class , class , class , class ,
            class , class , class , class >
  friend class phdmesh::Array ;
};

//----------------------------------------------------------------------
//----------------------------------------------------------------------

/** \endcond */

}

#undef ARRAY_CHECK

#endif

