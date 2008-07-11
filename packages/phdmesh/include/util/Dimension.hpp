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
 * @date   May 2008
 */

#ifndef util_Dimension_hpp
#define util_Dimension_hpp

#include <iosfwd>
#include <string>

namespace phdmesh {

//----------------------------------------------------------------------

class DimensionTag {
public:

  virtual ~DimensionTag() {}

  virtual const char * name() const = 0 ;

  /** Generate a descriptive text label for a dimension and index.
   *  Throw an exception for invalid arguments.
   */
  virtual std::string to_string( unsigned size , unsigned index ) const ;

  /** Generate an index from a dimension and text label.
   *  Throw an exception for invalid arguments.
   */
  virtual unsigned to_index( unsigned size , const std::string & ) const ;

  static const DimensionTag * descriptor() { return NULL ; }
};

//----------------------------------------------------------------------
/** @class  Dimension
 *  @brief  Multidimensional array specification and offset computation.
 *
 *  Given 'Dimension<...> & map'  the operators and methods are
 *  as follows.
 *
 *  Map an index to an offset:
 *
 *    unsigned offset = map( i1 , i2 , ... );
 *   
 *  Query if an index is valid:
 *
 *    bool ok = map.valid( i1 , i2 , ... );
 *
 *  Query the size, per dimension:
 *
 *    unsigned n1 , n2 , ... ;
 *    map.size( n1 , n2 , ... );
 *
 *  Query the total size, the product of the dimensions:
 *
 *    unsigned total_size = map.size();
 *
 */
template< class Tag1 = DimensionTag ,
          class Tag2 = DimensionTag ,
          class Tag3 = DimensionTag ,
          class Tag4 = DimensionTag ,
          class Tag5 = DimensionTag ,
          class Tag6 = DimensionTag ,
          class Tag7 = DimensionTag ,
          class Tag8 = DimensionTag >
class Dimension ;

template< class Dim , class Tag > class DimensionAppend ;
template< class Dim >             class DimensionTruncate ;

enum { DimensionMaximum = 8 };

void     stride_verify( unsigned , const unsigned *, const char * );
void     stride_copy( unsigned , unsigned *, const unsigned *, const char * );
unsigned stride_map(  unsigned , const unsigned * , const unsigned * );
void     stride_inv(  unsigned , const unsigned * , unsigned , unsigned * );
void     stride_size( unsigned , const unsigned * , unsigned * );
unsigned stride_size( unsigned , const unsigned * );

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
class Dimension< DimensionTag, DimensionTag,
                 DimensionTag, DimensionTag,
                 DimensionTag, DimensionTag,
                 DimensionTag, DimensionTag> {
public:
  enum { NumDim = 0 };

//  unsigned stride[ NumDim ];

  Dimension( const unsigned * const , const char * ) {}
};

//----------------------------------------------------------------------

template< class Tag1 >
class Dimension< Tag1 ,
                 DimensionTag ,
                 DimensionTag , DimensionTag ,
                 DimensionTag , DimensionTag ,
                 DimensionTag , DimensionTag > {
public:
  enum { NumDim = 1 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()( unsigned i1 ) const { return i1 ; }

  /** Map from offset to indices */
  void inverse( unsigned offset , unsigned & i1 ) const
    { i1 = offset ; }

  /** Query if the indices are valid */
  bool valid( unsigned i1 ) const { return i1 < stride[0] ; }

  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 ) const { n1 = stride[0] ; }

  explicit Dimension( unsigned n1 ) { stride[0] = n1 ; }

  Dimension()
    { stride[0] = 0 ; }

  Dimension( const Dimension & rhs )
    { stride[0] = rhs.stride[0] ; }

  Dimension & operator = ( const Dimension & rhs )
    { stride[0] = rhs.stride[0] ; return *this ; }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
    }

  //--------------------------------

  template< class T >
  explicit Dimension( const Dimension<Tag1,T> & rhs )
    { stride[0] = rhs.stride[0] ; }

  Dimension( const Dimension<> & , unsigned n )
    { stride[0] = n ; }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 >
class Dimension< Tag1 , Tag2 ,
                 DimensionTag , DimensionTag ,
                 DimensionTag , DimensionTag ,
                 DimensionTag , DimensionTag > {
public:
  enum { NumDim = 2 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()( unsigned i1 , unsigned i2 ) const
    { return i1 + i2 * stride[0] ; }

  /** Map from offset to indices */
  void inverse( unsigned offset , unsigned & i1 , unsigned & i2 ) const
    {
      i2 = offset / stride[0] ; offset %= stride[0] ;
      i1 = offset ;
    }

  /** Query if the indices are valid */
  bool valid( unsigned i1 , unsigned i2 ) const
    {
      return i1               < stride[0] &&
             i2 * stride[0] < stride[1] ;
    }


  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 , unsigned & n2 ) const
    {
      n1 = stride[0] ;
      n2 = stride[1] / stride[0] ;
    }

  Dimension( unsigned n1 , unsigned n2 )
    {
      stride[1] = n2 * (
      stride[0] = n1 );
    }

  Dimension()
    {
      stride[1] = 0 ;
      stride[0] = 0 ;
    }

  Dimension( const Dimension & rhs )
    {
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
    }

  Dimension & operator = ( const Dimension & rhs )
    {
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
      return *this ;
    }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
      stride[1] = n[1] ;
    }

  //--------------------------------

  template< class T >
  explicit Dimension( const Dimension<Tag1,Tag2,T> & rhs )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
    }

  Dimension( const Dimension<Tag1> & rhs , unsigned n )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = n * stride[0] ;
    }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 >
class Dimension< Tag1 , Tag2 , Tag3 ,
                 DimensionTag ,
                 DimensionTag , DimensionTag ,
                 DimensionTag , DimensionTag > {
public:
  enum { NumDim = 3 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 ) const
    { return i1 + i2 * stride[0] + i3 * stride[1] ; }

  /** Map from offset to indices */
  void inverse( unsigned offset ,
                unsigned & i1 , unsigned & i2 , unsigned & i3 ) const
    {
      i3 = offset / stride[1] ; offset %= stride[1] ;
      i2 = offset / stride[0] ; offset %= stride[0] ;
      i1 = offset ;
    }

  /** Query if the indices are valid */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 ) const
    {
      return i1               < stride[0] &&
             i2 * stride[0] < stride[1] &&
             i3 * stride[1] < stride[2] ;
    }


  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 , unsigned & n2 , unsigned & n3 ) const
    {
      n1 = stride[0] ;
      n2 = stride[1] / stride[0] ;
      n3 = stride[2] / stride[1] ;
    }

  Dimension( unsigned n1 , unsigned n2 , unsigned n3 )
    {
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 ));
    }

  Dimension()
    {
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = 0 ;
    }

  Dimension( const Dimension & rhs )
    {
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
    }

  Dimension & operator = ( const Dimension & rhs )
    {
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
      return *this ;
    }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
      stride[1] = n[1] ;
      stride[2] = n[2] ;
    }

  //--------------------------------

  template< class T >
  explicit Dimension( const Dimension<Tag1,Tag2,Tag3,T> & rhs )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
    }

  Dimension( const Dimension<Tag1,Tag2> & rhs , unsigned n )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = n * stride[1] ;
    }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
class Dimension< Tag1 , Tag2 , Tag3 , Tag4 ,
                 DimensionTag , DimensionTag ,
                 DimensionTag , DimensionTag > {
public:
  enum { NumDim = 4 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ) const
    {
      return i1 + i2 * stride[0] + i3 * stride[1] + i4 * stride[2] ;
    }

  /** Map from offset to indices */
  void inverse( unsigned offset ,
                unsigned & i1 , unsigned & i2 , unsigned & i3 , unsigned & i4 )
    const
    {
      i4 = offset / stride[2] ; offset %= stride[2] ;
      i3 = offset / stride[1] ; offset %= stride[1] ;
      i2 = offset / stride[0] ; offset %= stride[0] ;
      i1 = offset ;
    }

  /** Query if the indices are valid */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ) const
    {
      return i1               < stride[0] &&
             i2 * stride[0] < stride[1] &&
             i3 * stride[1] < stride[2] &&
             i4 * stride[2] < stride[3] ;
    }


  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 , unsigned & n2 , unsigned & n3 , unsigned & n4 )
    const
    {
      n1 = stride[0] ;
      n2 = stride[1] / stride[0] ;
      n3 = stride[2] / stride[1] ;
      n4 = stride[3] / stride[2] ;
    }

  Dimension( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 )
    {
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 )));
    }

  Dimension()
    {
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = 0 ;
    }

  Dimension( const Dimension & rhs )
    {
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
    }

  Dimension & operator = ( const Dimension & rhs )
    {
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
      return *this ;
    }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
      stride[1] = n[1] ;
      stride[2] = n[2] ;
      stride[3] = n[3] ;
    }

  //--------------------------------

  template< class T >
  explicit Dimension( const Dimension<Tag1,Tag2,Tag3,Tag4,T> & rhs )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
    }

  Dimension( const Dimension<Tag1,Tag2,Tag3> & rhs , unsigned n )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = n * stride[2] ;
    }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 >
class Dimension< Tag1 , Tag2 , Tag3 , Tag4 , Tag5 ,
                 DimensionTag , DimensionTag , DimensionTag > {
public:
  enum { NumDim = 5 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 ) const
    {
      return i1 + i2 * stride[0] + i3 * stride[1] + i4 * stride[2] +
                  i5 * stride[3] ;
    }

  /** Map from offset to indices */
  void inverse( unsigned offset ,
                unsigned & i1 , unsigned & i2 , unsigned & i3 , unsigned & i4 ,
                unsigned & i5 ) const
    {
      i5 = offset / stride[3] ; offset %= stride[3] ;
      i4 = offset / stride[2] ; offset %= stride[2] ;
      i3 = offset / stride[1] ; offset %= stride[1] ;
      i2 = offset / stride[0] ; offset %= stride[0] ;
      i1 = offset ;
    }

  /** Query if the indices are valid */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 ) const
    {
      return i1               < stride[0] &&
             i2 * stride[0] < stride[1] &&
             i3 * stride[1] < stride[2] &&
             i4 * stride[2] < stride[3] &&
             i5 * stride[3] < stride[4] ;
    }


  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 , unsigned & n2 , unsigned & n3 , unsigned & n4 ,
             unsigned & n5 ) const
    {
      n1 = stride[0] ;
      n2 = stride[1] / stride[0] ;
      n3 = stride[2] / stride[1] ;
      n4 = stride[3] / stride[2] ;
      n5 = stride[4] / stride[3] ;
    }

  Dimension( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
             unsigned n5 )
    {
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 ))));
    }

  Dimension()
    {
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = 0 ;
    }

  Dimension( const Dimension & rhs )
    {
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
    }

  Dimension & operator = ( const Dimension & rhs )
    {
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
      return *this ;
    }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
      stride[1] = n[1] ;
      stride[2] = n[2] ;
      stride[3] = n[3] ;
      stride[4] = n[4] ;
    }

  //--------------------------------

  template< class T >
  explicit Dimension(
    const Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,T> & rhs )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
      stride[4] = rhs.stride[4] ;
    }

  Dimension( const Dimension<Tag1,Tag2,Tag3,Tag4> & rhs , unsigned n )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
      stride[4] = n * stride[3] ;
    }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 >
class Dimension< Tag1 , Tag2 , Tag3 , Tag4 ,
                 Tag5 , Tag6 , DimensionTag , DimensionTag > {
public:
  enum { NumDim = 6 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 ) const
    {
      return i1 + i2 * stride[0] + i3 * stride[1] + i4 * stride[2] +
                  i5 * stride[3] + i6 * stride[4] ;
    }

  /** Map from offset to indices */
  void inverse( unsigned offset ,
                unsigned & i1 , unsigned & i2 , unsigned & i3 , unsigned & i4 ,
                unsigned & i5 , unsigned & i6 ) const
    {
      i6 = offset / stride[4] ; offset %= stride[4] ;
      i5 = offset / stride[3] ; offset %= stride[3] ;
      i4 = offset / stride[2] ; offset %= stride[2] ;
      i3 = offset / stride[1] ; offset %= stride[1] ;
      i2 = offset / stride[0] ; offset %= stride[0] ;
      i1 = offset ;
    }

  /** Query if the indices are valid */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 ) const
    {
      return i1               < stride[0] &&
             i2 * stride[0] < stride[1] &&
             i3 * stride[1] < stride[2] &&
             i4 * stride[2] < stride[3] &&
             i5 * stride[3] < stride[4] &&
             i6 * stride[4] < stride[5] ;
    }


  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 , unsigned & n2 , unsigned & n3 , unsigned & n4 ,
             unsigned & n5 , unsigned & n6 ) const
    {
      n1 = stride[0] ;
      n2 = stride[1] / stride[0] ;
      n3 = stride[2] / stride[1] ;
      n4 = stride[3] / stride[2] ;
      n5 = stride[4] / stride[3] ;
      n6 = stride[5] / stride[4] ;
    }

  Dimension( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
             unsigned n5 , unsigned n6 )
    {
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 )))));
    }

  Dimension()
    {
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = 0 ;
    }

  Dimension( const Dimension & rhs )
    {
      stride[5] = rhs.stride[5] ;
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
    }

  Dimension & operator = ( const Dimension & rhs )
    {
      stride[5] = rhs.stride[5] ;
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
      return *this ;
    }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
      stride[1] = n[1] ;
      stride[2] = n[2] ;
      stride[3] = n[3] ;
      stride[4] = n[4] ;
      stride[5] = n[5] ;
    }

  //--------------------------------

  template< class T >
  explicit Dimension(
    const Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,T> & rhs )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
      stride[4] = rhs.stride[4] ;
      stride[5] = rhs.stride[5] ;
    }

  Dimension(
    const Dimension<Tag1,Tag2,Tag3,Tag4,Tag5> & rhs , unsigned n )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
      stride[4] = rhs.stride[4] ;
      stride[5] = n * stride[4] ;
    }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
class Dimension< Tag1 , Tag2 , Tag3 , Tag4 ,
                 Tag5 , Tag6 , Tag7 , DimensionTag > {
public:
  enum { NumDim = 7 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 , unsigned i7 ) const
    {
      return i1 + i2 * stride[0] + i3 * stride[1] + i4 * stride[2] +
                  i5 * stride[3] + i6 * stride[4] + i7 * stride[5] ;
    }

  /** Map from offset to indices */
  void inverse( unsigned offset ,
                unsigned & i1 , unsigned & i2 , unsigned & i3 , unsigned & i4 ,
                unsigned & i5 , unsigned & i6 , unsigned & i7 ) const
    {
      i7 = offset / stride[5] ; offset %= stride[5] ;
      i6 = offset / stride[4] ; offset %= stride[4] ;
      i5 = offset / stride[3] ; offset %= stride[3] ;
      i4 = offset / stride[2] ; offset %= stride[2] ;
      i3 = offset / stride[1] ; offset %= stride[1] ;
      i2 = offset / stride[0] ; offset %= stride[0] ;
      i1 = offset ;
    }

  /** Query if the indices are valid */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 , unsigned i7 ) const
    {
      return i1               < stride[0] &&
             i2 * stride[0] < stride[1] &&
             i3 * stride[1] < stride[2] &&
             i4 * stride[2] < stride[3] &&
             i5 * stride[3] < stride[4] &&
             i6 * stride[4] < stride[5] &&
             i7 * stride[5] < stride[6] ;
    }


  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 , unsigned & n2 , unsigned & n3 , unsigned & n4 ,
             unsigned & n5 , unsigned & n6 , unsigned & n7 ) const
    {
      n1 = stride[0] ;
      n2 = stride[1] / stride[0] ;
      n3 = stride[2] / stride[1] ;
      n4 = stride[3] / stride[2] ;
      n5 = stride[4] / stride[3] ;
      n6 = stride[5] / stride[4] ;
      n7 = stride[6] / stride[5] ;
    }

  Dimension( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
             unsigned n5 , unsigned n6 , unsigned n7 )
    {
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 ))))));
    }

  Dimension()
    {
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = 0 ;
    }

  Dimension( const Dimension & rhs )
    {
      stride[6] = rhs.stride[6] ;
      stride[5] = rhs.stride[5] ;
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
    }

  Dimension & operator = ( const Dimension & rhs )
    {
      stride[6] = rhs.stride[6] ;
      stride[5] = rhs.stride[5] ;
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
      return *this ;
    }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
      stride[1] = n[1] ;
      stride[2] = n[2] ;
      stride[3] = n[3] ;
      stride[4] = n[4] ;
      stride[5] = n[5] ;
      stride[6] = n[6] ;
    }

  //--------------------------------

  template< class T >
  explicit Dimension(
    const Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,T> & rhs )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
      stride[4] = rhs.stride[4] ;
      stride[5] = rhs.stride[5] ;
      stride[6] = rhs.stride[6] ;
    }

  Dimension(
    const Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6> & rhs ,
    unsigned n )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
      stride[4] = rhs.stride[4] ;
      stride[5] = rhs.stride[5] ;
      stride[6] = n * stride[5] ;
    }
};

//----------------------------------------------------------------------

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
class Dimension {
public:
  enum { NumDim = 8 };

  unsigned stride[ NumDim ];

  /** Map dimension indices to an offset. */
  unsigned operator()(
    unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
    unsigned i5 , unsigned i6 , unsigned i7 , unsigned i8 ) const
    {
      return i1 + i2 * stride[0] + i3 * stride[1] + i4 * stride[2] +
                  i5 * stride[3] + i6 * stride[4] + i7 * stride[5] +
                  i8 * stride[6] ;
    }

  /** Map from offset to indices */
  void inverse( unsigned offset ,
                unsigned & i1 , unsigned & i2 , unsigned & i3 , unsigned & i4 ,
                unsigned & i5 , unsigned & i6 , unsigned & i7 , unsigned & i8 )
    const
    {
      i8 = offset / stride[6] ; offset %= stride[6] ;
      i7 = offset / stride[5] ; offset %= stride[5] ;
      i6 = offset / stride[4] ; offset %= stride[4] ;
      i5 = offset / stride[3] ; offset %= stride[3] ;
      i4 = offset / stride[2] ; offset %= stride[2] ;
      i3 = offset / stride[1] ; offset %= stride[1] ;
      i2 = offset / stride[0] ; offset %= stride[0] ;
      i1 = offset ;
    }

  /** Query if the indices are valid */
  bool valid( unsigned i1 , unsigned i2 , unsigned i3 , unsigned i4 ,
              unsigned i5 , unsigned i6 , unsigned i7 , unsigned i8 ) const
    {
      return i1               < stride[0] &&
             i2 * stride[0] < stride[1] &&
             i3 * stride[1] < stride[2] &&
             i4 * stride[2] < stride[3] &&
             i5 * stride[3] < stride[4] &&
             i6 * stride[4] < stride[5] &&
             i7 * stride[5] < stride[6] &&
             i8 * stride[6] < stride[7] ;
    }


  /** Query offset upper bound = the product of dimensions */
  unsigned size() const { return stride[ NumDim - 1 ]; }

  /** Query dimensions */
  void size( unsigned & n1 , unsigned & n2 , unsigned & n3 , unsigned & n4 ,
             unsigned & n5 , unsigned & n6 , unsigned & n7 , unsigned & n8 )
    const
    {
      n1 = stride[0] ;
      n2 = stride[1] / stride[0] ;
      n3 = stride[2] / stride[1] ;
      n4 = stride[3] / stride[2] ;
      n5 = stride[4] / stride[3] ;
      n6 = stride[5] / stride[4] ;
      n7 = stride[6] / stride[5] ;
      n8 = stride[7] / stride[6] ;
    }

  Dimension( unsigned n1 , unsigned n2 , unsigned n3 , unsigned n4 ,
             unsigned n5 , unsigned n6 , unsigned n7 , unsigned n8 )
    {
      stride[7] = n8 * (
      stride[6] = n7 * (
      stride[5] = n6 * (
      stride[4] = n5 * (
      stride[3] = n4 * (
      stride[2] = n3 * (
      stride[1] = n2 * (
      stride[0] = n1 )))))));
    }

  Dimension()
    {
      stride[7] = 0 ;
      stride[6] = 0 ;
      stride[5] = 0 ;
      stride[4] = 0 ;
      stride[3] = 0 ;
      stride[2] = 0 ;
      stride[1] = 0 ;
      stride[0] = 0 ;
    }

  Dimension( const Dimension & rhs )
    {
      stride[7] = rhs.stride[7] ;
      stride[6] = rhs.stride[6] ;
      stride[5] = rhs.stride[5] ;
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
    }

  Dimension & operator = ( const Dimension & rhs )
    {
      stride[7] = rhs.stride[7] ;
      stride[6] = rhs.stride[6] ;
      stride[5] = rhs.stride[5] ;
      stride[4] = rhs.stride[4] ;
      stride[3] = rhs.stride[3] ;
      stride[2] = rhs.stride[2] ;
      stride[1] = rhs.stride[1] ;
      stride[0] = rhs.stride[0] ;
      return *this ;
    }

  Dimension( const unsigned * const n , const char * error_source )
    {
      stride_verify( NumDim , n , error_source );
      stride[0] = n[0] ;
      stride[1] = n[1] ;
      stride[2] = n[2] ;
      stride[3] = n[3] ;
      stride[4] = n[4] ;
      stride[5] = n[5] ;
      stride[6] = n[6] ;
      stride[7] = n[7] ;
    }

  //--------------------------------

  struct Truncate {
    typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> type ;
  };

  Dimension(
    const Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> & rhs ,
    unsigned n )
    {
      stride[0] = rhs.stride[0] ;
      stride[1] = rhs.stride[1] ;
      stride[2] = rhs.stride[2] ;
      stride[3] = rhs.stride[3] ;
      stride[4] = rhs.stride[4] ;
      stride[5] = rhs.stride[5] ;
      stride[6] = rhs.stride[6] ;
      stride[7] = n * stride[6] ;
    }
};

//----------------------------------------------------------------------

template< class Tag1 >
struct DimensionAppend< Dimension<>, Tag1>
{ typedef Dimension<Tag1> type ; };

template< class Tag1 , class Tag2 >
struct DimensionAppend< Dimension<Tag1>, Tag2>
{ typedef Dimension<Tag1,Tag2> type ; };

template< class Tag1 , class Tag2 , class Tag3 >
struct DimensionAppend< Dimension<Tag1,Tag2>, Tag3>
{ typedef Dimension<Tag1,Tag2,Tag3> type ; };

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct DimensionAppend< Dimension<Tag1,Tag2,Tag3>, Tag4>
{ typedef Dimension<Tag1,Tag2,Tag3,Tag4> type ; };

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 , class Tag5 >
struct DimensionAppend< Dimension<Tag1,Tag2,Tag3,Tag4>, Tag5>
{ typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5> type ; };

template< class Tag1 , class Tag2 , class Tag3 ,
          class Tag4 , class Tag5 , class Tag6 >
struct DimensionAppend< Dimension<Tag1,Tag2,Tag3,Tag4,Tag5>, Tag6>
{ typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6> type ; };

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct DimensionAppend< Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6> , Tag7 >
{
  typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> type ;
};

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct DimensionAppend< Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> , Tag8 >
{
  typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> type ;
};

//----------------------------------------------------------------------

template< class Tag1 >
struct DimensionTruncate< Dimension<Tag1> >
{ typedef Dimension<> type ; };

template< class Tag1 , class Tag2 >
struct DimensionTruncate< Dimension<Tag1,Tag2> >
{ typedef Dimension<Tag1> type ; };

template< class Tag1 , class Tag2 , class Tag3 >
struct DimensionTruncate< Dimension<Tag1,Tag2,Tag3> >
{ typedef Dimension<Tag1,Tag2> type ; };

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 >
struct DimensionTruncate< Dimension<Tag1,Tag2,Tag3,Tag4> >
{ typedef Dimension<Tag1,Tag2,Tag3> type ; };

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 , class Tag5 >
struct DimensionTruncate< Dimension<Tag1,Tag2,Tag3,Tag4,Tag5> >
{ typedef Dimension<Tag1,Tag2,Tag3,Tag4> type ; };

template< class Tag1 , class Tag2 , class Tag3 ,
          class Tag4 , class Tag5 , class Tag6 >
struct DimensionTruncate< Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6> >
{ typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5> type ; };

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 >
struct DimensionTruncate< Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> >
{ typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6> type ; };

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
struct DimensionTruncate< Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> >
{ typedef Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7> type ; };

//----------------------------------------------------------------------

void print( std::ostream & ,
            const DimensionTag * ,
            const DimensionTag * ,
            const DimensionTag * ,
            const DimensionTag * ,
            const DimensionTag * ,
            const DimensionTag * ,
            const DimensionTag * ,
            const DimensionTag * ,
            const unsigned * );

template< class Tag1 , class Tag2 , class Tag3 , class Tag4 ,
          class Tag5 , class Tag6 , class Tag7 , class Tag8 >
inline
std::ostream & operator <<
  ( std::ostream & s ,
    const Dimension<Tag1,Tag2,Tag3,Tag4,Tag5,Tag6,Tag7,Tag8> & d )
{
  print( s , Tag1::descriptor() ,
             Tag2::Descriptor() ,
             Tag3::Descriptor() ,
             Tag4::Descriptor() ,
             Tag5::Descriptor() ,
             Tag6::Descriptor() ,
             Tag7::Descriptor() ,
             Tag8::Descriptor() ,
             d.stride() );
  return s ;
}

}

#endif


