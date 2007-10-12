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

#include <ctype.h>
#include <stdexcept>
#include <string>
#include <iostream>

#include <util/ValueIOS.hpp>

namespace phdmesh {

//----------------------------------------------------------------------

bool value_ios_get_token( std::istream & s , int t , bool required )
{
  static const char method_name[] = "phdmesh::value_ios_get_token" ;

  while ( s.good() && isspace( s.peek() ) ) { s.get(); }

  const bool found = s.peek() == t ;

  if ( found ) {
    s.get();
  }
  else if ( required ) {
    const char c[2] = { (char) t , 0 };
    std::string msg ;
    msg.append( method_name )
       .append( " FAILED to find '" )
       .append( c )
       .append( "'" );
    throw std::runtime_error(msg);
  }
  return found ;
}

bool value_ios_get_block_begin( std::istream & s , bool required )
{ return value_ios_get_token( s , '{' , required ); }

bool value_ios_get_block_end( std::istream & s , bool required )
{ return value_ios_get_token( s , '}' , required ); }

void value_ios_put_block_begin( std::ostream & s ) { s << "{" ; }
void value_ios_put_block_end(   std::ostream & s ) { s << "}" ; }

//----------------------------------------------------------------------

void value_ios_get_fixed_array( std::istream & s ,
                                const ValueIOS<void> & io ,
                                const size_t stride ,
                                const size_t length ,
                                void * v )
{
  unsigned char * p = reinterpret_cast<unsigned char*>(v);
  unsigned char * const e = p + stride * length ;
  value_ios_get_block_begin(s);
  for ( ; p < e ; p += stride ) { io.getp( s , p ); }
  value_ios_get_block_end(s);
}

void value_ios_put_fixed_array( std::ostream & s ,
                                unsigned indent ,
                                const ValueIOS<void> & io ,
                                const size_t stride ,
                                const size_t length ,
                                const void * v )
{
  const unsigned char * p = reinterpret_cast<const unsigned char*>(v);
  const unsigned char * const e = p + stride * length ;
  s << "{" ;
  for ( ; p < e ; p += stride ) { s << " " ; io.putp( s , indent , p ); }
  s << " }" ;
}

void value_ios_tell_fixed_array( std::ostream & s ,
                                 unsigned indent ,
                                 const ValueIOS<void> & io ,
                                 const size_t ,
                                 const size_t length ,
                                 const void * v )
{
  s << "[" << length << "]{ " ;
  io.tellp( s , indent , v );
  s << " }" ;
}

//----------------------------------------------------------------------

void ValueIOS<short>
  ::get( std::istream & s, short & v ) const
{ s >> v ; }

void ValueIOS<short int>
  ::put( std::ostream & s, unsigned, const short int & v) const
{ s << v ; }

void ValueIOS<short int>
  ::tell( std::ostream & s, unsigned, const short int & ) const
{ s << "short" ; }

void ValueIOS<short int>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<short int>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<short int>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<short int> & ValueIOS<short int>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS<unsigned short>
  ::get( std::istream & s, unsigned short & v ) const
{ s >> v ; }

void ValueIOS<unsigned short>
  ::put( std::ostream & s, unsigned, const unsigned short & v) const
{ s << v ; }

void ValueIOS<unsigned short>
  ::tell( std::ostream & s, unsigned, const unsigned short &  ) const
{ s << "unsigned short" ; }

void ValueIOS<unsigned short>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<unsigned short>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<unsigned short>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<unsigned short> & ValueIOS<unsigned short>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS<int>
  ::get( std::istream & s, int & v ) const
{ s >> v ; }

void ValueIOS<int>
  ::put( std::ostream & s, unsigned, const int & v) const
{ s << v ; }

void ValueIOS<int>
  ::tell( std::ostream & s, unsigned, const int &  ) const
{ s << "int" ; }

void ValueIOS<int>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<int>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<int>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<int> & ValueIOS<int>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS<unsigned>
  ::get( std::istream & s, unsigned & v ) const
{ s >> v ; }

void ValueIOS<unsigned>
  ::put( std::ostream & s, unsigned, const unsigned & v) const
{ s << v ; }

void ValueIOS<unsigned>
  ::tell( std::ostream & s, unsigned, const unsigned &  ) const
{ s << "unsigned" ; }

void ValueIOS<unsigned>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<unsigned>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<unsigned>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<unsigned> & ValueIOS<unsigned>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS<long>
  ::get( std::istream & s, long & v ) const
{ s >> v ; }

void ValueIOS<long>
  ::put( std::ostream & s, unsigned, const long & v) const
{ s << v ; }

void ValueIOS<long>
  ::tell( std::ostream & s, unsigned, const long &  ) const
{ s << "long" ; }

void ValueIOS<long>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<long>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<long>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<long> & ValueIOS<long>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS<unsigned long>
  ::get( std::istream & s, unsigned long & v ) const
{ s >> v ; }

void ValueIOS<unsigned long>
  ::put( std::ostream & s, unsigned, const unsigned long & v) const
{ s << v ; }

void ValueIOS<unsigned long>
  ::tell( std::ostream & s, unsigned, const unsigned long &  ) const
{ s << "unsigned long" ; }

void ValueIOS<unsigned long>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<unsigned long>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<unsigned long>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<unsigned long> & ValueIOS<unsigned long>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS<float>
  ::get( std::istream & s, float & v ) const
{ s >> v ; }

void ValueIOS<float>
  ::put( std::ostream & s, unsigned, const float & v) const
{ s << v ; }

void ValueIOS<float>
  ::tell( std::ostream & s, unsigned, const float &  ) const
{ s << "float" ; }

void ValueIOS<float>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<float>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<float>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<float> & ValueIOS<float>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS<double>
  ::get( std::istream & s, double & v ) const
{ s >> v ; }

void ValueIOS<double>
  ::put( std::ostream & s, unsigned, const double & v) const
{ s << v ; }

void ValueIOS<double>
  ::tell( std::ostream & s, unsigned, const double &  ) const
{ s << "double" ; }

void ValueIOS<double>
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS<double>
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS<double>
  ::tellp( std::ostream & s, unsigned i, const void * v) const
{ tell( s , i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS<double> & ValueIOS<double>::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

namespace {

template<typename T>
void get_array( std::istream & s , unsigned n , T * v )
{
  T * const e = v + n ;
  value_ios_get_block_begin( s , true );
  for ( ; v < e ; ++v ) { s >> *v ; }
  value_ios_get_block_end( s , true );
}

template<typename T>
void put_array( std::ostream & s , unsigned n , const T * v )
{
  const T * const e = v + n ;
  value_ios_put_block_begin( s );
  for ( ; v < e ; ++v ) { s << " " << *v ; }
  s << " " ;
  value_ios_put_block_end( s );
}

template<typename T>
void get_vector( std::istream & s, std::vector<T> & v )
{
  v.clear();
  if ( value_ios_get_block_begin( s ) ) {
    while ( ! value_ios_get_block_end( s ) ) {
      T tmp ; s >> tmp ; v.push_back( tmp );
    }
  }
}

template<typename T>
void put_vector( std::ostream & s, const std::vector<T> & v )
{
  typedef std::vector<T> VectorType ;
  typedef typename VectorType::const_iterator const_iterator ;
  value_ios_put_block_begin(s);
  for ( const_iterator j = v.begin() ; j != v.end() ; ++j ) {
    s << " " << *j ;
  }
  s << " " ;
  value_ios_put_block_end(s);
}

}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void value_ios_get_array( std::istream & s , unsigned n , short * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const short * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const short * )
{ s << "[" << n << "]{ short }" ; }

void value_ios_get_array( std::istream & s , unsigned n , unsigned short * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const unsigned short * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const unsigned short * )
{ s << "[" << n << "]{ unsigned short }" ; }

void value_ios_get_array( std::istream & s , unsigned n , int * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const int * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const int * )
{ s << "[" << n << "]{ int }" ; }

void value_ios_get_array( std::istream & s , unsigned n , unsigned * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const unsigned * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const unsigned * )
{ s << "[" << n << "]{ unsigned }" ; }

void value_ios_get_array( std::istream & s , unsigned n , long * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const long * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const long * )
{ s << "[" << n << "]{ long }" ; }

void value_ios_get_array( std::istream & s , unsigned n , unsigned long * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const unsigned long * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const unsigned long * )
{ s << "[" << n << "]{ unsigned long }" ; }

void value_ios_get_array( std::istream & s , unsigned n , float * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const float * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const float * )
{ s << "[" << n << "]{ float }" ; }

void value_ios_get_array( std::istream & s , unsigned n , double * v )
{ get_array( s , n , v ); }

void value_ios_put_array( std::ostream & s , unsigned n , const double * v )
{ put_array( s , n , v ); }

void value_ios_tell_array( std::ostream & s , unsigned n , const double * )
{ s << "[" << n << "]{ double }" ; }

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void ValueIOS< std::vector<short> >
  ::get( std::istream & s, std::vector<short> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<short> >
  ::put( std::ostream & s, unsigned, const std::vector<short> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<short> >
  ::tell( std::ostream & s, unsigned, const std::vector<short> &  ) const
{ s << "[*]{ short }" ; }

void ValueIOS< std::vector<short> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<short> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<short> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<short> > &
ValueIOS< std::vector<short> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::vector<unsigned short> >
  ::get( std::istream & s, std::vector<unsigned short> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<unsigned short> >
  ::put( std::ostream & s, unsigned, const std::vector<unsigned short> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<unsigned short> >
  ::tell( std::ostream & s, unsigned, const std::vector<unsigned short> &  ) const
{ s << "[*]{ unsigned short }" ; }

void ValueIOS< std::vector<unsigned short> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<unsigned short> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<unsigned short> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<unsigned short> > &
ValueIOS< std::vector<unsigned short> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::vector<int> >
  ::get( std::istream & s, std::vector<int> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<int> >
  ::put( std::ostream & s, unsigned, const std::vector<int> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<int> >
  ::tell( std::ostream & s, unsigned, const std::vector<int> &  ) const
{ s << "[*]{ int }" ; }

void ValueIOS< std::vector<int> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<int> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<int> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<int> > &
ValueIOS< std::vector<int> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::vector<unsigned> >
  ::get( std::istream & s, std::vector<unsigned> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<unsigned> >
  ::put( std::ostream & s, unsigned, const std::vector<unsigned> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<unsigned> >
  ::tell( std::ostream & s, unsigned, const std::vector<unsigned> &  ) const
{ s << "[*]{ unsigned }" ; }

void ValueIOS< std::vector<unsigned> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<unsigned> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<unsigned> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<unsigned> > &
ValueIOS< std::vector<unsigned> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::vector<long> >
  ::get( std::istream & s, std::vector<long> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<long> >
  ::put( std::ostream & s, unsigned, const std::vector<long> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<long> >
  ::tell( std::ostream & s, unsigned, const std::vector<long> &  ) const
{ s << "[*]{ long }" ; }

void ValueIOS< std::vector<long> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<long> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<long> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<long> > &
ValueIOS< std::vector<long> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::vector<unsigned long> >
  ::get( std::istream & s, std::vector<unsigned long> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<unsigned long> >
  ::put( std::ostream & s, unsigned, const std::vector<unsigned long> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<unsigned long> >
  ::tell( std::ostream & s, unsigned, const std::vector<unsigned long> &  ) const
{ s << "[*]{ unsigned long }" ; }

void ValueIOS< std::vector<unsigned long> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<unsigned long> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<unsigned long> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<unsigned long> > &
ValueIOS< std::vector<unsigned long> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::vector<float> >
  ::get( std::istream & s, std::vector<float> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<float> >
  ::put( std::ostream & s, unsigned, const std::vector<float> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<float> >
  ::tell( std::ostream & s, unsigned, const std::vector<float> &  ) const
{ s << "[*]{ float }" ; }

void ValueIOS< std::vector<float> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<float> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<float> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<float> > &
ValueIOS< std::vector<float> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::vector<double> >
  ::get( std::istream & s, std::vector<double> & v ) const
{ get_vector( s , v ); } 

void ValueIOS< std::vector<double> >
  ::put( std::ostream & s, unsigned, const std::vector<double> & v) const
{ put_vector( s , v ); }

void ValueIOS< std::vector<double> >
  ::tell( std::ostream & s, unsigned, const std::vector<double> &  ) const
{ s << "[*]{ double }" ; }

void ValueIOS< std::vector<double> >
  ::getp( std::istream & s, void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::vector<double> >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::vector<double> >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

const ValueIOS< std::vector<double> > &
ValueIOS< std::vector<double> >::singleton()
{ static const ValueIOS<ValueType> io ; return io ; }

//----------------------------------------------------------------------

void ValueIOS< std::string >
  ::get( std::istream & s , std::string & v ) const
{ s >> v ; }

void ValueIOS< std::string >
  ::put( std::ostream & s , unsigned, const std::string & v ) const
{ s << v ; }

void ValueIOS< std::string >
  ::tell( std::ostream & s , unsigned, const std::string & ) const
{ s << "<non-white-space-string>" ; }

void ValueIOS< std::string >
  ::getp( std::istream & s , void * v ) const
{ get( s , * reinterpret_cast<ValueType*>(v) ); }

void ValueIOS< std::string >
  ::putp( std::ostream & s, unsigned i, const void * v) const
{ put( s , i, * reinterpret_cast<const ValueType*>(v) ); }

void ValueIOS< std::string >
  ::tellp( std::ostream & s , unsigned i, const void * v) const
{ tell( s, i, * reinterpret_cast<const ValueType*>(v) ); }

//----------------------------------------------------------------------

void ValueIOS_Quoted::get( std::istream & s , std::string & v ) const
{
  value_ios_get_token( s , '"' , true );
  v.clear();
  while ( s.good() && s.peek() != '"' ) {
    char tmp[2] = { (char) s.get() , 0 };
    v.append( tmp );
  }
  if ( s.good() ) { s.get(); }
}

void ValueIOS_Quoted
  ::put( std::ostream & s, unsigned, const std::string & v ) const
{ s << "\"" << v << "\"" ; }

void ValueIOS_Quoted
  ::tell( std::ostream & s, unsigned, const std::string & ) const
{ s << "\"<non-quote-char-string>\"" ; }

//----------------------------------------------------------------------

}

