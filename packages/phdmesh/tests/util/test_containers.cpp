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

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <set>
#include <util/Parallel.hpp>
#include <util/FixedPoolAlloc.hpp>
#include <util/NamedValue.hpp>

using namespace phdmesh ;

//----------------------------------------------------------------------

namespace {

enum { NBYTE = 4096 };

typedef std::set< int , std::less<int> , FixedPoolAllocator<NBYTE> > TestSet ;

void test_fixed_pool_buffer()
{
  FixedPoolBuffer<NBYTE> buffer ;
  TestSet::allocator_type allocator( buffer );
  std::less<int> compare ;

  TestSet test_set_1( compare , allocator );

  std::set<int> test_set_2 ;

  for ( int i = 0 ; i < 10 ; ++i ) {
    test_set_1.insert( i );
  }
  for ( int i = 0 ; i < 10 ; ++i ) {
    test_set_2.insert( i );
  }

  std::cout << "test_set_1 = " ;
  for ( TestSet::iterator
        i = test_set_1.begin() ; test_set_1.end() != i ; ++i ) {
    std::cout << *i << " " ;
  }
  std::cout << std::endl ;
}

}

//----------------------------------------------------------------------

namespace std {

template<typename T>
ostream & operator << ( ostream & s , const vector<T> & v )
{
  const unsigned N = v.size();
  if ( v.size() ) {
    s << v[0] ;
    for ( unsigned i = 1 ; i < N ; ++i ) { s << " " << v[i] ; }
  }
  return s ;
}

template<typename T>
istream & operator >> ( istream & s , vector<T> & v )
{
  const unsigned N = v.size();
  for ( unsigned i = 0 ; i < N ; ++i ) { s >> v[i] ; }
  return s ;
}

}

namespace {

void test_named_value()
{
  NamedValue<double> pd("d");
  NamedValue<std::vector<double> > pdv("d_vec") ;
  NamedValue<int> pi("i");
  NamedValue<std::vector<int> > piv("i_vec") ;
  NamedValue<int[3]> pi3("pi3");
  NamedValue<unsigned[3]> pui3("pui3");

  NamedValue<NamedValueSet> ps( "my_values" ) ;

  pd.value = 0 ;
  pi.value = 0 ;
  pi3.value[0] = 0 ;
  pi3.value[1] = 1 ;
  pi3.value[2] = 2 ;
  pui3.value[0] = 3 ;
  pui3.value[1] = 4 ;
  pui3.value[2] = 5 ;

  pdv.value.resize(2);
  piv.value.resize(2);

  pdv.value[0] = 0 ;
  pdv.value[1] = 1 ;
  piv.value[0] = 0 ;
  piv.value[1] = 1 ;

  ps.value.insert( pd );
  ps.value.insert( pdv );
  ps.value.insert( pi );
  ps.value.insert( piv );

  ps.value.insert( pi3 );
  ps.value.insert( pui3 );

  NamedValue<> * r_pi  = ps.value.find( pi.name );
  NamedValue<> * r_pdv = ps.value.find( pdv.name );

  r_pi->put<int>() = 10 ;
  r_pdv->put<double>(0) = 10 ;
  r_pdv->put<double>(1) = 20 ;

  {
    std::ofstream ofile( "scratch_test_container" , std::ios::out );
    ps.write( ofile );
  }
  {
    pi.value = -1 ;
    r_pi->put<int>() = -1 ;
    r_pdv->put<double>(0) = -1 ;
    r_pdv->put<double>(1) = -1 ;

    std::cout << ps.name << " = " ;
    ps.write( std::cout );
    std::cout << std::endl ;

    std::ifstream ifile( "scratch_test_container" );
    ps.read( ifile );

    std::cout << ps.name << " = " ;
    ps.write( std::cout );
    std::cout << std::endl ;
  }

  try {
    NamedValue<NamedValueSet> ps_nest1( "nest1" );
    NamedValue<NamedValueSet> ps_nest2( "nest2" );
    NamedValue<NamedValueSet> ps_nest3( "nest3" );
    ps.value.insert( ps_nest1 );
    ps_nest1.value.insert( ps_nest2 );
    ps_nest2.value.insert( pi3 );

    std::string path("nest1.nest2.pi3");
    NamedValue<> * test_find = ps.value.find( path , '.' );
    if ( test_find != & pi3 ) {
      std::cout << "  FAILED to find " << path << std::endl ;
    }
    else {
      std::cout << "  SUCCESSFULLY found " << path << std::endl ;
    }

    std::cout << ps.name << " = " ;
    ps.write( std::cout );
    std::cout << std::endl ;

    ps.tell( std::cout );
    std::cout << std::endl ;

    ps_nest3.value.insert( ps );
    ps_nest2.value.insert( ps_nest3 );
  }
  catch( const std::exception & x ) {
    std::cout << "correctly caught:" << std::endl
              << "  " << x.what() << std::endl << std::endl ;
  }

  std::cout << ps.name << " = " ;
  ps.write( std::cout );
  std::cout << std::endl ;
}

}

//----------------------------------------------------------------------

void test_containers( ParallelMachine , std::istream & )
{
  try {
    test_fixed_pool_buffer();
    test_named_value();
    std::cout << "TEST_CONTAINERS PASSED" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "TEST_CONTAINERS FAILED: " << x.what() << std::endl ;
  }
  catch( ... ) {
    std::cout << "TEST_CONTAINERS FAILED: <unknown>" << std::endl ;
  }

}


