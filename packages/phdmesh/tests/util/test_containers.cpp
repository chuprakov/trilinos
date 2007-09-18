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
#include <stdexcept>
#include <set>
#include <util/Parallel.hpp>
#include <util/FixedPoolAlloc.hpp>

//----------------------------------------------------------------------

enum { NBYTE = 4096 };

typedef std::set< int , std::less<int> , phdmesh::FixedPoolAllocator<NBYTE> > TestSet ;

void test_containers( phdmesh::ParallelMachine , std::istream & )
{
  try {
    phdmesh::FixedPoolBuffer<NBYTE> buffer ;
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

    std::cout << "TEST_CONTAINERS PASSED" << std::endl ;
  }
  catch( const std::exception & x ) {
    std::cout << "TEST_CONTAINERS FAILED: " << x.what() << std::endl ;
  }
  catch( ... ) {
    std::cout << "TEST_CONTAINERS FAILED: <unknown>" << std::endl ;
  }

}


