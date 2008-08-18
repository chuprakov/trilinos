/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2004 Sandia Corporation and Argonne National
    Laboratory.  Under the terms of Contract DE-AC04-94AL85000 
    with Sandia Corporation, the U.S. Government retains certain 
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 
    diachin2@llnl.gov, djmelan@sandia.gov, mbrewer@sandia.gov, 
    pknupp@sandia.gov, tleurent@mcs.anl.gov, tmunson@mcs.anl.gov      
   
  ***************************************************************** */
//
//    AUTHOR: Thomas Leurent <tleurent@mcs.anl.gov>
//       ORG: Argonne National Laboratory
//    E-MAIL: tleurent@mcs.anl.gov
//
// ORIG-DATE: 13-Nov-02 at 18:05:56
//  LAST-MOD:  9-Jun-04 at 14:50:51 by Thomas Leurent
//
// DESCRIPTION:
// ============
/*! \file PatchDataTest.cpp

Unit testing of various functions in the PatchData class. 

 */
// DESCRIP-END.
//



#include "Mesquite.hpp"
#include "PatchData.hpp"
#include "PatchDataInstances.hpp"
#include "UnitUtil.hpp"

#include "cppunit/extensions/HelperMacros.h"

#include <algorithm>

using namespace Mesquite;

using std::cout;
using std::cerr;
using std::endl;

class PatchDataTest : public CppUnit::TestFixture
{
private:
  CPPUNIT_TEST_SUITE(PatchDataTest);
  CPPUNIT_TEST (test_num_corners);
  CPPUNIT_TEST (test_get_element_vertex_indices);
  CPPUNIT_TEST (test_get_vertex_element_indices);
  CPPUNIT_TEST (test_get_element_vertex_coordinates);
  CPPUNIT_TEST (test_move_free_vertices_constrained);
  CPPUNIT_TEST (test_movement_function);
  CPPUNIT_TEST (test_get_adj_elems_2d);
  CPPUNIT_TEST (test_get_minmax_element_area);
  CPPUNIT_TEST (test_get_barrier_delta);
  CPPUNIT_TEST (test_sub_patch);
  CPPUNIT_TEST (test_fill);
  CPPUNIT_TEST (test_reorder);
  CPPUNIT_TEST_SUITE_END();
   
private:

   MsqVertex vtx_0_0;
   MsqVertex vtx_0_1;
   MsqVertex vtx_1_0;
   MsqVertex vtx_1_1;
   MsqVertex vtx_2_0;
   MsqVertex vtx_2_1;

   MsqMeshEntity tri1;
   MsqMeshEntity tri2;
   MsqMeshEntity quad1;   
   
   PatchData mPatch2D;

public:
  void setUp()
  {
     MsqPrintError err(cout);

     /* our 2D set up: 2 triangles and one quad are available
       1___3___5
        |\1|   |
        |0\| 2 |
       0---2---4
     */
     vtx_0_0.set(0,0,0);
     vtx_0_1.set(0,1,0);
     vtx_1_0.set(1,0,0);
     vtx_1_1.set(1,1,0);
     vtx_2_0.set(2,0,0);
     vtx_2_1.set(2,1,0);
     
     double coords[] = { 0,0,0,
                         0,1,0,
                         1,0,0,
                         1,1,0,
                         2,0,0,
                         2,1,0 };
     
     EntityTopology types[] = { TRIANGLE, TRIANGLE, QUADRILATERAL };
     
     size_t connectivity[] = { 0, 2, 1,
                               1, 2, 3,
                               3, 2, 4, 5 };
     
     size_t counts[] = { 3, 3, 4 };
     
     mPatch2D.fill( 6, coords,
                    3, types,
                    counts, connectivity,
                    0, err );
  }
  
  void tearDown()
  {
  }
  
public:
  PatchDataTest()
    {}
  
   void test_num_corners()
   {
     MsqPrintError err(cout);
     size_t n = mPatch2D.num_corners();
     CPPUNIT_ASSERT(n==10);
   }

   void test_get_element_vertex_indices()
   {

      MsqPrintError err(cout);
      
      std::vector<size_t> vtx_ind;
      std::vector<size_t> res;

      // test we get the right vertices for element 1 (tri)
      mPatch2D.get_element_vertex_indices(1, vtx_ind, err); CPPUNIT_ASSERT(!err);
      res.push_back(1); res.push_back(2); res.push_back(3);
      CPPUNIT_ASSERT( vtx_ind==res );

      // test we get the right vertices for element 2 (quad)
      vtx_ind.clear(); res.clear();
      mPatch2D.get_element_vertex_indices(2, vtx_ind, err); CPPUNIT_ASSERT(!err);
      res.push_back(3); res.push_back(2); res.push_back(4); res.push_back(5);
      CPPUNIT_ASSERT( vtx_ind==res );
   }

   void test_get_vertex_element_indices()
   {
     /*  1___3___5
         |\1|   |
         |0\| 2 |
         0---2---4   */
     MsqPrintError err(cout);
     
     std::vector<size_t> elem_ind;
     std::vector<size_t> res;
     
     mPatch2D.generate_vertex_to_element_data();
     
     // test we get the elements contiguous to vertex 3
     mPatch2D.get_vertex_element_indices(3, elem_ind,err); CPPUNIT_ASSERT(!err);
     res.push_back(2); res.push_back(1);
     CPPUNIT_ASSERT(res==elem_ind);
     
     // test we get the elements contiguous to vertex 2
     elem_ind.clear(); res.clear();
     mPatch2D.get_vertex_element_indices(2, elem_ind,err); CPPUNIT_ASSERT(!err);
     res.push_back(2); res.push_back(1); res.push_back(0);
     CPPUNIT_ASSERT(res==elem_ind);
   }

   void test_get_element_vertex_coordinates()
   {
      MsqPrintError err(cout);

      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(1, coords,err); CPPUNIT_ASSERT(!err);
      
      CPPUNIT_ASSERT( coords[0]==vtx_0_1 );
      CPPUNIT_ASSERT( coords[1]==vtx_1_0 );
      CPPUNIT_ASSERT( coords[2]==vtx_1_1 );
   }

   /* This tests the move_vertices() function as well as the
      PatchDataCoordsMemento functionality
      */
   void test_move_free_vertices_constrained()
   {
      MsqPrintError err(cout);

      // gets a memento of the patch coordinates.
      PatchDataVerticesMemento* coords_mem = mPatch2D.create_vertices_memento(err);
      CPPUNIT_ASSERT(!err);
      
      // Move the two first vertices in direction dk by step size s;
      Vector3D dk[6];
      dk[0].set(-1,-2,0);
      dk[1].set(-1, 2,0);
      double s = 0.3;
      mPatch2D.move_free_vertices_constrained(dk, 6, s, err); CPPUNIT_ASSERT(!err);

      // gets the new coordinates and  checks the vertices were displaced as expected.
      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      Vector3D new_vtx_0_0 = vtx_0_0 + s*dk[0];
      Vector3D new_vtx_0_1 = vtx_0_1 + s*dk[1];
      CPPUNIT_ASSERT(coords[0] == new_vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == new_vtx_0_1);

      // restore the PatchData to previous coords.
      mPatch2D.set_to_vertices_memento(coords_mem, err); CPPUNIT_ASSERT(!err);

      // gets the new coordinates and  checks the vertices are back to original.
      coords.clear();
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      CPPUNIT_ASSERT(coords[0] == vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == vtx_0_1);

      delete coords_mem;
   }

   void test_movement_function()
   {
      MsqPrintError err(cout);
      // gets a memento of the patch coordinates.
      PatchDataVerticesMemento* coords_mem = mPatch2D.create_vertices_memento(err);
      CPPUNIT_ASSERT(!err);
      
      // Move the two first vertices in direction dk by step size s;
      Vector3D dk[6];
      dk[0].set(0,-2,0);
      dk[1].set(-1,0,0);
      double s = 1;
      mPatch2D.move_free_vertices_constrained(dk, 6, 1, err); CPPUNIT_ASSERT(!err);
      // gets the new coordinates and  checks the vertices were displaced as expected.
      std::vector< Vector3D > coords;
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      Vector3D new_vtx_0_0 = vtx_0_0 + s*dk[0];
      Vector3D new_vtx_0_1 = vtx_0_1 + s*dk[1];
      CPPUNIT_ASSERT(coords[0] == new_vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == new_vtx_0_1);
      double m_dist=mPatch2D.get_max_vertex_movement_squared(coords_mem,err);
      CPPUNIT_ASSERT(m_dist==4.0);
      // restore the PatchData to previous coords.
      mPatch2D.set_to_vertices_memento(coords_mem, err); CPPUNIT_ASSERT(!err);
      // gets the new coordinates and  checks the vertices are back to original.
      coords.clear();
      mPatch2D.get_element_vertex_coordinates(0, coords,err);
      CPPUNIT_ASSERT(coords[0] == vtx_0_0);
      CPPUNIT_ASSERT(coords[2] == vtx_0_1);

      delete coords_mem;
   }
  
/*Tests the function PatchData::get_adjacent_entities_via_n_dim()
  which finds the elements adjacent to a given element.  If 'n'
  equals 0 the elements must share a vertex; if 'n' equals 1 the
  elements must share an edge; and if 'n' equals 2 the elements
  must share a face.*/
   void test_get_adj_elems_2d()
   {
     MsqPrintError err(cout);
     std::vector<size_t> elems_0;
       //find elements sharing an edge with oth elem (should be 1)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 0, elems_0, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_0.back() == 1);
     std::vector<size_t> elems_1;
       //find elements sharing an edge with 1st elem (should be 0 and 2)
     mPatch2D.get_adjacent_entities_via_n_dim(1, 1, elems_1, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_1.size() == 2);
     std::vector<size_t> elems_2;
       //find elements sharing an vert with 0th elem (should be 1 and 2).
     mPatch2D.get_adjacent_entities_via_n_dim(0, 0, elems_2, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_2.size() == 2);
     std::vector<size_t> elems_3;
     //find elements sharing an face with 0th elem (should be empty).
     mPatch2D.get_adjacent_entities_via_n_dim(2, 0, elems_3, err);
     CPPUNIT_ASSERT(!err);
     CPPUNIT_ASSERT(elems_3.size() == 0);
   }

  
   void test_get_minmax_element_area()
   {
     MsqPrintError err(cout);
     double min, max;
     mPatch2D.get_minmax_element_unsigned_area(min, max, err); CPPUNIT_ASSERT(!err);

     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.5, min, 0.0001 );
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 1.0, max, 0.0001 );
   }        
     
   void test_get_barrier_delta()
   {
     MsqPrintError err(cout);
     
     PatchData pd1;
     create_six_quads_patch_with_domain(pd1, err); CPPUNIT_ASSERT(!err);
     pd1.clear_computed_info();
     double b = pd1.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b, 0.00001 );
     destroy_patch_with_domain(pd1);

     PatchData pd2;
     create_six_quads_patch_inverted_with_domain(pd2, err); CPPUNIT_ASSERT(!err);
     pd2.clear_computed_info();
     b = pd2.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.003, b, 0.00001 );
     destroy_patch_with_domain(pd2);

     PatchData pd3;
     create_twelve_hex_patch(pd3, err); CPPUNIT_ASSERT(!err);
     pd3.clear_computed_info();
     b = pd3.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0, b, 0.00001 );

     PatchData pd4;
     create_twelve_hex_patch_inverted(pd4, err); CPPUNIT_ASSERT(!err);
     pd4.clear_computed_info();
     b = pd4.get_barrier_delta(err); CPPUNIT_ASSERT(!err);
//     cout << "b : " <<b<<endl;     
     CPPUNIT_ASSERT_DOUBLES_EQUAL( 0.0025, b, 0.00001 );
   } 

  void check_sub_patch( unsigned vtx, unsigned layers, PatchData& pd, PatchData& sub );
  
  void test_sub_patch( );
  
  void test_patch_contents( bool reorder );
  
  void test_fill() { test_patch_contents(false); }
  void test_reorder() { test_patch_contents(true); }
};


CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "PatchDataTest");
CPPUNIT_TEST_SUITE_NAMED_REGISTRATION(PatchDataTest, "Unit");

void PatchDataTest::check_sub_patch( unsigned vtx, unsigned layers, PatchData& pd, PatchData& sub )
{
  unsigned i, j;
  msq_std::set<size_t> seen;
  msq_std::vector<size_t> vtx_map, elem_map;

    // test vertex list consistency 
  vtx_map.resize( sub.num_nodes() );
  for (i = 0; i < sub.num_nodes(); ++i) {
      // get index in old patch for this vertex
    Mesh::VertexHandle h = sub.get_vertex_handles_array()[i];
    Mesh::VertexHandle* end = pd.get_vertex_handles_array() + pd.num_nodes();
    Mesh::VertexHandle* ptr = msq_std::find( pd.get_vertex_handles_array(), end, h );
    CPPUNIT_ASSERT(ptr != end);
    size_t idx = ptr - pd.get_vertex_handles_array();
    CPPUNIT_ASSERT( idx < pd.num_nodes() );
// put handle in map
    vtx_map[i] = idx;
      // make sure we don't have duplicates of vertices
    CPPUNIT_ASSERT( seen.insert(idx).second );
      // make sure vertices have same coords
    CPPUNIT_ASSERT_VECTORS_EQUAL( pd.vertex_by_index(idx), sub.vertex_by_index(i), 1e-12 );
  }

    // test element list consistency 
  seen.clear();
  elem_map.resize( sub.num_elements() );
  for (i = 0; i < sub.num_elements(); ++i) {
      // get index in old patch for element
    Mesh::ElementHandle h = sub.get_element_handles_array()[i];
    Mesh::ElementHandle* end = pd.get_element_handles_array() + pd.num_nodes();
    Mesh::ElementHandle* ptr = msq_std::find( pd.get_element_handles_array(), end, h );
    CPPUNIT_ASSERT(ptr != end);
    size_t idx = ptr - pd.get_element_handles_array();
    CPPUNIT_ASSERT( idx < pd.num_elements() );
// put handle in map
    elem_map[i] = idx;
      // make sure we don't have duplicate elements
    CPPUNIT_ASSERT( seen.insert(idx).second );
      // get elements
    MsqMeshEntity& elem1 = pd.element_by_index(idx);
    MsqMeshEntity& elem2 = sub.element_by_index(i);
      // compare element data
    CPPUNIT_ASSERT_EQUAL( elem1.get_element_type(), elem2.get_element_type() );
    CPPUNIT_ASSERT_EQUAL( elem1.node_count(), elem2.node_count() );
      // get connectivity for elements
    msq_std::vector<size_t> vtx1, vtx2;
    elem1.get_node_indices( vtx1 );
    elem2.get_node_indices( vtx2 );
    CPPUNIT_ASSERT_EQUAL( vtx1.size(), vtx2.size() );
      // compare connectivity
    for (j = 0; j < vtx1.size(); ++j) {
      CPPUNIT_ASSERT( vtx1[j] < pd.num_nodes() );
      CPPUNIT_ASSERT( vtx2[j] < sub.num_nodes() );
      CPPUNIT_ASSERT_EQUAL( vtx1[j], vtx_map[vtx2[j]] );
    }
  }

    // test that the subpatch has the elements adjacent to the specified 
    // vertex.

    // first get list of adjacent elements in original patch
  seen.clear();
  for (i = 0; i < pd.num_elements(); ++i) {
    msq_std::vector<size_t> vtx_list;
    pd.element_by_index(i).get_node_indices( vtx_list );
    if (msq_std::find( vtx_list.begin(), vtx_list.end(), vtx ) != vtx_list.end())
      seen.insert(i);
  }

    // if 1 layer, then should match element count
  if (layers == 1) {
    CPPUNIT_ASSERT_EQUAL( seen.size(), sub.num_elements() );
  }

    // remove from the set each element in the subpatch
  for (i = 0; i < sub.num_elements(); ++i) {
    size_t idx = elem_map[i];
    msq_std::set<size_t>::iterator it = seen.find( idx );
    if (it != seen.end()) {
      seen.erase(it);
    }
    else {
      CPPUNIT_ASSERT( layers > 1 );
    }
  }
  CPPUNIT_ASSERT( seen.empty() );
}

void PatchDataTest::test_sub_patch( )
{
  MsqPrintError err( msq_stdio::cout );
  PatchData pd, sub;
  create_twelve_hex_patch( pd, err );
  CPPUNIT_ASSERT(!err);

  for (unsigned i = 0; i < pd.num_free_vertices(); ++i) {
    unsigned layers = i % 2 ? 2 : 1;
    pd.get_subpatch( i, layers, sub, err );
    CPPUNIT_ASSERT(!err);

    check_sub_patch( i, layers, pd, sub );
  }
}



void PatchDataTest::test_patch_contents( bool reorder )
{
  const unsigned NUM_VERTEX = 15;
  const unsigned NUM_ELEMENT = 9;

    // Mesh data used to populate patch
    // Use a relatively randomized order for verices
    // so patch reordering will result in a changed
    // vertex ordering.
  double coords[3*NUM_VERTEX] = { 6, 6, 3, // 0
                                  0, 0, 0,
                                  0, 6, 3,
                                  4, 2, 2, 
                                  2, 4, 2,
                                  4, 4, 2, // 5
                                  0, 6, 3,
                                  2, 2, 1, 
                                  2, 6, 3, 
                                  4, 0, 2, 
                                  6, 3, 3, // 10
                                  0, 4, 2,
                                  2, 0, 1,
                                  6, 2, 2,
                                  0, 2, 1 }; // 14
  size_t conn[] = { 3, 5, 4, 7,
                    7, 4, 11, 14,
                    5, 0, 8,
                    11, 5, 8,
                    13, 3, 9, 6, 
                    10, 5, 3, 13,
                    12, 7, 14, 1,
                    10, 0, 5,
                    7, 12, 9, 3 };
  size_t conn_len[NUM_ELEMENT] = { 4, 4, 3, 3, 4, 4, 4, 3, 4 };
  EntityTopology types[NUM_ELEMENT] = { QUADRILATERAL,
                                        QUADRILATERAL,
                                        TRIANGLE,
                                        TRIANGLE,
                                        QUADRILATERAL,
                                        QUADRILATERAL,
                                        QUADRILATERAL,
                                        TRIANGLE,
                                        QUADRILATERAL };
    // mark vertices along X and Y axis as fixed
  bool fixed[NUM_VERTEX] = { false, true, true, false, false, false, 
                             true, false, false, true, false, true,
                             true, false, false };

    // populate patch data
  PatchData pd;
  MsqPrintError err(msq_stdio::cout);
  pd.fill( NUM_VERTEX, coords, NUM_ELEMENT, types, conn_len, conn, fixed, err );
  CPPUNIT_ASSERT(!MSQ_CHKERR(err));
  
  if (reorder)
    pd.reorder();
    
    // count free vertices
  unsigned i, j;
  size_t num_free = 0;
  for (i = 0; i < NUM_VERTEX; ++i)
    if (!fixed[i])
      ++num_free;
  CPPUNIT_ASSERT_EQUAL( num_free, pd.num_free_vertices() );
  
    // NOTE: PatchData will reorder contents either because reorder()
    //       was called or to group vertices by fixed/free status.
    //       We will assume that the handles arrays for vertices and
    //       elements contain the initial positions in the input
    //       arrays used to populate the patch data.
  
    // Test vertex handles
  msq_std::vector<bool> seen( NUM_VERTEX, false );
  for (i = 0; i < pd.num_nodes(); ++i) {
    size_t h = (size_t)(pd.get_vertex_handles_array()[i]);
    CPPUNIT_ASSERT( h < NUM_VERTEX );
    CPPUNIT_ASSERT( !seen[h] );
    seen[h] = true;
  }
  
    // Test vertex coordinates
  for (i = 0; i < pd.num_nodes(); ++i) {
    size_t h = (size_t)(pd.get_vertex_handles_array()[i]);
    MsqVertex vtx = pd.vertex_by_index(i);
    CPPUNIT_ASSERT_DOUBLES_EQUAL( vtx[0], coords[3*h  ], DBL_EPSILON );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( vtx[1], coords[3*h+1], DBL_EPSILON );
    CPPUNIT_ASSERT_DOUBLES_EQUAL( vtx[2], coords[3*h+2], DBL_EPSILON );
  }
  
    // Test vertex fixed flags
  for (i = 0; i < pd.num_nodes(); ++i) {
    size_t h = (size_t)(pd.get_vertex_handles_array()[i]);
    if (fixed[h]) {
      CPPUNIT_ASSERT( i >= pd.num_free_vertices() );
      CPPUNIT_ASSERT( !pd.vertex_by_index(i).is_free_vertex() );
    }
    else {
      CPPUNIT_ASSERT( i < pd.num_free_vertices() );
      CPPUNIT_ASSERT( pd.vertex_by_index(i).is_free_vertex() );
    }
  }
  
    // Test element handles
  seen.clear();
  seen.resize( NUM_ELEMENT, false );
  for (i = 0; i < pd.num_elements(); ++i) {
    size_t h = (size_t)(pd.get_element_handles_array()[i]);
    CPPUNIT_ASSERT( h < NUM_ELEMENT );
    CPPUNIT_ASSERT( !seen[h] );
    seen[h] = true;
  }
  
    // Test element types 
  for (i = 0; i < pd.num_elements(); ++i) {
    size_t h = (size_t)(pd.get_element_handles_array()[i]);
    CPPUNIT_ASSERT_EQUAL( types[h], pd.element_by_index(i).get_element_type() );
  }
  
    // Test element connectivity 
  for (i = 0; i < pd.num_elements(); ++i) {
    size_t h = (size_t)(pd.get_element_handles_array()[i]);
    MsqMeshEntity& elem = pd.element_by_index(i);
    CPPUNIT_ASSERT_EQUAL( conn_len[h], elem.vertex_count() );
    CPPUNIT_ASSERT_EQUAL( conn_len[h], elem.node_count() );
    
      // calculate offset in input list for element connectivity
    unsigned conn_pos = 0;
    for (j = 0; j < h; ++j)
      conn_pos += conn_len[j];

    const size_t* elem_conn = elem.get_vertex_index_array();
    for (unsigned j = 0; j < elem.vertex_count(); ++j) {
      size_t vh = (size_t)(pd.get_vertex_handles_array()[elem_conn[j]]);
      CPPUNIT_ASSERT_EQUAL( vh, conn[conn_pos] );
      ++conn_pos;
    }
  }
}

  
