#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_topology/topology.hpp>

#include <vector>
#include <string>
#include <sstream>
#include <algorithm>


STKUNIT_UNIT_TEST( stk_topology, lexicographical_smallest_permutation)
{
  using stk::topology;

  topology t = topology::TRI_3;

  const char nodes[]="bac";

  char permutation_nodes[4] = "";

  unsigned permutation_index = t.lexicographical_smallest_permutation(nodes);
  t.permutation_nodes(nodes,permutation_index,permutation_nodes);

  STKUNIT_EXPECT_EQ( std::string("abc"), std::string(permutation_nodes));

  permutation_index = t.lexicographical_smallest_permutation(nodes,true);
  t.permutation_nodes(nodes,permutation_index,permutation_nodes);

  STKUNIT_EXPECT_EQ( std::string("acb"), std::string(permutation_nodes));
}

STKUNIT_UNIT_TEST( stk_topology, side_node_ordinals)
{
  using stk::topology;

  const char nodes[] = "12345678";

  {
    topology t = topology::QUAD_4_2D;
    std::cout << "QUAD_4_2D side_nodes\n";
    for (unsigned s=0; s<t.num_sides(); ++s) {
      char side_nodes[9] = {};
      t.side_nodes( nodes, s, side_nodes );
      std::cout << "  " << side_nodes << std::endl;
    }
  }

  {
    topology t = topology::HEX_8;
    std::cout << "HEX_8 side_nodes\n";
    for (unsigned s=0; s<t.num_sides(); ++s) {
      char side_nodes[9] = {};
      t.side_nodes( nodes, s, side_nodes );
      std::cout << "  " << side_nodes << std::endl;
    }
  }

};


STKUNIT_UNIT_TEST( stk_topology, arbitrary_topology )
{
  using stk::topology;

  topology t = stk::create_superelement_topology(6);

  STKUNIT_EXPECT_EQ( t.num_nodes(), 6u);
  STKUNIT_EXPECT_EQ( t.rank(), topology::ELEMENT_RANK);

  STKUNIT_EXPECT_EQ( true, t.is_superelement());
  {
    std::ostringstream name;
    name << t ;
    std::string goldName("SUPERELEMENT_TOPOLOGY_6");
    STKUNIT_EXPECT_EQ( goldName, name.str() );
  }

  topology notSuper = topology::HEX_8;
  STKUNIT_EXPECT_FALSE( notSuper.is_superelement());

  topology newT = stk::create_superelement_topology(8);

  STKUNIT_EXPECT_NE( newT.num_nodes(), 6u);
  STKUNIT_EXPECT_EQ( newT.rank(), topology::ELEMENT_RANK);

  STKUNIT_EXPECT_EQ( true, newT.is_superelement());
  {
    std::ostringstream name;
    name << newT ;
    std::string goldName("SUPERELEMENT_TOPOLOGY_6");
    STKUNIT_EXPECT_NE( goldName, name.str() );
  }

  topology anotherT = stk::create_superelement_topology(6);
  STKUNIT_EXPECT_EQ(t, anotherT);

  topology badT = stk::create_superelement_topology(-2);

  STKUNIT_EXPECT_EQ( badT.rank(), topology::INVALID_RANK);
}
