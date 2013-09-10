
#include <utility>
#include <iostream>

#include <KokkosCore_config.h>

#include <BoxElemPart.hpp>

#include <Kokkos_Threads.hpp>
#if defined( KOKKOS_HAVE_CUDA )
#include <Kokkos_Cuda.hpp>
#endif

namespace Kokkos {
namespace Example {
template< class > void test_fixture();
}
}

void test_box()
{
  const unsigned global_size = 2500 ;
  const unsigned global_box[3][2] = { { 0 , 1000 } , { 0 , 1100 } , { 0 , 1200 } };

  size_t global_count = 0 ;
  size_t global_max = 0 ;
  size_t global_min = Kokkos::Example::box_count( global_box );
  size_t intersect_error = 0 ;
  size_t neighbor_max = 0 ;

  for ( unsigned global_rank = 0 ; global_rank < global_size ; ++global_rank ) {
    unsigned box[3][2] = { { 0 , global_box[0][1] } , { 0 , global_box[1][1] } , { 0 , global_box[2][1] } };
    unsigned ghost_box[3][2] ;
    size_t neighbor_count = 0 ;

    Kokkos::Example::box_partition( global_size , global_rank , global_box , box );

    Kokkos::Example::box_ghost_layer( global_box , box , 1 , ghost_box );

    {
      const size_t n = Kokkos::Example::box_count( box );

      global_max = std::max( global_max , n );
      global_min = std::min( global_min , n );
      global_count += n ;
    }

    for ( unsigned other_rank = 0 ; other_rank  < global_size ; ++other_rank ) {

      if ( other_rank == global_rank ) continue ;

      unsigned other_box[3][2] = { { 0 , global_box[0][1] } , { 0 , global_box[1][1] } , { 0 , global_box[2][1] } };
      unsigned intersect_box[3][2] ;

      Kokkos::Example::box_partition( global_size , other_rank , global_box , other_box );

      Kokkos::Example::box_intersect( intersect_box , box , other_box );

      const size_t n = Kokkos::Example::box_count( intersect_box );

      intersect_error += n ;

      Kokkos::Example::box_intersect( intersect_box , ghost_box , other_box );

      neighbor_count += Kokkos::Example::box_count( intersect_box ) ? 1 : 0 ;

      if ( n ) {
        std::cout << "box = {"
                  << " [ " << box[0][0] << " , " << box[0][1] << " )"
                  << " [ " << box[1][0] << " , " << box[1][1] << " )"
                  << " [ " << box[2][0] << " , " << box[2][1] << " )"
                  << " }" << std::endl ;
        std::cout << "other_box = {"
                  << " [ " << other_box[0][0] << " , " << other_box[0][1] << " )"
                  << " [ " << other_box[1][0] << " , " << other_box[1][1] << " )"
                  << " [ " << other_box[2][0] << " , " << other_box[2][1] << " )"
                  << " }" << std::endl ;
        return ;
      }
    }

    neighbor_max = std::max( neighbor_max , neighbor_count );
  }

  std::cout << "count( global_box ) = " << Kokkos::Example::box_count( global_box ) << std::endl ;
  std::cout << "sum partition( global_box ) = " << global_count << std::endl ;
  std::cout << "avg partition( global_box ) = " << size_t( double(global_count) / double(global_size)) << std::endl ;
  std::cout << "min partition( global_box ) = " << global_min << std::endl ;
  std::cout << "max partition( global_box ) = " << global_max << std::endl ;
  std::cout << "sum intersect( global_box ) = " << intersect_error << std::endl ;
  std::cout << "max neighbor = " << neighbor_max << std::endl ;
}

void test_elem()
{
  const Kokkos::Example::BoxElemPart::Decompose
    decompose = Kokkos::Example::BoxElemPart:: DecomposeElem ; // DecomposeElem | DecomposeNode ;
  const unsigned global_size = 256 ;
  const unsigned global_nx = 100 ;
  const unsigned global_ny = 120 ;
  const unsigned global_nz = 140 ;

  double node_count_avg = 0 ;
  size_t node_count_max = 0 ;
  size_t node_count_min = ( global_nx + 1 ) * ( global_ny + 1 ) * ( global_nz + 1 );
  double elem_count_avg = 0 ;
  size_t elem_count_max = 0 ;
  size_t elem_count_min = global_nx * global_ny * global_nz ;
  double recv_count_avg = 0 ;
  size_t recv_count_max = 0 ;
  size_t recv_count_min = global_size ;
  double send_count_avg = 0 ;
  size_t send_count_max = 0 ;
  size_t send_count_min = global_size ;

  for ( unsigned r = 0 ; r < global_size ; ++r ) {
    const Kokkos::Example::BoxElemPart
       fixture( Kokkos::Example::BoxElemPart::ElemLinear ,
                decompose , global_size , r , global_nx , global_ny , global_nz );

    // Print a sample:

    // if ( r == global_size * 2 / 3 ) fixture.print( std::cout );

    // Verify recv/send alignment:

    {
      unsigned recv_lid = fixture.owns_node_count();

      for ( unsigned i = 0 ; i < fixture.recv_node_msg_count() ; ++i ) {
        const unsigned recv_rank  = fixture.recv_node_rank( i );
        const unsigned recv_count = fixture.recv_node_count( i );

        const Kokkos::Example::BoxElemPart other_fixture(
           Kokkos::Example::BoxElemPart::ElemLinear ,
           decompose , global_size , recv_rank , global_nx , global_ny , global_nz );

        unsigned send_item = 0 ;

        unsigned j = 0 ;
        while ( j < other_fixture.send_node_msg_count() && other_fixture.send_node_rank(j) != r ) {
          send_item += other_fixture.send_node_count( j );
           ++j ;
        }

        if ( recv_count != other_fixture.send_node_count(j) ) {
          std::cout << "Error P[" << r << "].recv(" << recv_count << ") != "
                    << "P[" << recv_rank << "].send(" << other_fixture.send_node_count(j) << ")"
                    << std::endl ;
        }
        else {

          for ( unsigned k = 0 ; k < recv_count ; ++k , ++send_item , ++recv_lid ) {

            const unsigned send_lid = other_fixture.send_node_id( send_item );

            unsigned recv_coord[3] , send_coord[3] ;

            fixture.local_node_coord( recv_lid , recv_coord );

            other_fixture.local_node_coord( send_lid , send_coord );

            if ( recv_coord[0] != send_coord[0] ||
                 recv_coord[1] != send_coord[1] ||
                 recv_coord[2] != send_coord[2] ) {
              std::cout << "Error P[" << r << "].recv[" << recv_lid << "]{ "
                        << recv_coord[0] << " , "
                        << recv_coord[1] << " , "
                        << recv_coord[2] << " } != "
                        << "P[" << recv_rank << "].send[" << send_lid << "]{ "
                        << send_coord[0] << " , "
                        << send_coord[1] << " , "
                        << send_coord[2] << " }"
                        << std::endl ;
            }
          }
        }
      }
    }

    node_count_avg += fixture.owns_node_count();
    elem_count_avg += fixture.uses_elem_count();
    recv_count_avg += fixture.recv_node_msg_count();
    send_count_avg += fixture.send_node_msg_count();

    elem_count_min = std::min( (size_t) fixture.uses_elem_count() , elem_count_min );
    elem_count_max = std::max( (size_t) fixture.uses_elem_count() , elem_count_max );
    node_count_min = std::min( (size_t) fixture.owns_node_count() , node_count_min );
    node_count_max = std::max( (size_t) fixture.owns_node_count() , node_count_max );

    recv_count_max = std::max( (size_t) fixture.recv_node_msg_count() , recv_count_max );
    recv_count_min = std::min( (size_t) fixture.recv_node_msg_count() , recv_count_min );
    send_count_max = std::max( (size_t) fixture.send_node_msg_count() , send_count_max );
    send_count_min = std::min( (size_t) fixture.send_node_msg_count() , send_count_min );
  }

  node_count_avg /= double(global_size);
  elem_count_avg /= double(global_size);
  recv_count_avg /= double(global_size);
  send_count_avg /= double(global_size);

  std::cout << "Elem min(" << elem_count_min << ") avg(" << elem_count_avg << ") max(" << elem_count_max << ") " << std::endl
            << "Node min(" << node_count_min << ") avg(" << node_count_avg << ") max(" << node_count_max << ") " << std::endl
            << "Recv min(" << recv_count_min << ") avg(" << recv_count_avg << ") max(" << recv_count_max << ") " << std::endl
            << "Send min(" << send_count_min << ") avg(" << send_count_avg << ") max(" << send_count_max << ") " << std::endl
            ;
}

int main()
{
//  test_box();
//  test_elem();
  {
    std::cout << "test_fixture< Threads >" << std::endl ;
    Kokkos::Threads::initialize( std::pair<unsigned,unsigned>( 1 , 1 ) );
    Kokkos::Example::test_fixture< Kokkos::Threads >();
    Kokkos::Threads::finalize();
  }
#if defined( KOKKOS_HAVE_CUDA )
  {
    std::cout << "test_fixture< Cuda >" << std::endl ;
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
    Kokkos::Example::test_fixture< Kokkos::Cuda >();
    Kokkos::Cuda::finalize();
  }
#endif
}

