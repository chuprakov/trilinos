// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestHelpers.hpp"
#include "Stokhos_UnitTestHelpers.hpp"

#include "Stokhos_Sacado_Kokkos_UQ_PCE.hpp"
#include "Stokhos_LegendreBasis.hpp"
#include "Stokhos_CompletePolynomialBasis.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

//
// Tests various View< Sacado::UQ::PCE<...>,...> operations work
// as expected
//

// Helper functions

template <typename kokkos_cijk_type, typename ordinal_type>
kokkos_cijk_type build_cijk(ordinal_type stoch_dim,
                            ordinal_type poly_ord)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Array;

  typedef typename kokkos_cijk_type::value_type value_type;
  typedef typename kokkos_cijk_type::device_type device_type;
  typedef Stokhos::OneDOrthogPolyBasis<ordinal_type,value_type> one_d_basis;
  typedef Stokhos::LegendreBasis<ordinal_type,value_type> legendre_basis;
  typedef Stokhos::CompletePolynomialBasis<ordinal_type,value_type> product_basis;
  typedef Stokhos::Sparse3Tensor<ordinal_type,value_type> Cijk;

  // Create product basis
  Array< RCP<const one_d_basis> > bases(stoch_dim);
  for (ordinal_type i=0; i<stoch_dim; i++)
    bases[i] = rcp(new legendre_basis(poly_ord, true));
  RCP<const product_basis> basis = rcp(new product_basis(bases));

  // Triple product tensor
  RCP<Cijk> cijk = basis->computeTripleProductTensor();

  // Kokkos triple product tensor
  kokkos_cijk_type kokkos_cijk =
    Stokhos::create_product_tensor<device_type>(*basis, *cijk);

  return kokkos_cijk;
}

template <typename scalar, typename ordinal>
inline
scalar generate_pce_coefficient( const ordinal nFEM,
                                    const ordinal nStoch,
                                    const ordinal iColFEM,
                                    const ordinal iStoch )
{
  const scalar X_fem = 100.0 + scalar(iColFEM) / scalar(nFEM);
  const scalar X_stoch =  1.0 + scalar(iStoch) / scalar(nStoch);
  return X_fem + X_stoch;
  //return 1.0;
}

template <typename ViewType>
bool
checkPCEView(const ViewType& v, Teuchos::FancyOStream& out) {
  typedef ViewType view_type;
  typedef typename view_type::size_type size_type;
  typedef typename view_type::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename host_array_type::value_type scalar_type;

  // Copy to host
  host_view_type h_v = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(h_v, v);
  host_array_type h_a = h_v;

  size_type num_rows, num_cols;

  // For layout left, sacado dimension becomes first dimension
  // instead of last
  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right || !view_type::is_contiguous) {
    num_rows = h_a.dimension_0();
    num_cols = h_a.dimension_1();
  }
  else {
    num_rows = h_a.dimension_1();
    num_cols = h_a.dimension_0();
  }
  bool success = true;
  if (is_right || !view_type::is_contiguous) {
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        scalar_type val = h_a(i,j);
        scalar_type val_expected =
          generate_pce_coefficient<scalar_type>(
            num_rows, num_cols, i, j);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }
  else {
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        scalar_type val = h_a(j,i);
        scalar_type val_expected =
          generate_pce_coefficient<scalar_type>(
            num_rows, num_cols, i, j);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }

  return success;
}

template <typename ViewType>
bool
checkConstantPCEView(const ViewType& v,
                     const typename ViewType::intrinsic_scalar_type& constant_val_expected,
                     Teuchos::FancyOStream& out) {
  typedef ViewType view_type;
  typedef typename view_type::size_type size_type;
  typedef typename view_type::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename host_array_type::value_type scalar_type;

  // Copy to host
  host_view_type h_v = Kokkos::create_mirror_view(v);
  Kokkos::deep_copy(h_v, v);
  host_array_type h_a = h_v;

  size_type num_rows, num_cols;

  // For layout left, sacado dimension becomes first dimension
  // instead of last
  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right || !view_type::is_contiguous) {
    num_rows = h_a.dimension_0();
    num_cols = h_a.dimension_1();
  }
  else {
    num_rows = h_a.dimension_1();
    num_cols = h_a.dimension_0();
  }

  bool success = true;
  if (is_right || !view_type::is_contiguous) {
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        scalar_type val = h_a(i,j);
        scalar_type val_expected =
          j == 0 ? constant_val_expected : scalar_type(0);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }
  else {
    for (size_type i=0; i<num_rows; ++i) {
      for (size_type j=0; j<num_cols; ++j) {
        scalar_type val = h_a(j,i);
        scalar_type val_expected =
          j == 0 ? constant_val_expected : scalar_type(0);
        TEUCHOS_TEST_EQUALITY(val, val_expected, out, success);
      }
    }
  }

  return success;
}

template <typename DataType, typename LayoutType, typename DeviceType>
struct ApplyView {
  typedef Kokkos::View<DataType,LayoutType,DeviceType> type;
};

struct NoLayout {};
template <typename DataType, typename DeviceType>
struct ApplyView<DataType,NoLayout,DeviceType> {
  typedef Kokkos::View<DataType,DeviceType> type;
};

//
// Tests
//

const int global_num_rows = 11;
const int global_num_cols = 9;

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v("view", cijk, num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;

  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right || !ViewType::is_contiguous) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_v);

  success = checkPCEView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_ConstantScalar, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v("view", cijk, num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, val );

  success = checkConstantPCEView(v, val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_ConstantPCE, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v("view", cijk, num_rows, num_cols);
  Scalar val = 1.2345;

  Kokkos::deep_copy( v, PCE(val) );

  success = checkConstantPCEView(v, val, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_HostArray, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v("view", cijk, num_rows, num_cols);
  host_array_type h_a = Kokkos::create_mirror_view(v);

  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right || !ViewType::is_contiguous) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_a);

  success = checkPCEView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeepCopy_DeviceArray, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename ViewType::array_type array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v("view", cijk, num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;
  array_type a = v;

  for (size_type i=0; i<num_rows; ++i)
    for (size_type j=0; j<num_cols; ++j)
      h_a(i,j) = generate_pce_coefficient<Scalar>(
        num_rows, num_cols, i, j);
  Kokkos::deep_copy(a, h_v);

  success = checkPCEView(v, out);
}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, Unmanaged, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v("view", cijk, num_rows, num_cols);
  host_view_type h_v = Kokkos::create_mirror_view(v);
  host_array_type h_a = h_v;

  bool is_right = Kokkos::Impl::is_same< typename ViewType::array_layout,
                                         Kokkos::LayoutRight >::value;
  if (is_right || !ViewType::is_contiguous) {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(i,j) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  else {
    for (size_type i=0; i<num_rows; ++i)
      for (size_type j=0; j<num_cols; ++j)
        h_a(j,i) = generate_pce_coefficient<Scalar>(
          num_rows, num_cols, i, j);
  }
  Kokkos::deep_copy(v, h_v);

  // Create unmanaged view
  ViewType v2(Kokkos::view_without_managing, v.ptr_on_device(), cijk,
              num_rows, num_cols);

  success = checkPCEView(v2, out);
}


namespace Test {

template< class ViewType >
struct PCEAtomicFunctor {
  typedef typename ViewType::device_type device_type ;

  typedef typename ViewType::value_type pce_type ;
  typedef typename pce_type::storage_type::value_type scalar_type ;

  scalar_type m_s ;
  ViewType m_v ;

  PCEAtomicFunctor( const ViewType & v , const scalar_type & s ) :
    m_v( v ), m_s( s )
  {
    Kokkos::parallel_for( m_v.dimension_0() , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( int i ) const
  {
    pce_type v( m_s );
    atomic_assign( & m_v(i) , v );
  }
};

}

TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( Kokkos_View_PCE, DeviceAtomic, Storage, Layout )
{
  typedef typename Storage::device_type Device;
  typedef typename Storage::value_type Scalar;
  typedef Sacado::UQ::PCE<Storage> PCE;
  typedef typename ApplyView<PCE*,Layout,Device>::type ViewType;
  typedef typename ViewType::size_type size_type;
  typedef typename ViewType::HostMirror host_view_type;
  typedef typename host_view_type::array_type host_array_type;
  typedef typename ViewType::array_type array_type;
  typedef typename PCE::cijk_type Cijk;

  // Build Cijk tensor
  const int stoch_dim = 2;
  const int poly_ord = 3;
  Cijk cijk = build_cijk<Cijk>(stoch_dim, poly_ord);

  const size_type num_rows = 11;
  const size_type num_cols = cijk.dimension();
  ViewType v("view", cijk, num_rows, num_cols);
  Scalar val = 1.2345;

  (void) Test::PCEAtomicFunctor<ViewType>( v , val );

  success = checkConstantPCEView(v, val, out);
}


#define VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT( STORAGE, LAYOUT )             \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy, STORAGE, LAYOUT )                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_ConstantScalar, STORAGE, LAYOUT )         \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_ConstantPCE, STORAGE, LAYOUT )            \
  TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, Unmanaged, STORAGE, LAYOUT )

// Some tests the fail, or fail to compile

  /*
    // These don't compile as there are no deep_copy overloads between
    // static view spec and its array_type.  That could be done, but
    // would require some additional work to ensure the shapes match.
    // It is simple enough to create an array_type view, so deep copying
    // between matching views isn't much more trouble.
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                                 \
    Kokkos_View_PCE, DeepCopy_HostArray, STORAGE, LAYOUT )               \
    TEUCHOS_UNIT_TEST_TEMPLATE_2_INSTANT(                       \
    Kokkos_View_PCE, DeepCopy_DeviceArray, STORAGE, LAYOUT )
  */

#define VIEW_UQ_PCE_TESTS_STORAGE( STORAGE )                            \
  using Kokkos::LayoutLeft;                                             \
  using Kokkos::LayoutRight;                                            \
  VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT(STORAGE, NoLayout)                   \
  VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT(STORAGE, LayoutLeft)                 \
  VIEW_UQ_PCE_TESTS_STORAGE_LAYOUT(STORAGE, LayoutRight)

#define VIEW_UQ_PCE_TESTS_ORDINAL_SCALAR_DEVICE( ORDINAL, SCALAR, DEVICE ) \
  typedef Stokhos::DynamicStorage<ORDINAL,SCALAR,DEVICE> DS;            \
  VIEW_UQ_PCE_TESTS_STORAGE( DS )

#define VIEW_UQ_PCE_TESTS_DEVICE( DEVICE )                              \
  VIEW_UQ_PCE_TESTS_ORDINAL_SCALAR_DEVICE( int, double, DEVICE )
