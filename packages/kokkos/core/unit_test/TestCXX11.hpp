/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <Kokkos_Parallel.hpp>
namespace TestCXX11 {

template<class DeviceType>
struct FunctorAddTest{
  typedef Kokkos::View<double**,DeviceType> view_type;
  view_type a_, b_;
  typedef DeviceType device_type;
  FunctorAddTest(view_type & a, view_type &b):a_(a),b_(b) {}
  void operator() (const int& i) const {
    b_(i,0) = a_(i,1) + a_(i,2);
    b_(i,1) = a_(i,0) - a_(i,3);
    b_(i,2) = a_(i,4) + a_(i,0);
    b_(i,3) = a_(i,2) - a_(i,1);
    b_(i,4) = a_(i,3) + a_(i,4);
  }
};

template<class DeviceType>
double AddTestFunctor() {
  Kokkos::View<double**,DeviceType> a("A",100,5);
  Kokkos::View<double**,DeviceType> b("B",100,5);
  Kokkos::View<double**,typename DeviceType::host_mirror_device_type>
     h_a = Kokkos::create_mirror_view(a);
  Kokkos::View<double**,typename DeviceType::host_mirror_device_type>
     h_b = Kokkos::create_mirror_view(b);

  for(int i=0;i<100;i++) {
    for(int j=0;j<5;j++)
       h_a(i,j) = 0.1*i/(1.1*j+1.0) + 0.5*j;
  }
  Kokkos::deep_copy(a,h_a);

  Kokkos::parallel_for(100,FunctorAddTest<DeviceType>(a,b));
  Kokkos::deep_copy(h_b,b);

  double result = 0;
  for(int i=0;i<100;i++) {
      for(int j=0;j<5;j++)
         result += h_b(i,j);
    }

  return result;
}

template<class DeviceType>
double TestVariantFunctor(int test) {
  switch (test) {
    case 1: return AddTestFunctor<DeviceType>();
  }
  return 0;
}

#if defined (KOKKOS_HAVE_CXX11)
template<class DeviceType>
double AddTestLambda() {
  Kokkos::View<double**,DeviceType> a("A",100,5);
  Kokkos::View<double**,DeviceType> b("B",100,5);
  Kokkos::View<double**,typename DeviceType::host_mirror_device_type>
     h_a = Kokkos::create_mirror_view(a);
  Kokkos::View<double**,typename DeviceType::host_mirror_device_type>
     h_b = Kokkos::create_mirror_view(b);

  for(int i=0;i<100;i++) {
    for(int j=0;j<5;j++)
       h_a(i,j) = 0.1*i/(1.1*j+1.0) + 0.5*j;
  }
  Kokkos::deep_copy(a,h_a);

  Kokkos::parallel_for(100,[=](const int& i)  {
    b(i,0) = a(i,1) + a(i,2);
    b(i,1) = a(i,0) - a(i,3);
    b(i,2) = a(i,4) + a(i,0);
    b(i,3) = a(i,2) - a(i,1);
    b(i,4) = a(i,3) + a(i,4);
  });
  Kokkos::deep_copy(h_b,b);

  double result = 0;
  for(int i=0;i<100;i++) {
      for(int j=0;j<5;j++)
         result += h_b(i,j);
    }

  return result;
}
#else
template<class DeviceType>
double AddTestLambda() {
  return AddTestFunctor<DeviceType>();
}
#endif


template<class DeviceType>
double TestVariantLambda(int test) {
  switch (test) {
    case 1: return AddTestLambda<DeviceType>();
  }
  return 0;
}

template<class DeviceType>
bool Test(int test) {
  double res_functor = TestVariantFunctor<DeviceType>(test);
  double res_lambda = TestVariantLambda<DeviceType>(test);

  bool passed = true;

  if ( res_functor != res_lambda ) {
    passed = false;

    std::cout << "CXX11 ( test = "
              << test << " FAILED : "
              << res_functor << " != " << res_lambda
              << std::endl ;
  }

  return passed ;
}

}
