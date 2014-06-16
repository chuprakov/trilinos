/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
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

#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <cstdio>
#include <cstdlib>

typedef Kokkos::Impl::DefaultDeviceType Device;
typedef Device::host_mirror_device_type Host;
#define TS 16

struct find_2_tuples {
  int chunk_size;
  Kokkos::View<const int*> data;
  Kokkos::View<int**> histogram;

  find_2_tuples(int chunk_size_, Kokkos::DualView<int*> data_,
                Kokkos::DualView<int**> histogram_):chunk_size(chunk_size_),
                data(data_.d_view),histogram(histogram_.d_view) {
      data_.sync<Device>();
      histogram_.sync<Device>();
      histogram_.modify<Device>();
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (Device dev) const {
    Kokkos::View<int**,Kokkos::MemoryUnmanaged> l_histogram(dev,TS,TS);
    Kokkos::View<int*,Kokkos::MemoryUnmanaged> l_data(dev,chunk_size+1);

    const int i = dev.league_rank() * chunk_size;
    for(int j = dev.team_rank(); j<chunk_size+1; j+=dev.team_size())
      l_data(j) = data(i+j);

    for(int k = dev.team_rank(); k < TS; k+=dev.team_size())
      for(int l = 0; l < TS; l++)
        l_histogram(k,l) = 0;
    dev.team_barrier();

    for(int j = 0; j<chunk_size; j++) {
      for(int k = dev.team_rank(); k < TS; k+=dev.team_size())
        for(int l = 0; l < TS; l++) {
          if((l_data(j) == k) && (l_data(j+1)==l))
            l_histogram(k,l)++;
        }
    }

    for(int k = dev.team_rank(); k < TS; k+=dev.team_size())
      for(int l = 0; l < TS; l++) {
        Kokkos::atomic_fetch_add(&histogram(k,l),l_histogram(k,l));
      }
    dev.team_barrier();
  }
  size_t shmem_size() const { return sizeof(int)*(chunk_size+2 + TS*TS); }
};

int main(int narg, char* args[]) {
  Kokkos::initialize(narg,args);
  
  int chunk_size = 1024;
  int nchunks = 100000; //1024*1024;
  Kokkos::DualView<int*> data("data",nchunks*chunk_size+1);

  srand(1231093);

  for(int i = 0; i < data.dimension_0(); i++) {
    data.h_view(i) = rand()%TS;
  }
  data.modify<Host>();
  data.sync<Device>();

  Kokkos::DualView<int**> histogram("histogram",TS,TS);


  Kokkos::Impl::Timer timer;
  Kokkos::parallel_for(
      Kokkos::ParallelWorkRequest(nchunks,TS<Device::team_max()?TS:Device::team_max()),
      find_2_tuples(chunk_size,data,histogram));
  Kokkos::fence();
  double time = timer.seconds();

  histogram.sync<Host>();

  printf("Time: %lf \n\n",time);
  int sum = 0;
  for(int k=0; k<TS; k++) {
    for(int l=0; l<TS; l++) {
      printf("%i ",histogram.h_view(k,l));
      sum += histogram.h_view(k,l);
    }
    printf("\n");
  }
  printf("Result: %i %i\n",sum,chunk_size*nchunks);
  Kokkos::finalize();
}

