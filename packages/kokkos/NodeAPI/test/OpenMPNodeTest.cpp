/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER
*/

#include "NodeTest.hpp"

#include "Kokkos_OpenMPNode.hpp"
#include <Teuchos_ScalarTraits.hpp>
#include <omp.h>

namespace {

  using Kokkos::OpenMPNode;
  RCP<OpenMPNode> ompNode_;

  template <>
  RCP<OpenMPNode> getNode<OpenMPNode>() {
    return ompNode_;
  }

  template <>
  void initNode<OpenMPNode>() {
    Teuchos::ParameterList plist;

    int useNumThreads = 1;
    if (NodeTest::numThreads != -1) {
      useNumThreads = NodeTest::numThreads;
    }
    plist.set<int>("Num Threads",useNumThreads);
    plist.set<int>("Verbose",NodeTest::verbose);
    ompNode_ = rcp(new OpenMPNode(plist));
  }

  template <>
  std::pair<double,double> nativeTimings<float,OpenMPNode>(int N, int numIters, float &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("float,OpenMPNode init"), sTime("float,OpenMPNode sum");
    Teuchos::ArrayRCP<float> buff = Teuchos::arcp<float>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        float *bptr = buff.getRawPtr();
#pragma omp parallel for default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          bptr[i] = 1.0f;
        }
      }
    }
    float sum = 0.0f;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        const float *bptr = buff.getRawPtr();
        sum = 0.0f;
#pragma omp parallel for reduction (+:sum) default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          sum += bptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  template <>
  std::pair<double,double> nativeTimings<int,OpenMPNode>(int N, int numIters, int &result) {
    std::pair<double,double> ret;
    Teuchos::Time iTime("int,OpenMPNode init"), sTime("int,OpenMPNode sum");
    Teuchos::ArrayRCP<int> buff = Teuchos::arcp<int>(N);
    {
      Teuchos::TimeMonitor localTimer(iTime);
      for (int t=0; t < numIters; ++t) {
        int *bptr = buff.getRawPtr();
#pragma omp parallel for default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          bptr[i] = 1.0f;
        }
      }
    }
    int sum = 0;
    {
      Teuchos::TimeMonitor localTimer(sTime);
      for (int t=0; t < numIters; ++t) {
        const int *bptr = buff.getRawPtr();
        sum = 0;
#pragma omp parallel for reduction (+:sum) default(none) shared(bptr,N)
        for (int i=0; i < N; ++i) {
          sum += bptr[i];
        }
      }
    }
    result = sum;
    ret.first  = iTime.totalElapsedTime();
    ret.second = sTime.totalElapsedTime(); 
    return ret;
  }

  TEST_NODE(OpenMPNode)

}
