//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#include <Teuchos_UnitTestHarness.hpp>

#include "Kokkos_ConfigDefs.hpp"
#include "Kokkos_DefaultNode.hpp"
#include "Kokkos_DefaultSparseOps.hpp"
#include "Kokkos_Version.hpp"
#include "DefaultSparseOps_TestSparseOps.hpp"

namespace {
  using Kokkos::DefaultNode;
  using Teuchos::arcp;
  using Teuchos::ArrayRCP;
  using Teuchos::null;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef double scalar_type;
  typedef int ordinal_type;
  typedef Kokkos::DefaultNode::DefaultNodeType node_type;

  // Values of command-line arguments.  CommandLineProcessor only
  // accepts double for floating-point values, and int for integer
  // values.
  double tol = 1000 * Teuchos::ScalarTraits<double>::eps ();
  int N = 100;

  // Set up the extra command-line arguments.
  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP ();
    clp.addOutputSetupOptions (true);
    clp.setOption ("numRows", &N, "Number of rows (and columns) in the matrices to test.");
    clp.setOption ("tol", &tol, "Tolerance for relative error.");
  }

  //
  // UNIT TESTS
  //

  // Test sparse matrix-(multi)vector multiply and sparse triangular solve.
  TEUCHOS_UNIT_TEST( DefaultHostSparseOps, TestSparseOps )
  {
    typedef Kokkos::DefaultHostSparseOps<scalar_type, ordinal_type, node_type> sparse_ops_type;
    RCP<Node> node = Kokkos::DefaultNode::getDefaultNode();
    TestSparseOps<sparse_ops_type> tester;
    tester.testSparseOps (N, tol);
  }
}
