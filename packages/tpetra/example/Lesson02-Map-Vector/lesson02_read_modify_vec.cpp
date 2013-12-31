// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// @HEADER

/*!
\example lesson02_read_modify_vec.cpp
\brief Read and modify the entries of a vector (Tpetra::Vector),
  using local indices.

\ref Tpetra_Lesson02 explains this example in detail.
*/

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

void
exampleRoutine (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                std::ostream& out)
{
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::reduceAll;

  // Print out the Tpetra software version information.
  out << Tpetra::version() << endl << endl;

  // The "Scalar" type is the type of the values stored in the Tpetra::Vector.
  typedef double scalar_type;

  // The "LocalOrdinal" (LO) type is the type of "local" indices.
  typedef int local_ordinal_type;

  // The "GlobalOrdinal" (GO) type is the type of "global" indices.
  typedef long global_ordinal_type;

  // The Kokkos "Node" type describes the type of shared-memory
  // parallelism that Tpetra will use _within_ an MPI process.
  typedef Kokkos::DefaultNode::DefaultNodeType node_type;

  // Type of the Map used to describe a parallel distribution.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  // Get a pointer to the default Kokkos Node.  This creates an
  // instance of the Node if it hasn't yet been created.  We'll
  // need this when creating the Tpetra::Map objects.
  RCP<node_type> node = Kokkos::DefaultNode::getDefaultNode ();

  //////////////////////////////////////////////////////////////////////
  // Create a Tpetra Map
  //////////////////////////////////////////////////////////////////////

  // The total (global, i.e., over all MPI processes) number of
  // entries in the Map.
  //
  // For this example, we scale the global number of entries in the
  // Map with the number of MPI processes.  That way, you can run this
  // example with any number of MPI processes and every process will
  // still have a positive number of entries.
  const Tpetra::global_size_t numGlobalEntries = comm->getSize() * 5;

  // Index base of the Map.  We choose zero-based (C-style) indexing.
  const global_ordinal_type indexBase = 0;

  // Construct a Map that puts the same number of equations on each
  // MPI process.  Pass in the Kokkos Node, so that this line of code
  // will work with any Kokkos Node type.
  RCP<const map_type> contigMap =
    rcp (new map_type (numGlobalEntries, indexBase, comm,
                       Tpetra::GloballyDistributed, node));

  //////////////////////////////////////////////////////////////////////
  // Create a Tpetra Vector
  //////////////////////////////////////////////////////////////////////

  // The type of the Tpetra::Vector specialization.
  typedef Tpetra::Vector<scalar_type, local_ordinal_type,
    global_ordinal_type, node_type> vector_type;

  // Create a Vector with the Map we created above.
  // This version of the constructor will fill in the vector with zeros.
  RCP<vector_type> x = rcp (new vector_type (contigMap));

  //////////////////////////////////////////////////////////////////////
  // Fill the Vector with a single number, or with random numbers
  //////////////////////////////////////////////////////////////////////

  // Set all entries of x to 42.0.
  x->putScalar (42.0);

  // Print the norm of x.
  out << "Norm of x (all entries are 42.0): " << x->norm2 () << endl;

  // Set the entries of x to (pseudo)random numbers.  Please don't
  // consider this a good parallel pseudorandom number generator.
  x->randomize ();

  // Print the norm of x.
  out << "Norm of x (random numbers): " << x->norm2 () << endl;

  //////////////////////////////////////////////////////////////////////
  // Read the entries of the Vector
  //////////////////////////////////////////////////////////////////////

  {
    // Get a const persisting view of the entries in the Vector.
    // "Const" means that you cannot modify the entries.
    // "Persisting" means that the view persists beyond the lifetime of the Vector.
    // Even after the Vector's destructor is called, the view won't go away.
    // If the Vector's entries change during the lifetime of a persisting view,
    // the entries of the persisting view are undefined.
    //
    // An ArrayRCP acts like an array or std::vector,
    // but is reference-counted like std::shared_ptr or Teuchos::RCP.
    // You may decrement the reference count manually by assigning Teuchos::null to it.
    // We put this code in an inner scope (in an extra pair of {})
    // so that the ArrayRCP will fall out of scope before the next example,
    // which modifies the entries of the Vector.

    ArrayRCP<const scalar_type> x_data = x->get1dView ();

    // x_data.size () may be longer than the number of local rows in the Vector,
    // so be sure to ask the Vector for its dimensions, rather than the ArrayRCP.
    const size_t localLength = x->getLocalLength ();

    // Count the local number of entries less than 0.5.
    // Use local indices to access the entries of x_data.
    size_t localCount = 0;
    for (size_t k = 0; k < localLength; ++k) {
      if (x_data[k] < 0.5) {
        ++localCount;
      }
    }

    // "reduceAll" is a type-safe templated version of MPI_Allreduce.
    // "outArg" is like taking the address using &, but it makes it
    // more clear that its argument is an output argument of a function.
    size_t globalCount = 0;
    reduceAll<int, size_t> (*comm, Teuchos::REDUCE_SUM, localCount, outArg (globalCount));

    // Find the total number of entries less than 0.5,
    // over all processes in the Vector's communicator.
    out << "x has " << globalCount << " entr" << (globalCount != 1 ? "ies" : "y")
        << " less than 0.5." << endl;
  }

  //////////////////////////////////////////////////////////////////////
  // Modify the entries of the Vector
  //////////////////////////////////////////////////////////////////////

  {
    // Get a nonconst persisting view of the entries in the Vector.
    // "Nonconst" means that you may modify the entries.
    // "Persisting" means that the view persists beyond the lifetime of the Vector.
    // Even after the Vector's destructor is called, the view won't go away.
    // If you create two nonconst persisting views of the same Vector,
    // and modify the entries of one view during the lifetime of the other view,
    // the entries of the other view are undefined.
    ArrayRCP<scalar_type> x_data = x->get1dViewNonConst ();

    // Use local indices to access the entries of x_data.
    // x_data.size () may be longer than the number of local rows in the Vector,
    // so be sure to ask the Vector for its dimensions, rather than the ArrayRCP.
    const size_t localLength = x->getLocalLength ();
    for (size_t k = 0; k < localLength; ++k) {
      // Add k (the local index) to every entry of x.
      x_data[k] += scalar_type (k);
    }

    // Print the norm of x.
    out << "Norm of x (modified random numbers): " << x->norm2 () << endl;
  }
}

//
// The same main() driver routine as in the first Tpetra lesson.
//
int
main (int argc, char *argv[])
{
  using std::endl;
  using Teuchos::RCP;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  const int myRank = comm->getRank();
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // We have a communicator and an output stream.
  // Let's do something with them!
  exampleRoutine (comm, out);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}
