// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_DefaultPlatform.hpp>

#include <Kokkos_ConfigDefs.hpp>
#include <Kokkos_SerialNode.hpp>

#if defined(HAVE_KOKKOS_TBB)
#  include <Kokkos_TBBNode.hpp>
#endif // defined(HAVE_KOKKOS_TBB)

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <algorithm>

using std::endl;

namespace Tpetra {
  namespace MatrixMarket {
    namespace Test {

      /// \fn getNode
      /// \brief Return an RCP to a Kokkos Node
      ///
      template<class NodeType>
      Teuchos::RCP<NodeType>
      getNode() {
	throw std::runtime_error ("This Kokkos Node type not supported (compile-time error)");
      }

      template<>
      Teuchos::RCP<Kokkos::SerialNode>
      getNode() {
	Teuchos::ParameterList defaultParams;
	return Teuchos::rcp (new Kokkos::SerialNode (defaultParams));
      }

#if defined(HAVE_KOKKOS_TBB)
      template<>
      Teuchos::RCP<Kokkos::TBBNode>
      getNode() {
	// "Num Threads" specifies the number of threads.  Defaults to an
	// automatically chosen value.
	Teuchos::ParameterList defaultParams;
	return Teuchos::rcp (new Kokkos::TBBNode (defaultParams));
      }
#endif // defined(HAVE_KOKKOS_TBB)

      /// Test Tpetra::MatrixMarket::Reader::readSparseFile()
      ///
      /// \param filename [in] Name of the Matrix Market format sparse
      ///   matrix file to read (on MPI Rank 0 only).
      /// \param pComm [in] Communicator, over whose MPI ranks to
      ///   distribute the returned Tpetra::CrsMatrix.
      /// \param echo [in] Whether or not to echo the resulting 
      ///   matrix to cout in Matrix Market format.
      /// \param tolerant [in] Whether or not to parse the file 
      ///   tolerantly.
      /// \param verbose [in] Whether to print verbose output.
      /// \param debug [in] Whether to print debugging output.
      ///
      void
      testReadSparseFile (const std::string& filename, 
			  const Teuchos::RCP<const Teuchos::Comm<int> >& pComm,
			  const bool echo,
			  const bool tolerant, 
			  const bool verbose,
			  const bool debug)
      {
	using Teuchos::RCP;
	using std::cerr;
	using std::cout;
	using std::endl;

	typedef double scalar_type;
	typedef int local_ordinal_type;
	typedef int global_ordinal_type;
	// typedef size_t global_ordinal_type;
	// #if defined(HAVE_KOKKOS_TBB)
	//       typedef Kokkos::TBBNode node_type;
	// #else
	typedef Kokkos::SerialNode node_type;
	// #endif // defined(HAVE_KOKKOS_TBB)

	typedef Teuchos::ScalarTraits<scalar_type> STS;

	// Get a Kokkos Node instance for the particular Node type.
	RCP<node_type> pNode = getNode<node_type>();
	const int myRank = Teuchos::rank (*pComm);

	if (verbose && myRank == 0)
	  cout << "About to read Matrix Market file \"" << filename 
	       << "\":" << endl;

	// Read the sparse matrix from the given Matrix Market file.
	// This routine acts like an MPI barrier.
	const bool callFillComplete = true;
	typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, 
	  global_ordinal_type, node_type> sparse_matrix_type;
	typedef Tpetra::MatrixMarket::Reader<sparse_matrix_type> reader_type;
	RCP<sparse_matrix_type> pMatrix =
	  reader_type::readSparseFile (filename, pComm, pNode, 
				       callFillComplete, tolerant, debug);
	if (! pMatrix.is_null())
	  {
	    if (verbose && myRank == 0)
	      cout << "Successfully read Matrix Market file \"" << filename 
		   << "\"." << endl;
	  }
	else 
	  {
	    if (verbose && myRank == 0)
	      cout << "Failed to read Matrix Market file \"" << filename 
		   << "\"." << endl;
	  }
	if (echo)
	  {
	    using Tpetra::MatrixMarket::Writer;
	    typedef Writer<sparse_matrix_type> writer_type;
	    writer_type::writeSparse (cout, pMatrix);
	  }
      }
    } // namespace Test
  } // namespace MatrixMarket
} // namespace Tpetra

/// \fn main
/// \brief Benchmark driver
int 
main (int argc, char *argv[]) 
{
  using Teuchos::Comm;
  using Teuchos::CommandLineProcessor;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::cout;
  using std::endl;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &cout);
  RCP<const Comm<int> > pComm = 
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  std::string filename;  // Matrix Market file to read
  bool tolerant = false; // Parse the file tolerantly?
  bool echo = false;     // Echo the read-in matrix back?
  bool verbose = false;  // Verbosity of output
  bool debug = false;    // Print debugging info?

  CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("verbose", "quiet", &verbose,
		  "Print messages and results.");
  cmdp.setOption ("debug", "nodebug", &debug,
  		  "Print debugging information.");
  cmdp.setOption ("filename", &filename,
		  "Name of the Matrix Market sparse matrix file to read.");
  cmdp.setOption ("tolerant", "strict", &tolerant, 
		  "Whether to parse the Matrix Market file tolerantly.");
  cmdp.setOption ("echo", "noecho", &echo,
		  "Whether to echo the read-in matrix back to stdout on Rank 0 "
		  "in Matrix Market format.  Symmetric storage will have been "
		  "expanded, so the result will not be identical to the input "
		  "file, though the matrix represented will be the same.");

  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult = 
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, or does not
    // explicitly say to run the benchmark, we let this "test" pass
    // trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED)
      {
	if (Teuchos::rank(*pComm) == 0)
	  cout << "End Result: TEST PASSED" << endl;
	return EXIT_SUCCESS;
      }
    TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL, 
		       std::invalid_argument, 
		       "Failed to parse command-line arguments");
  }

  // Test reading in the sparse matrix.  If no filename or an empty
  // filename is specified, we don't invoke the test and report a
  // "TEST PASSED" message.
  if (filename != "")
    {
      using Tpetra::MatrixMarket::Test::testReadSparseFile;
      testReadSparseFile (filename, pComm, echo, tolerant, verbose, debug);
    }

  // Only Rank 0 gets to write to cout.
  if (Teuchos::rank(*pComm) == 0)
    std::cout << "End Result: TEST PASSED" << endl;
  return EXIT_SUCCESS;
}



