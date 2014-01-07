// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>
#ifdef HAVE_TPETRA_KOKKOSCOMPAT
#include <KokkosCore_config.h>
#ifdef KOKKOS_USE_CUDA_BUILD
  #define DO_COMPILATION
#else
  #ifndef KOKKOS_HAVE_CUDA
    #define DO_COMPILATION
  #endif
#endif
#else
  #define DO_COMPILATION
#endif

#ifdef DO_COMPILATION

#include <Teuchos_UnitTestHarness.hpp>
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_ETIHelperMacros.h>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>

using Teuchos::RCP;
using Teuchos::rcp;

typedef int LO; // Local Ordinal
typedef int GO; // Global Ordinal
typedef double scalar_type; // scalar type

typedef Tpetra::DefaultPlatform::DefaultPlatformType           platform_type;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType node_type;
typedef Tpetra::Map<LO,GO,node_type>                           map_type;
typedef Tpetra::CrsMatrix<scalar_type, LO, GO, node_type>      matrix_type;

TEUCHOS_UNIT_TEST( CrsMatrix, Bug5978 )
{
  using std::endl;

  platform_type& platform (Tpetra::DefaultPlatform::getDefaultPlatform ());
  RCP<const Teuchos::Comm<int> > comm = platform.getComm ();
  RCP<node_type> node = platform.getNode ();
  const int numProc = comm->getSize ();
  const int myRank = comm->getRank ();

  if (numProc != 2) {
    out << "This test must be run with exactly 2 MPI processes, but you ran "
        << "it with " << numProc << " process" << (numProc != 1 ? "es" : "")
        << "." << endl;
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::runtime_error, "This test must be run with exactly "
      "2 MPI processes, but you ran it with " << numProc << " process"
      << (numProc != 1 ? "es" : "") << ".");
  }
  out << "Proc " << myRank << ": Running test" << endl;
  comm->barrier ();

  Teuchos::Array<GO> rowGIDs, colGIDs;
  if (myRank == 0) {
    rowGIDs.push_back (0);
    rowGIDs.push_back (1);
    colGIDs.push_back (0);
    colGIDs.push_back (1);
  }

  out << "Proc " << myRank << ": Creating row and column Maps" << endl;
  comm->barrier ();

  const GO INVALID = Teuchos::OrdinalTraits<GO>::invalid ();
  RCP<const map_type> rowMap =
    rcp (new map_type (INVALID, rowGIDs (), 0, comm, node));
  RCP<const map_type> colMap =
    rcp (new map_type (INVALID, colGIDs (), 0, comm, node));

  out << "Proc " << myRank << ": Creating matrix" << endl;
  comm->barrier ();

  Teuchos::ArrayRCP<size_t> count (rowGIDs.size ());
  count.assign (rowGIDs.size (), 1);
  matrix_type A (rowMap, colMap, count);

  out << "Proc " << myRank << ": Filling matrix" << endl;
  comm->barrier ();

  Teuchos::Array<LO> column (1);
  Teuchos::Array<scalar_type> entry (1, 7);
  for (int i = 0; i < int (rowGIDs.size ()); ++i) {
    column[0] = i;
    A.insertLocalValues (i, column (), entry ());
  }

  out << "Proc " << myRank << ": Calling fillComplete" << endl;
  comm->barrier ();

  A.fillComplete (colMap, rowMap);

  comm->barrier ();
  out << "Proc " << myRank << ": Done with test" << endl;
}


#endif  //DO_COMPILATION

