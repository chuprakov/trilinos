#
# We list MueLu as a dependency, but we only use Xpetra.  So
#   when Xpetra is split out as its own package, we will
#   depend on Xpetra only.
#
SET(LIB_REQUIRED_DEP_PACKAGES Tpetra Kokkos Epetra Teuchos MueLu)
SET(LIB_OPTIONAL_DEP_PACKAGES)
SET(TEST_REQUIRED_DEP_PACKAGES Tpetra Kokkos Epetra Teuchos EpetraExt MueLu)
SET(TEST_OPTIONAL_DEP_PACKAGES)
SET(LIB_REQUIRED_DEP_TPLS)
SET(LIB_OPTIONAL_DEP_TPLS METIS PaToH ParMETIS Scotch CCOLAMD OVIS)
SET(TEST_REQUIRED_DEP_TPLS) 
SET(TEST_OPTIONAL_DEP_TPLS METIS PaToH ParMETIS Scotch CCOLAMD OVIS)

