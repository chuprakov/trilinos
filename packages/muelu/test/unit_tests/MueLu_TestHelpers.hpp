#ifndef MUELU_TEST_HELPERS_H
#define MUELU_TEST_HELPERS_H

#include <iostream>

// Teuchos
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"

// Cthulhu
#include "Cthulhu_ConfigDefs.hpp"
#include "Cthulhu_DefaultPlatform.hpp"
#include "Cthulhu_Parameters.hpp"
#include "Cthulhu_MapFactory.hpp"
#include "Cthulhu_CrsOperator.hpp"

// MueLu
#include "MueLu_Exceptions.hpp"

// Gallery
#include "MueLu_MatrixFactory.hpp"
#include "MueLu_MatrixTypes.hpp"

namespace MueLu {

  namespace TestHelpers {
    
    using Cthulhu::global_size_t;
    using Teuchos::RCP;

    class Parameters {

    private:
      Parameters() {} // static class
      
    public:
      
      static Cthulhu::Parameters cthulhuParameters;
      
      inline static RCP<const Teuchos::Comm<int> > getDefaultComm() {
        return Cthulhu::DefaultPlatform::getDefaultPlatform().getComm();
      }
      
      inline static Cthulhu::UnderlyingLib getLib() {
        return MueLu::TestHelpers::Parameters::cthulhuParameters.GetLib();
      }
    };

    template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
    class Factory {
#include "MueLu_UseShortNames.hpp"

    private: 
      Factory() {} // static class

    public:

      //
      // Method that creates a map containing a specified number of local elements per process.
      //
      static const RCP<const Map> BuildMap(LO numElementsPerProc) { 

        RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();
  
        const global_size_t INVALID = Teuchos::OrdinalTraits<global_size_t>::invalid();
  
        return MapFactory::Build(TestHelpers::Parameters::getLib(), INVALID, numElementsPerProc, 0, comm);
  
      } // BuildMap()

      // Create a matrix as specified by parameter list options
      static RCP<CrsOperator> BuildMatrix(Teuchos::ParameterList &matrixList) {
        RCP<const Teuchos::Comm<int> > comm = TestHelpers::Parameters::getDefaultComm();

        GO nx,ny,nz;
        nx = ny = nz = 5;
        nx = matrixList.get("nx",nx);
        ny = matrixList.get("ny",ny);
        nz = matrixList.get("nz",nz);

        std::string matrixType = matrixList.get("matrixType","Laplace1D");
        GO numGlobalElements;
        if (matrixType == "Laplace1D")
          numGlobalElements = nx;
        else if (matrixType == "Laplace2D")
          numGlobalElements = nx*ny;
        else if (matrixType == "Laplace3D")
          numGlobalElements = nx*ny*nz;
        else {
          std::string msg = matrixType + " is unsupported (in unit testing)";
          throw(MueLu::Exceptions::RuntimeError(msg));
        }

        RCP<const Map> map = MapFactory::Build(TestHelpers::Parameters::getLib(), numGlobalElements, 0, comm);

        RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>(matrixType,map,matrixList);
        return Op;
      } // BuildMatrix()

      // Create a 1D Poisson matrix with the specified number of rows
      // nx: global number of rows
      static RCP<CrsOperator> Build1DPoisson(GO nx) {
        Teuchos::ParameterList matrixList;
        matrixList.set("nx", nx);
        matrixList.set("matrixType","Laplace1D");
        RCP<CrsOperator> A = BuildMatrix(matrixList);
        return A;
      } // Build1DPoisson()
 
    }; // class Factory

  } // namespace TestHelpers

} // namespace MueLu


// Macro to skip a test when UnderlyingLib==Epetra or Tpetra 
#define MUELU_TEST_ONLY_FOR(UnderlyingLib) \
  if (MueLu::TestHelpers::Parameters::getLib() == UnderlyingLib)

// Macro to skip a test when Epetra is used with Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal) \
  if (!(MueLu::TestHelpers::Parameters::getLib() == Cthulhu::UseEpetra && (Teuchos::OrdinalTraits<LocalOrdinal>::name() != string("int") || Teuchos::OrdinalTraits<GlobalOrdinal>::name() != string("int"))))

// Macro to skip a test when Epetra is used with Scalar != double or Ordinal != int
#define MUELU_TEST_EPETRA_ONLY_FOR_DOUBLE_AND_INT(Scalar, LocalOrdinal, GlobalOrdinal) \
  if (!(MueLu::TestHelpers::Parameters::getLib() == Cthulhu::UseEpetra && Teuchos::ScalarTraits<Scalar>::name() != string("double"))) \
    MUELU_TEST_EPETRA_ONLY_FOR_INT(LocalOrdinal, GlobalOrdinal)






//

//TODO: add directly to Teuchos ?
#include "../cthulhu/test/Cthulhu_UnitTestHelpers.hpp" // declaration of TEUCHOS_UNIT_TEST_TEMPLATE_5_DECL

#endif // ifndef MUELU_TEST_HELPERS_H
