#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#define CTHULHU_USE_EPETRA
#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_EpetraCrsMatrix.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu.hpp>

#define CTHULHU_ENABLED //TODO
#include <MueLu_MatrixFactory.hpp>
/**********************************************************************************/

#include <MueLu_AggAlgorithm.hpp>

// For the moment, this file is just a modified version of ML_Linker.hpp

int main(int argc, char *argv[]) {
  
  std::cout << "Hello World !" << std::endl;

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  LO numThreads=1;
  GO nx=4;
  GO ny=4;
  GO nz=4;
  Teuchos::CommandLineProcessor cmdp(false,true);
  std::string matrixType("Laplace1D");
  cmdp.setOption("nt",&numThreads,"number of threads.");
  cmdp.setOption("nx",&nx,"mesh points in x-direction.");
  cmdp.setOption("ny",&ny,"mesh points in y-direction.");
  cmdp.setOption("nz",&nz,"mesh points in z-direction.");
  cmdp.setOption("matrixType",&matrixType,"matrix type: Laplace1D, Laplace2D, Star2D, Laplace3D");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  std::cout << "#threads = " << numThreads << std::endl;
  std::cout << "problem size = " << nx*ny << std::endl;
  std::cout << "matrix type = " << matrixType << std::endl;

  Teuchos::ParameterList pl;
  pl.set("Num Threads",numThreads);

  GO numGlobalElements = nx*ny;
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  LO indexBase = 0;

  RCP<const Map > map;
  map = rcp( new MyMap(numGlobalElements, indexBase, comm) );

  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>(matrixType,map,matrixList); //TODO: Operator vs. CrsOperator

  RCP<const Epetra_CrsMatrix> A;

  { // Get the underlying Epetra Mtx (Wow ! It's paintful ! => I should create a function to do that)
    RCP<const CrsMatrix> tmp_CrsMtx = Op->get_CrsMatrix();
    const RCP<const Cthulhu::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Cthulhu::EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null) { std::cout << "Error !" << std::endl; return 1; }
    A = tmp_ECrsMtx->getEpetra_CrsMatrix();
  }
  
  Graph *graph;
  std::string name = "Uncoupled";
  graph = new Graph(*A, name.c_str());
  
  int printFlag=6;
  if (graph->GetComm().MyPID() == 0 && printFlag < MueLu_PrintLevel())
    printf("main() Aggregate_CoarsenUncoupled : \n");
  
  AggregationOptions aggOptions;
  
  aggOptions.SetPrintFlag(printFlag);      
  aggOptions.SetMinNodesPerAggregate(2);  
  aggOptions.SetMaxNeighAlreadySelected(5);
  aggOptions.SetOrdering(1);
  aggOptions.SetPhase3AggCreation(0.5);
  
  Aggregates *aggregates = NULL;
  
  aggregates = MueLu_Aggregate_CoarsenUncoupled(&aggOptions,graph);

  name = "UC_CleanUp";
  MueLu_AggregateLeftOvers(&aggOptions, aggregates, name.c_str(), graph);
  
  Epetra_IntVector Final( aggregates->GetVertex2AggId()->Map() );

  for (int i = 0; i < aggregates->GetVertex2AggId()->Map().NumMyElements(); i++) 
    Final[i] = (*(aggregates->GetVertex2AggId()))[i] + (*(aggregates->GetProcWinner()))[i]*1000;
  printf("finals\n");
  cout << Final << endl; sleep(2);
  
  delete aggregates; 
  delete graph;
  
  return EXIT_SUCCESS;

}
