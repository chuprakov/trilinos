#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Cthulhu
#include <Cthulhu_Parameters.hpp>
#include <Cthulhu_Map.hpp>
#include <Cthulhu_MapFactory.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu.hpp>

// Gallery
#define CTHULHU_ENABLED // == Gallery have to be build with the support of Cthulhu matrices.
#include <MueLu_GalleryParameters.hpp>
#include <MueLu_MatrixFactory.hpp>

// Aggregation
#include "MueLu_UCAggregationFactory.hpp"

#include "MueLu_UseShortNames.hpp"

int main(int argc, char *argv[]) {
  
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor cmdp(false);
  
  MueLu::Gallery::Parameters matrixParameters(cmdp);   // manage parameters of the test case
  Cthulhu::Parameters cthulhuParameters(cmdp);  // manage parameters of cthulhu
  
  switch (cmdp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  cthulhuParameters.check();

  matrixParameters.print();
  cthulhuParameters.print();

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  const RCP<const Map> map = MapFactory::Build(cthulhuParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList()); //TODO: Operator vs. CrsOperator
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  int currentPrintLevel=10;
  int printFlag=6;
  
  MueLu::AggregationOptions aggOptions;
  
  aggOptions.SetPrintFlag(printFlag);      
  aggOptions.SetMinNodesPerAggregate(2);  
  aggOptions.SetMaxNeighAlreadySelected(5);
  // aggOptions.SetOrdering(1); //TODO: RandomReorder()
  aggOptions.SetOrdering(2);
  aggOptions.SetPhase3AggCreation(0.5);

  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  
  if (comm->getRank() == 0 && printFlag < currentPrintLevel)
    printf("main() Aggregate_CoarsenUncoupled : \n");
 
  RCP<Graph> graph = rcp(new Graph(Op->getCrsGraph(), "Uncoupled"));
  
  RCP<UCAggregationFactory> AggFact = rcp(new UCAggregationFactory());
  RCP<Aggregates> aggregates = AggFact->Build(*graph, aggOptions);
  
  /**********************************************************************************/
  /*                                                                                */
  /**********************************************************************************/
  
  RCP<Cthulhu::Vector<int> > Final_ = Cthulhu::VectorFactory<int>::Build( aggregates->GetVertex2AggId()->getMap() );

  {
    Teuchos::ArrayRCP<int> Final = Final_->getDataNonConst(0);
    Teuchos::ArrayRCP<const int> vertex2AggId = aggregates->GetVertex2AggId()->getData(0);
    Teuchos::ArrayRCP<const int> procWinner   = aggregates->GetProcWinner()->getData(0);

    for (size_t i = 0; i < aggregates->GetVertex2AggId()->getMap()->getNodeNumElements(); i++) 
      Final[i] = vertex2AggId[i] + procWinner[i]*1000;
  }
  
  printf("finals\n");
  cout << *Final_ << endl; sleep(2);

  return EXIT_SUCCESS;
}
