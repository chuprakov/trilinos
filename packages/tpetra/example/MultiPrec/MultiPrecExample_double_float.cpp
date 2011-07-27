#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_HybridPlatform.hpp>
#include <TpetraExt_TypeStack.hpp>

#include <iostream>

#include "MultiPrecDriver.hpp"

/** \file MultiPrecExample_double_float.cpp
    \brief An example of a multi-precision algorithm, using a flexible preconditioned CG with recursive precision.
 */

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Comm;
  using Teuchos::ParameterList;

  // 
  // Get the communicator
  //
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  auto comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
  const int myImageID = comm->getRank();

  //
  // Get example parameters from command-line processor
  //  
  bool verbose = (myImageID==0);
  bool unfused = false;
  std::string matfile;
  std::string xmlfile;
  std::string machineFile;
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("matrix-file",&matfile,"Filename for matrix");
  cmdp.setOption("param-file", &xmlfile,"XML file for solver parameters");
  cmdp.setOption("machine-file",&machineFile,"Filename for XML machine description file.");
  cmdp.setOption("unfused","no-unfused",&unfused,"Test unfused iteration.");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // 
  // read machine file and initialize platform
  // 
  RCP<Teuchos::ParameterList> machinePL = Teuchos::parameterList();
  std::string defaultMachine(
    " <ParameterList>                                                               "
    "   <ParameterList name='%1=0'>                                                 "
    "     <Parameter name='NodeType'     type='string' value='Kokkos::SerialNode'/> "
    "   </ParameterList>                                                            "
    " </ParameterList>                                                              "
  );
  Teuchos::updateParametersFromXmlString(defaultMachine,machinePL.getRawPtr());
  if (machineFile != "") Teuchos::updateParametersFromXmlFile(machineFile,machinePL.getRawPtr());

  // 
  // create the platform object
  // 
  Tpetra::HybridPlatform platform(comm,*machinePL);

  // 
  // Define the type stack
  // 
  TPETRAEXT_TYPESTACK2(MPStack, double, float)

  //
  // instantiate a driver on the scalar stack
  //
  MultiPrecDriver<MPStack> driver;
  // hand output stream to driver
  if (verbose) driver.out = Teuchos::getFancyOStream(Teuchos::rcp(&std::cout,false));
  else         driver.out = Teuchos::getFancyOStream(Teuchos::rcp(new Teuchos::oblackholestream()));
  // hand matrix file to driver
  driver.matrixFile = matfile;
  // other params
  driver.unfusedTest = unfused;

  //
  // get the solver parameters
  // 
  RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
  // default solver stack parameters
  std::string xmlString(
    " <ParameterList>                                                    \n"
    "   <Parameter name='tolerance' value='1e-15' type='double'/>        \n"
    "   <Parameter name='verbose' value='2' type='int'/>                 \n"
    "   <ParameterList name='child'>                                     \n"
    "     <Parameter name='tolerance' value='1e-7' type='double'/>       \n"
    "     <Parameter name='verbose' value='1' type='int'/>               \n"
    "     <Parameter name='Extract Diagonal' value='true' type='bool'/>  \n" 
    "   </ParameterList>                                                 \n"
    " </ParameterList>                                                   \n"
  );
  Teuchos::updateParametersFromXmlString(xmlString,params.getRawPtr());
  if (xmlfile != "") Teuchos::updateParametersFromXmlFile(xmlfile,params.getRawPtr());
  // hand solver parameters to driver
  driver.params = params;

  // 
  // run the driver
  // 
  platform.runUserCode(driver);

  //
  // Print result
  if (driver.testPassed) {
    *driver.out << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}

/** \example MultiPrecExample_double_float.cpp 
    Demonstrate using Tpetra::RTI and a multi-precision flexible preconditioned CG, Tpetra::TypeStack and related utilities.
  */
