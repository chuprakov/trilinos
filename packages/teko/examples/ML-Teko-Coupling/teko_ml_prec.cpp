// Teuchos includes 
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Epetra includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_InverseLibrary.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

#include <iostream>

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

Teko::LinearOp readNS(Epetra_Comm & Comm);

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // build global communicator
   Epetra_MpiComm Comm(MPI_COMM_WORLD);

   Teko::LinearOp A = readNS(Comm);

   // Build the preconditioner
   /////////////////////////////////////////////////////////

   Teuchos::RCP<Teuchos::ParameterList> tekoPL = Teuchos::getParametersFromXmlFile("ml_teko.xml"); 
   Teuchos::RCP<Teuchos::ParameterList> invLibPL = Teuchos::rcpFromRef(tekoPL->sublist("Inverse Library"));

   // build an InverseLibrary
   RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*invLibPL);
   RCP<Teko::InverseFactory> inverse = invLib->getInverseFactory("ML-Teko");
   Teko::LinearOp invA = Teko::buildInverse(*inverse,A);

   return 0;
}

Teko::ModifiableLinearOp readOp(Epetra_Comm & Comm,const std::string & fileName)
{
   Epetra_CrsMatrix * crsMat = 0;
   Teko::ModifiableLinearOp output;

   int finfo = EpetraExt::MatrixMarketFileToCrsMatrix(fileName.c_str(), Comm, crsMat);

   if(finfo==0) {
      output = Thyra::nonconstEpetraLinearOp(Teuchos::rcp(crsMat));
   }
   else
      TEUCHOS_ASSERT(false);

   return output;
}

Teko::LinearOp readNS(Epetra_Comm & Comm)
{
  Teko::BlockedLinearOp blo = Teko::createBlockedOp();
  Teko::beginBlockFill(blo,2,2);
     blo->setNonconstBlock(0,0,readOp(Comm,"data/F.mm"));
     blo->setNonconstBlock(1,1,readOp(Comm,"data/C.mm"));
     blo->setNonconstBlock(0,1,readOp(Comm,"data/Bt.mm"));
     blo->setNonconstBlock(1,0,readOp(Comm,"data/B.mm"));
  Teko::endBlockFill(blo);

  return blo;
}
