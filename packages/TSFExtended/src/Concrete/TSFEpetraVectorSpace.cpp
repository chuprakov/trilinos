#include "TSFEpetraVectorSpace.hpp"
#include "Teuchos_Utils.hpp"

using namespace TSFExtended;
using namespace Teuchos;


EpetraVectorSpace::EpetraVectorSpace()
	: TSFCore::EpetraVectorSpace()
{}

EpetraVectorSpace::EpetraVectorSpace(const RefCountPtr<const Epetra_Map>& localMap)
	: TSFCore::EpetraVectorSpace(localMap)
{}

string EpetraVectorSpace::describe() const 
{
	string rtn = "EpetraVectorSpace[";
  rtn += "nLocal=" 
    + Teuchos::toString(epetra_map()->NumMyElements())
    + " nGlobal=" 
    + Teuchos::toString(epetra_map()->NumGlobalElements()) 
    + "]";


  return rtn;
}


