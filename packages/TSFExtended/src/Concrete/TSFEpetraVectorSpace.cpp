#include "TSFEpetraVectorSpace.hpp"
#include "TSFEpetraVector.hpp"
#include "Teuchos_Utils.hpp"

using namespace TSFExtended;
using namespace Teuchos;


EpetraVectorSpace::EpetraVectorSpace()
	: TSFCore::EpetraVectorSpace()
{}

EpetraVectorSpace::EpetraVectorSpace(const RefCountPtr<const Epetra_Map>& localMap)
	: TSFCore::EpetraVectorSpace(localMap)
{}

RefCountPtr<TSFCore::Vector<double> > EpetraVectorSpace::createMember() const
{
  RefCountPtr<Epetra_Vector> vec 
    = rcp(new Epetra_Vector(*epetra_map(), false));
  //  RefCountPtr<const TSFCore::EpetraVectorSpace> me = rcp(this, false);
  RefCountPtr<const EpetraVectorSpace> me = rcp(this, false);
  return rcp(new EpetraVector(vec, me));
}

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


