#include "TSFEpetraVectorSpace.h"

using namespace TSFExtended;
using namespace Teuchos;

EpetraVectorSpace::EpetraVectorSpace()
	: TSFCore::EpetraVectorSpace(),
    Handleable<TSFCore::VectorSpace<double> >(),
		Describable()
{}

EpetraVectorSpace::EpetraVectorSpace(const TSFSmartPtr<Epetra_Map>& localMap)
	: TSFCore::EpetraVectorSpace(localMap),
    Handleable<TSFCore::VectorSpace<double> >(),
		Describable()
{}

string EpetraVectorSpace::describe() const 
{
	string rtn = "EpetraVectorSpace[";
	if (ghostMap_.get()==0) 
		{
			rtn += "nLocal=" 
				+ TSFUtils::toString(localMap_->NumMyElements())
				+ " nGlobal=" 
				+ TSFUtils::toString(localMap_->NumGlobalElements());
		}
	else
		{
			rtn += "nLocal=" 
				+ TSFUtils::toString(localMap_->NumMyElements() )
				+ " nRemote=" 
				+ TSFUtils::toString(ghostMap_->NumMyElements() - localMap_->NumMyElements())
				+ " nGlobal=" 
				+ TSFUtils::toString(ghostMap_->NumGlobalElements());
		}
	rtn += "]";

  return rtn;
}



