#include "TSFConfig.h"
#include "PetraVectorSpace.h"

#include "PetraVector.h"
#include "TSFUtils.h"


#if HAVE_PETRA

#include "Epetra_MpiComm.h"


using namespace TSF;



PetraVectorSpace::PetraVectorSpace(const TSFSmartPtr<Epetra_Map>& localMap,
																	 const TSFSmartPtr<Epetra_Map>& ghostMap,
																	 const TSFSmartPtr<Epetra_Import>& import)
	: TSFMPIVectorSpace(),
		localMap_(localMap), ghostMap_(ghostMap), importer_(import)
{
#if HAVE_MPI
	const Epetra_Comm& epetraComm = localMap->Comm();
	const Epetra_MpiComm& comm = dynamic_cast<const Epetra_MpiComm&>(epetraComm);
	setMPIComm(comm.Comm());
#endif
}

PetraVectorSpace::PetraVectorSpace(const TSFSmartPtr<Epetra_Map>& localMap)
	: TSFMPIVectorSpace(),
		localMap_(localMap), ghostMap_(0), importer_(0)
{
	
#if HAVE_MPI
	const Epetra_Comm& epetraComm = localMap->Comm();
	const Epetra_MpiComm& comm = dynamic_cast<const Epetra_MpiComm&>(epetraComm);
	setMPIComm(comm.Comm());
#endif
}

TSFVectorSpaceBase* PetraVectorSpace::deepCopy() const 
{
	return new PetraVectorSpace(*this);
}

int PetraVectorSpace::dim() const 
{
	return localMap_->NumGlobalElements();
}

TSFVectorBase* PetraVectorSpace::createMember(const TSFVectorSpace& handle) const
{
	if (ghostMap_.get()==0)
		{
			Epetra_Vector* localValues = new Epetra_Vector(*localMap_);
			return new PetraVector(handle, localValues);
		}
	else
		{
			Epetra_Vector* allValues = new Epetra_Vector(*ghostMap_);
			double* values;
			allValues->ExtractView(&values);
			Epetra_Vector* localValues 
				= new Epetra_Vector(View, *localMap_, values);
			return new PetraVector(handle, localValues, allValues);
		}
}


void PetraVectorSpace::print(ostream& os) const 
{
	string rtn = "PetraVectorSpace[";
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

	os << rtn;
}


TSFSmartPtr<Epetra_Map> PetraVectorSpace::getMap(const TSFVectorSpace& space)
{
	if (getPtr(space)->ghostMap_.get()==0) return getPtr(space)->localMap_;
	return getPtr(space)->ghostMap_;
}

TSFSmartPtr<Epetra_Map> PetraVectorSpace::getLocalMap(const TSFVectorSpace& space)
{
	return getPtr(space)->localMap_;
}

TSFSmartPtr<Epetra_Import> PetraVectorSpace::getImporter(const TSFVectorSpace& space)
{
	return getPtr(space)->importer_;
}

const PetraVectorSpace* PetraVectorSpace::getPtr(const TSFVectorSpace& space)
{
	const PetraVectorSpace* ptr 
		= dynamic_cast<const PetraVectorSpace*>(space.ptr());
	if (ptr==0) TSFError::raise("bad cast in PetraVectorSpace::getPtr");
	return ptr;
}



#endif
