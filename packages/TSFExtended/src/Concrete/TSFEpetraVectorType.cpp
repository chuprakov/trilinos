#include "EpetraVectorType.h"

using namespace TSFExtended;
using namespace Teuchos;


RefCountPtr<TSFCore::VectorSpace> 
EpetraVectorType::createVecSpc(int dimension,
                               int nLocal,
                               int firstLocal) const
{
	RefCountPtr<Epetra_Map> map = rcp(new Epetra_Map(dimension, nLocal,
                                               0, *epetra_comm()));
	return rcp(new EpetraVectorSpace(map));
}

RefCountPtr<TSFCore::VectorSpace> 
EpetraVectorType::createVecSpc(int dimension,
                               int nLocal,
                               const int* localIndices) const
{
	RefCountPtr<Epetra_Map> map = rcp(new Epetra_Map(dimension, nLocal,
                                               (int*) localIndices,
                                               0, *epetra_comm()));
	return rcp(new EpetraVectorSpace(map));
}


RefCountPtr<TSFCore::VectorSpace> 
EpetraVectorType::createVecSpc(int dimension,
                               int nLocal,
                               const int* localIndices,
                               int nGhost,
                               const int* ghostIndices) const
{
	RefCountPtr<Epetra_Map> localMap = rcp(new Epetra_Map(dimension, nLocal,
                                                    (int*) localIndices,
                                                    0, *epetra_comm()));	

	int* allIndices = new int[nLocal + nGhost];
	for (int i=0; i<nLocal; i++)
		{
			allIndices[i] = localIndices[i];
		}
	for (int i=0; i<nGhost; i++)
		{
			allIndices[i+nLocal] = ghostIndices[i];
		}

	int totalGhosts;
	epetra_comm()->SumAll(&nGhost, &totalGhosts, 1);
	
	
	RefCountPtr<Epetra_Map> ghostMap = rcp(new Epetra_Map(dimension + totalGhosts, 
                                                        nLocal+nGhost,
                                                        (int*) allIndices,
                                                        0, *epetra_comm()));
	delete [] allIndices;
	RefCountPtr<Epetra_Import> importer = rcp(new Epetra_Import(*ghostMap, *localMap));
	return rcp(new EpetraVectorSpace(localMap, ghostMap, importer));
}


