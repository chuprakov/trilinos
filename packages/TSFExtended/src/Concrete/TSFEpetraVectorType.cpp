#include "TSFEpetraVectorType.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_RefCountPtr.hpp"

using namespace TSFExtended;
using namespace Teuchos;

EpetraVectorType::EpetraVectorType()
#ifdef HAVE_MPI
  : TSFCore::EpetraVectorSpaceFactory(rcp(new Epetra_MpiComm(MPI_COMM_WORLD)))
#else
    : TSFCore::EpetraVectorSpaceFactory()
#endif
{;}

RefCountPtr<const TSFCore::VectorSpace<double> > 
EpetraVectorType::createVecSpc(int dimension) const
{
	RefCountPtr<Epetra_Map> map = rcp(new Epetra_Map(dimension,
                                               0, *epetra_comm()));
	return rcp(new EpetraVectorSpace(map));
}


RefCountPtr<const TSFCore::VectorSpace<double> > 
EpetraVectorType::createVecSpc(int dimension,
                               int nLocal,
                               int firstLocal) const
{
	RefCountPtr<Epetra_Map> map = rcp(new Epetra_Map(dimension, nLocal,
                                               0, *epetra_comm()));
	return rcp(new EpetraVectorSpace(map));
}

RefCountPtr<const TSFCore::VectorSpace<double> > 
EpetraVectorType::createVecSpc(int dimension,
                               int nLocal,
                               const int* localIndices) const
{
	RefCountPtr<Epetra_Map> map = rcp(new Epetra_Map(dimension, nLocal,
                                               (int*) localIndices,
                                               0, *epetra_comm()));
	return rcp(new EpetraVectorSpace(map));
}




