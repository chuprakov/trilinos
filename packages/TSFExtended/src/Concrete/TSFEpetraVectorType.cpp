#include "TSFEpetraVectorType.hpp"
#include "TSFEpetraVectorSpace.hpp"
#include "Epetra_Map.h"
#include "Epetra_Import.h"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "Teuchos_RefCountPtr.hpp"
#include "TSFEpetraMatrix.hpp"

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
EpetraVectorType::createSpace(int dimension,
                              int nLocal,
                              const int* localIndices) const
{
	RefCountPtr<Epetra_Map> map = rcp(new Epetra_Map(dimension, nLocal,
                                                   (int*) localIndices,
                                                   0, *epetra_comm()));
	return rcp(new EpetraVectorSpace(map));
}

LinearOperator<double>
EpetraVectorType::createMatrix(const VectorSpace<double>& domain,
                               const VectorSpace<double>& range) const
{
  RefCountPtr<const EpetraVectorSpace> pd 
    = rcp_dynamic_cast<const EpetraVectorSpace>(domain.ptr());

  RefCountPtr<const EpetraVectorSpace> pr 
    = rcp_dynamic_cast<const EpetraVectorSpace>(range.ptr());


  TEST_FOR_EXCEPTION(pd.get()==0, runtime_error, 
                     "incompatible domain space given to "
                     "EpetraVectorType::createMatrix()");

  TEST_FOR_EXCEPTION(pr.get()==0, runtime_error, 
                     "incompatible range space given to "
                     "EpetraVectorType::createMatrix()");

  RefCountPtr<TSFCore::LinearOp<double> > A = rcp(new EpetraMatrix(pd, pr));

  return A;
}




