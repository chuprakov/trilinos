#include "TSFConfig.h"

#if HAVE_PETRA

#include "PetraVectorType.h"

#include "TSFError.h"
#include "PetraVectorSpace.h"
#include "TSFMatrixOperator.h"
#include "PetraMatrix.h"
#include "ILUKPreconditionerFactory.h"
#include "BICGSTABSolver.h"

#if HAVE_PETRA_MPI
#include <mpi.h>
#endif

#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"

using namespace TSF;


Epetra_Comm& PetraVectorType::comm()
{
#if HAVE_PETRA_MPI
	static TSFSmartPtr<Epetra_Comm> rtn = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
	static TSFSmartPtr<Epetra_Comm> rtn = new Epetra_SerialComm();
#endif
	return *rtn;
}

TSFVectorSpace PetraVectorType::createSpace(int dimension) const
{
	TSFSmartPtr<Epetra_Map> map = new Epetra_Map(dimension, 0, comm());
	return new PetraVectorSpace(map);
}

TSFVectorSpace PetraVectorType::createSpace(int dimension,
																										int nLocal,
																										int firstLocal) const
{
	TSFSmartPtr<Epetra_Map> map = new Epetra_Map(dimension, nLocal,
																						 0, comm());
	return new PetraVectorSpace(map);
}

TSFVectorSpace PetraVectorType::createSpace(int dimension,
																										int nLocal,
																										const int* localIndices) const
{
	TSFSmartPtr<Epetra_Map> map = new Epetra_Map(dimension, nLocal,
																						 (int*) localIndices,
																						 0, comm());
	return new PetraVectorSpace(map);
}


TSFVectorSpace PetraVectorType::createSpace(int dimension,
																										int nLocal,
																										const int* localIndices,
																										int nGhost,
																										const int* ghostIndices) const
{
	TSFSmartPtr<Epetra_Map> localMap = new Epetra_Map(dimension, nLocal,
																									(int*) localIndices,
																									0, comm());	

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
	comm().SumAll(&nGhost, &totalGhosts, 1);
	
	
	TSFSmartPtr<Epetra_Map> ghostMap = new Epetra_Map(dimension + totalGhosts, 
																									nLocal+nGhost,
																									(int*) allIndices,
																									0, comm());
	delete [] allIndices;
	TSFSmartPtr<Epetra_Import> importer = new Epetra_Import(*ghostMap, *localMap);
	return new PetraVectorSpace(localMap, ghostMap, importer);
}

TSFMatrixOperator* PetraVectorType::createMatrix(const TSFVectorSpace& domain,
																												 const TSFVectorSpace& range) const
{
	return new PetraMatrix(domain, range);
}

TSFLinearSolver PetraVectorType::defaultSolver() const
{
	TSFPreconditionerFactory precond = new ILUKPreconditionerFactory(2);
	TSFLinearSolver solver = new BICGSTABSolver(precond, 1.e-14, 2000);
	return solver;
}

#endif
