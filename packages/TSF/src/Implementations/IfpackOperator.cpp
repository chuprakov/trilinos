#include "TSFConfig.h"
#include "TSFOut.h"

#if HAVE_PETRA

#define PETRA_BOOL_SUPPORTED

#include "IfpackOperator.h"

#include "PetraVectorSpace.h"
#include "PetraMatrix.h"
#include "TSFUtils.h"
#include "Ifpack_CrsRiluk.h"
#include "TSFLeftPreconditioner.h"
#include "TSFPreconditioner.h"

#if HAVE_PETRA_MPI
#include "mpi.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"

using namespace TSF;

IfpackOperator::IfpackOperator(const TSFVectorSpace& domain,
			       const TSFVectorSpace& range,
			       Ifpack_CrsRiluk* prec,
			       Ifpack_IlukGraph* graph)
	: TSFLinearOperatorBase(domain, range),
	  precondGraph_(graph),
	  precond_(prec)
{}

IfpackOperator::~IfpackOperator()
{
  delete precond_;
  delete precondGraph_;
}

void IfpackOperator::apply(const TSFVector& argument, 
													 TSFVector& result) const 
{

	// get petra vectors from the abstract arguments
	const Epetra_Vector& in = PetraVector::getLocalValues(argument);
	Epetra_Vector& out = PetraVector::getLocalValues(result);

	// ifpack's solve is logically const but declared non-const since 
	// internal data changes. So, do a const_cast.
	Ifpack_CrsRiluk* p = const_cast<Ifpack_CrsRiluk*>(precond_);

	int ierr;
	{
		TSFTimeMonitor t(opTimer());
		ierr = p->Solve((int) false, in, out);
	}
	if (ierr < 0)
		{
			TSFError::raise("IfpackOperator::apply detected ierr=" 
											+ TSFUtils::toString(ierr));
		}
}

void IfpackOperator::applyAdjoint(const TSFVector& argument, 
													 TSFVector& result) const 
{

	// get petra vectors from the abstract arguments
	const Epetra_Vector& in = PetraVector::getLocalValues(argument);
	Epetra_Vector& out = PetraVector::getLocalValues(result);

	// ifpack's solve is logically const but declared non-const since 
	// internal data changes. So, do a const_cast.
	Ifpack_CrsRiluk* p = const_cast<Ifpack_CrsRiluk*>(precond_);

	int ierr;
	{
		TSFTimeMonitor t(opTimer());
		ierr = p->Solve((int) true, in, out);
	}
	if (ierr < 0)
		{
			TSFError::raise("IfpackOperator::apply detected ierr=" 
											+ TSFUtils::toString(ierr));
		}
}

TSFTimer& IfpackOperator::opTimer()
{
	static TSFSmartPtr<TSFTimer> timer= TSFTimer::getNewTimer("IFPACK mvmult");
	return *timer;
}


#endif /* HAVE_PETRA */
