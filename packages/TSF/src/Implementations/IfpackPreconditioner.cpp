#include "IfpackPreconditioner.h"

IfpackLeftPreconditioner::
IfpackLeftPreconditioner(const TSFVectorSpace& domain,
												 const TSFVectorSpace& range,
												 const TSFSmartPtr<Ifpack_CrsRiluk>& prec,
												 const TSFSmartPtr<Ifpack_IlukGraph>& graph)
	: 
	TSFLinearOperatorBase(domain, range), 
	precond_(prec),
	precondGraph_(graph)
{}

void IfpackLeftPreconditioner::apply(const TSFVector& in, 
																		 TSFVector& out) const
{
	const Epetra_Vector& in 
		= TSFVectorAdapter<Epetra_Vector>::getConcrete(in);
	Epetra_Vector& out 
		= TSFVectorAdapter<Epetra_Vector>::getConcrete(out);

	int ierr = precond_->Solve(false, in, out);
}
													 
