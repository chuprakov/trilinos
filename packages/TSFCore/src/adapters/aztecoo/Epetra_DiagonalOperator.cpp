// ///////////////////////////////////////////////
// Epetra_DiagonalOperator.cpp

#include "Epetra_DiagonalOperator.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

namespace Epetra {

// Constructors / initializers / accessors

DiagonalOperator::DiagonalOperator()
{}

DiagonalOperator::DiagonalOperator(
	const Teuchos::RefCountPtr<const Epetra_Map>      &map
	,const Teuchos::RefCountPtr<const Epetra_Vector>  &diag
	)
{
	initialize(map,diag);
}

void DiagonalOperator::initialize(
	const Teuchos::RefCountPtr<const Epetra_Map>      &map
	,const Teuchos::RefCountPtr<const Epetra_Vector>  &diag
	)
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION( map.get()==NULL, std::invalid_argument, "Error!" );
	TEST_FOR_EXCEPTION( diag.get()==NULL, std::invalid_argument, "Error!" );
	// ToDo: Validate that the maps are compatible!
#endif // _DEBUG
	map_ = map;
	diag_ = diag;
	UseTranspose_ = false;
}

void DiagonalOperator::uninitialize(
	Teuchos::RefCountPtr<const Epetra_Map>      *map
	,Teuchos::RefCountPtr<const Epetra_Vector>  *diag
	)
{
	if(map) *map = map_;
	if(diag) *diag = diag_;
	map_ = Teuchos::null;
	diag_ = Teuchos::null;
}

// Overridden from Epetra_Operator

int DiagonalOperator::SetUseTranspose(bool UseTranspose)
{
	UseTranspose_ = UseTranspose; // Totally ignored!
}

int DiagonalOperator::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
	Y.Multiply( 1.0, *diag_, X, 0.0 );  // Y = Diag * X (element-wise multiplication)
	return 0;
}

int DiagonalOperator::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
	Y.ReciprocalMultiply( 1.0, X, *diag_, 0.0 );  // Y = X / Diag (element-wise division)
	return 0;
}

double DiagonalOperator::NormInf() const
{
	double nrm[1];
	diag_->NormInf(nrm);
	return nrm[0];
}

char* DiagonalOperator::Label() const
{
	return NULL;
}

bool DiagonalOperator::UseTranspose() const
{
	return UseTranspose_;
}

bool DiagonalOperator::HasNormInf() const
{
	return true;
}

const Epetra_Comm&
DiagonalOperator::Comm() const
{
	return map_->Comm();
}

const Epetra_Map&
DiagonalOperator::OperatorDomainMap() const
{
	return *map_;
}

const Epetra_Map&
DiagonalOperator::OperatorRangeMap() const
{
	return *map_;
}

} // namespace Epetra
