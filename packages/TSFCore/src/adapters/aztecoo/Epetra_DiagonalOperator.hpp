// ///////////////////////////////////////////////
// Epetra_DiagonalOperator.hpp

#ifndef EPETRA_DIAGONAL_OPERATOR_HPP
#define EPETRA_DIAGONAL_OPERATOR_HPP

#include "Epetra_Operator.h"
#include "Teuchos_RefCountPtr.hpp"

class Epetra_Vector;

namespace Epetra {

///
/** Implements <tt>Epetra_Operator</tt> for a digonal matrix.
 *
 */
class DiagonalOperator : public Epetra_Operator {
public:

	/** @name Constructors / initializers / accessors */
	//@{

	/// Construct to uninitialized
	DiagonalOperator();

	/// Calls <tt>initialize()</tt>.
	DiagonalOperator(
		const Teuchos::RefCountPtr<const Epetra_Map>      &map
		,const Teuchos::RefCountPtr<const Epetra_Vector>  &diag
		);

	///
	void initialize(
		const Teuchos::RefCountPtr<const Epetra_Map>      &map
		,const Teuchos::RefCountPtr<const Epetra_Vector>  &diag
		);

	///
 	void uninitialize(
		Teuchos::RefCountPtr<const Epetra_Map>      *map   = NULL
		,Teuchos::RefCountPtr<const Epetra_Vector>  *diag  = NULL
		);

	//@}

	/** @name Overridden from Epetra_Operator */
	//@{
	
	///
	int SetUseTranspose(bool UseTranspose);
	///
	int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
	///
	int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
	///
	double NormInf() const;
	///
	char * Label() const;
	///
	bool UseTranspose() const;
	///
	bool HasNormInf() const;
	///
	const Epetra_Comm & Comm() const;
	///
	const Epetra_Map & OperatorDomainMap() const;
	///
	const Epetra_Map & OperatorRangeMap() const;

	//@}

private:

	// /////////////////////
	// Private data members

	bool                                        UseTranspose_;
	Teuchos::RefCountPtr<const Epetra_Map>      map_;
	Teuchos::RefCountPtr<const Epetra_Vector>   diag_;

}; // class DiagonalOperator

} // namespace Epetra

#endif // EPETRA_DIAGONAL_OPERATOR_HPP
