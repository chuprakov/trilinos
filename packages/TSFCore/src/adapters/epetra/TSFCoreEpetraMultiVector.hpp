// /////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraMultiVector.hpp

#ifndef TSFCORE_EPETRA_MULTI_VECTOR_HPP
#define TSFCORE_EPETRA_MULTI_VECTOR_HPP

#include <vector>

#include "TSFCoreMPIMultiVectorBase.hpp"

///
class Epetra_MultiVector;

namespace TSFCore {

///
class EpetraVectorSpace;

///
/** Optimized <tt>MultiVector</tt> subclass for <tt>Epetra_MultiVector</tt>.
 *
 * ToDo: Finish documentation!
 */
class EpetraMultiVector : virtual public MPIMultiVectorBase<RTOp_value_type> {
public:
	
	///
	typedef RTOp_value_type Scalar;

	/** @name Constructors/Initializers */
	//@{

	///
	/** Construct to initialized.
	 *
	 * Postconditions:<ul>
	 * <tt> <tt>this->range().get() == NULL</tt>
	 * <tt> <tt>this->domain().get() == NULL</tt>
	 * </ul>
	 */
	EpetraMultiVector();

	/// Calls <tt>initalize()</tt>.
	EpetraMultiVector(
		const Teuchos::RefCountPtr<Epetra_MultiVector>         &epetra_multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_range       = Teuchos::null
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_domain      = Teuchos::null
		);
	
	///
	/** Initialize given the <tt>Epetra_multi_vec</tt>.
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_MultiVector>         &epetra_multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_range       = Teuchos::null
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_domain      = Teuchos::null
		);
	
	/// Set uninitalized.
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_MultiVector>        *epetra_multi_vec = NULL
		,Teuchos::RefCountPtr<const EpetraVectorSpace>  *epetra_range     = NULL
		,Teuchos::RefCountPtr<const EpetraVectorSpace>  *epetra_domain    = NULL
		);

	//@}

	/** @name Overridden from OpBase */
	//@{
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > domain() const;
	//@}

	/** @name Overridden from LinearOp */
	//@{
	///
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;
	//@}

	/** @name Overridden from MultiVector */
	//@{
	///
	Teuchos::RefCountPtr<Vector<Scalar> > col(Index j);
	///
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const Range1D& col_rng );
	///
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const int numCols, const int cols[] );
	//@}

	/** @name Overridden from MPIMultiVectorBase */
	//@{
	///
	Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const;
	///
	void getLocalData( Scalar **values, Index *leadingDim );
	//@}

private:
	
#ifdef DOXYGEN_COMPILE
	Epetra_MultiVector                              *epetra_multi_vec;
	EpetraVectorSpace                               *epetra_range;
	EpetraVectorSpace                               *epetra_domain;
#else
	Teuchos::RefCountPtr<Epetra_MultiVector>        epetra_multi_vec_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>   epetra_range_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>   epetra_domain_;
#endif
	
}; // end class EpetraMultiVector

} // end namespace TSFCore

#endif // TSFCORE_EPETRA_MULTI_VECTOR_HPP
