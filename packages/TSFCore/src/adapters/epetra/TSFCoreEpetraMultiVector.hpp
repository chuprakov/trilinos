// /////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraMultiVector.hpp

#ifndef TSFCORE_EPETRA_MULTI_VECTOR_HPP
#define TSFCORE_EPETRA_MULTI_VECTOR_HPP

#include <vector>

#include "TSFCoreMPIMultiVectorBase.hpp"

// Define this to use Epetra_MultiVector::Multiply(...) to implement apply(...)
//#define TSFCORE_EPETRA_USE_EPETRA_MULTI_VECTOR_MULTIPLY

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
	 * <li> <tt>this->range().get() == NULL</tt>
	 * <li> <tt>this->domain().get() == NULL</tt>
	 * </ul>
	 */
	EpetraMultiVector();

	/// Calls <tt>initalize()</tt>.
	EpetraMultiVector(
		const Teuchos::RefCountPtr<Epetra_MultiVector>         &epetra_multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_range
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_domain      = Teuchos::null
		);
	
	///
	/** Initialize given the <tt>Epetra_multi_vec</tt>.
   *
   * ToDo: Finish documentation!
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_MultiVector>         &epetra_multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_range
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_domain      = Teuchos::null
		);
	
	///
  /** Set uninitalized.
   *
   * ToDo: Finish documentation!
   */
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_MultiVector>        *epetra_multi_vec = NULL
		,Teuchos::RefCountPtr<const EpetraVectorSpace>  *epetra_range     = NULL
		,Teuchos::RefCountPtr<const EpetraVectorSpace>  *epetra_domain    = NULL
		);

	///
	/** Get a smart pointer to non-<tt>const</tt> <tt>Epetra_MultiVector</tt> object.
	 *
	 * Note: Only the numerical values of the Epetra vector should be changed by the
	 * client.  Any other type of change will invalidate <tt>this</tt>.
	 */
 	Teuchos::RefCountPtr<Epetra_MultiVector> epetra_multi_vec();

	///
	/** Get a smart pointer to <tt>const</tt> <tt>Epetra_MultiVector</tt> object.
	 *
	 * This gives access to read the <tt>Epetra_MultiVector</tt> object only
	 * so this should be very safe.
	 */
	Teuchos::RefCountPtr<const Epetra_MultiVector> epetra_multi_vec() const;

	//@}

	/** @name Overridden from OpBase */
	//@{
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > domain() const;
	//@}

	/** @name Overridden from LinearOp */
	//@{
#ifdef TSFCORE_EPETRA_USE_EPETRA_MULTI_VECTOR_MULTIPLY
	///
	void apply(
		const ETransp                 M_trans
		,const MultiVector<Scalar>    &X
		,MultiVector<Scalar>          *Y
		,const Scalar                 alpha
		,const Scalar                 beta
		) const;
#endif
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
	void getLocalData( const Scalar **values, Index *leadingDim ) const;
	///
	void freeLocalData( const Scalar *values ) const;
	///
	void getLocalData( Scalar **values, Index *leadingDim );
	///
	void commitLocalData( Scalar *values );
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

// ///////////////////////////////////
// Inline members

inline
Teuchos::RefCountPtr<Epetra_MultiVector>
EpetraMultiVector::epetra_multi_vec()
{
	return epetra_multi_vec_;
}

inline
Teuchos::RefCountPtr<const Epetra_MultiVector>
EpetraMultiVector::epetra_multi_vec() const
{
	return epetra_multi_vec_;
}

} // end namespace TSFCore

#endif // TSFCORE_EPETRA_MULTI_VECTOR_HPP
