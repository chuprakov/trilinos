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
/** Default subclass for <tt>MultiVector</tt> implemented using columns
 * of separate abstract vectors.
 *
 * This is a very bad implementation of a multi-vector but this will
 * work in situations where you need a multi-vector but some
 * underlying linear algebra library does not directly support them.
 *
 * This subclass can be used to represent a <tt>%MultiVector</tt>
 * wrapper around a single <tt>Vector</tt> object so that a single
 * vector can be passed to a method that expects a <tt>%MultiVector</tt>
 * object.
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
		const Teuchos::RefCountPtr<Epetra_MultiVector>         multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &range
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &domain   = Teuchos::null
		);
	
	///
	/** Initialize given the spaces for the columns and rows.
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_MultiVector>         multi_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &range
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &domain   = Teuchos::null
		);
	
	/// Set uninitalized.
	void setUninitialized();

	//@}

	/** @name Overridden from LinearOp */
	//@{
	///
	Teuchos::RefCountPtr<const VectorSpace<Scalar> > domain() const;
	//@}

	/** @name Overridden from MultiVector */
	//@{
	///
	Teuchos::RefCountPtr<Vector<Scalar> > col(Index j);
	///
	Teuchos::RefCountPtr<MultiVector<Scalar> > subView( const Range1D& col_rng );
	//@}

	/** @name Overridden from MPIMultiVectorBase */
	//@{
	///
	Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const;
	//@}

private:
	
#ifdef DOXYGEN_COMPILE
	const EpetraVectorSpace                                *range;
	const EpetraVectorSpace                                *domain;
#else
	Teuchos::RefCountPtr<Epetra_MultiVector>          multi_vec_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>     range_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>     domain_;
#endif
	
}; // end class EpetraMultiVector

} // end namespace TSFCore

#endif // TSFCORE_EPETRA_MULTI_VECTOR_HPP
