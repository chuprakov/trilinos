// /////////////////////////////////////////////////////////////////
// TSFCoreEpetraVector.hpp

#ifndef TSFCORE_EPETRA_VECTOR_HPP
#define TSFCORE_EPETRA_VECTOR_HPP

#include "TSFCoreMPIVectorBase.hpp"

///
class Epetra_Vector;

namespace TSFCore {

///
class EpetraVectorSpace;

///
/** Wrapper for Epetra vectors.
 *
 * The implementation of this class is quite trivial once you
 * understand the implications of the base class
 * <tt>MPIVectorBase</tt> that it derives from.  This vector
 * implementation will be seamlessly compatible with every other
 * MPI-based vector implementation.
 *
 * Note that an <tt>EpetraVectorSpace</tt> object is returned from
 * <tt>mpiSpace()</tt>.
 *
 * Note: the default copy constructor is allowed but the default
 * assignment operator is not.
 */
class EpetraVector : public MPIVectorBase<RTOp_value_type> {
public:

	///
	typedef RTOp_value_type Scalar;

	/** @name Constructors/initializers */
	//@{

	///
	/** Construct to uninitialized.
	 *
	 * Has the same postconditions as <tt>setUninitialized()</tt>.
	 */
	EpetraVector();

	///
	/** Calls <tt>this->initialize()</tt>.
	 */
	EpetraVector(
		const Teuchos::RefCountPtr<Epetra_Vector>              &epetra_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_vec_spc
		);

	///
	/** Initialize given <tt>Epetra_Vector</tt> and <tt>EpetraVectorSpace</tt> objects.
	 *
	 * @param  epetra_vec
   *              [in] Smart pointer to the <tt>Epetra_Vector</tt> object that
	 *              <tt>*this</tt> wraps.
	 * @param  epetra_vec_spc
   *              [in] Smart pointer to the <tt>EpetraVectorSpace</tt> object that
	 *              wraps a <tt>Epetra_Map</tt> defining this vector space.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>epetra_vec.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * <li> <tt>epetra_vec_spc.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_vec().get() == epetra_vec.get()</tt>
	 * <li> <tt>this->mpiSpace().get() != NULL</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<Epetra_Vector>              &epetra_vec
		,const Teuchos::RefCountPtr<const EpetraVectorSpace>   &epetra_vec_spc
		);

	///
	/** Set uninitialized and return the (possibly) owned
	 * <tt>Epetra_Vector</tt> object.
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_vec().get() == NULL</tt>.
	 * <li> <tt>this->mpiSpace().get() == NULL</tt>
	 * </ul>
	 */
	void setUninitialized(
		Teuchos::RefCountPtr<Epetra_Vector>              *epetra_vec     = NULL
		,Teuchos::RefCountPtr<const EpetraVectorSpace>   *epetra_vec_spc = NULL
		);

	///
	/** Return a smart pointer object to a <tt>const</tt> <tt>Epetra_Vector</tt> object.
	 */
	Teuchos::RefCountPtr<const Epetra_Vector> epetra_vec() const;
	
	//@}

	/** @name Overridden from MPIVectorBase */
	//@{

	///
	Teuchos::RefCountPtr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const;
	///
	void getLocalData( Scalar** values, ptrdiff_t* stride );

	//@}

private:

	// ///////////////////////////////////////
	// Private data members

	Teuchos::RefCountPtr<Epetra_Vector>            epetra_vec_;
	Teuchos::RefCountPtr<const EpetraVectorSpace>  epetra_vec_spc_;

	// ///////////////////////////////////////
	// Private member functions

	// Not defined and not to be called
	EpetraVector& operator=(const EpetraVector&);

}; // end class EpetraVector

// //////////////////////////////////////
// Inline members

inline
Teuchos::RefCountPtr<const Epetra_Vector>
EpetraVector::epetra_vec() const
{
	return epetra_vec_;
}

} // end namespace TSFCore

#endif // TSFCORE_EPETRA_VECTOR_HPP
