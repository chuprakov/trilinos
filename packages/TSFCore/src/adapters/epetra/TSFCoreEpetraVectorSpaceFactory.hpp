// ///////////////////////////////////////////////////////////////
// TSFCoreEpetraVectorSpaceFactory.hpp

#ifndef TSFCORE_EPETRA_VECTOR_SPACE_FACTORY_HPP
#define TSFCORE_EPETRA_VECTOR_SPACE_FACTORY_HPP

#include "TSFCoreEpetraTypes.hpp"
#include "TSFCoreVectorSpaceFactory.hpp"

namespace TSFCore {

///
/** Subclass of <tt>VectorSpaceFactory</tt> that creates <tt>EpetraVectorSpace</tt>
 * objects given a dimension.
 *
 * ToDo: Finish documentation!
 */
class EpetraVectorSpaceFactory : public VectorSpaceFactory<double> {
public:

	///
	typedef double Scalar;

	/** @name Constructors / initializers */
	//@{

	///
	/** Constructs to an uninitialized state.
	 *
	 * See the postconditions for <tt>setUninitialized()</tt>.
	 */
	EpetraVectorSpaceFactory();

	///
	/** Calls <tt>initialize()</tt>.
	 */
	EpetraVectorSpaceFactory(
		const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
		);

	///
	/** Initialize given a communicator.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>epetra_comm.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_com.get()==epetra_comm.get()</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<const Epetra_Comm>  &epetra_comm
		);

	///
	/** Set uninitialized and return the underlying <tt>Epetra_BlockMap</tt>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_map.get()==NULL</tt>
	 * </ul>
	 */
	void setUninitialized(
		Teuchos::RefCountPtr<const Epetra_Comm> *epetra_comm = NULL
		);

	///
	/** Return a smart pointer to the underlying <tt>Epetra_BlockMap</tt> object.
	 */
	Teuchos::RefCountPtr<const Epetra_Comm> epetra_comm() const;

	//@}

	/** @name Overridden from VectorSpaceFactory */
	//@{

	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > createVecSpc(int dim) const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	Epetra_Comm                                     *epetra_comm;
#else	
	Teuchos::RefCountPtr<const Epetra_Comm>         epetra_comm_;
#endif

}; // end class EpetraVectorSpaceFactory

} // end namespace TSFCore

#endif  // TSFCORE_EPETRA_VECTOR_SPACE_FACTORY_HPP
