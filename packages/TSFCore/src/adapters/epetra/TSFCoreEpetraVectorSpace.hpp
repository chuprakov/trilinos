// ////////////////////////////////////////////////////////////////////////////
// TSFCoreEpetraVectorSpace.hpp

#ifndef TSFCORE_EPETRA_VECTOR_SPACE_HPP
#define TSFCORE_EPETRA_VECTOR_SPACE_HPP

#include "TSFCoreEpetraTypes.hpp"
#include "TSFCoreMPIVectorSpaceBase.hpp"

namespace TSFCore {

///
/** <tt>%VectorSpace</tt> Subclass for Epetra vectors and multi-vectors.
 *
 * This uses an <tt>Epetra_BlockMap</tt> object to implement the
 * <tt>MPIVectorSpaceBase</tt> interface.  By implementing the
 * <tt>MPIVectorSpaceBase</tt> interface, this implementation allows
 * the seemless collaboration of different vectors and multi-vectors
 * through the interface classes <tt>MPIVectorBase</tt> and
 * <tt>MPIMultiVectorBase</tt>.
 *
 * This class works properly even if Epetra is not compiled with
 * support for MPI (i.e. <tt>PETRA_COMM_MPI</tt> is not defined when
 * compiling and linking).  If MPI support is not compiled into
 * Epetra, then the dummy implementation defined in
 * <tt>RTOp_mpi.h</tt> is used instead.
 *
 * The default copy constructor and assignment operators are allowed
 * since they have the correct semantics.
 */
class EpetraVectorSpace : public MPIVectorSpaceBase<RTOp_value_type> {
public:

	///
	typedef RTOp_value_type Scalar;

	/** @name Constructors / initializers */
	//@{

	///
	/** Constructs to an uninitialized state.
	 *
	 * See the postconditions for <tt>setUninitialized()</tt>.
	 */
	EpetraVectorSpace();

	///
	/** Calls <tt>initialize()</tt>.
	 */
	EpetraVectorSpace(
		const Teuchos::RefCountPtr<const Epetra_BlockMap>  &epetra_map
		);

	///
	/** Initialize given an <tt>Epetra_BlockMap</tt>.
	 *
	 * Preconditions:<ul>
	 * <li> <tt>epetra_map.get()!=NULL</tt> (throw <tt>std::invalid_argument</tt>)
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_map.get()==epetra_map.get()</tt>
	 * <li> <tt>this->dim()==epetra_map->NumGlobalElements()</tt>
	 * <li> if (<tt>PETRA_COMM_MPI</tt> is defined and <tt>dynamic_cast<const Epetra_MpiComm*>(epetra_map.get()) != NULL</tt>)
	 *      <ul><li><tt>this->mpiComm() == dynamic_cast<const Epetra_MpiComm&>(epetra_map)->Comm()</tt></ul>
	 *      else
	 *      <ul><li><tt>this->mpiComm() == MPI_COMM_NULL</tt></ul>
	 * <li> <tt>this->localOffset() == epetra_map->MinMyGID() - epetra_map->IndexBase()</tt>
	 * <li> <tt>this->localSubdim() == epetra_map->NumMyElements()</tt>
	 * </ul>
	 */
	void initialize(
		const Teuchos::RefCountPtr<const Epetra_BlockMap>  &epetra_map
		);

	///
	/** Set uninitialized and return the underlying <tt>Epetra_BlockMap</tt>
	 *
	 * Postconditions:<ul>
	 * <li> <tt>this->epetra_map.get()==NULL</tt>
	 * <li> <tt>this->dim()==0</tt>
	 * <li> <tt>this->mpiComm() == MPI_COMM_NULL</tt>
	 * <li> <tt>this->localOffset() == -1</tt>
	 * <li> <tt>this->localSubdim() == -1)</tt>
	 * </ul>
	 */
	void setUninitialized(
		Teuchos::RefCountPtr<const Epetra_BlockMap> *epetra_map = NULL
		);

	///
	/** Return a smart pointer to the underlying <tt>Epetra_BlockMap</tt> object.
	 */
	Teuchos::RefCountPtr<const Epetra_BlockMap> epetra_map() const;

	//@}

	/** @name Overridden from VectorSpece */
	//@{

	/// Returns 0 if uninitialized.
	Index dim() const;
	/// Returns an allocated <tt>EpetraVector</tt> object. 
	Teuchos::RefCountPtr<Vector<Scalar> > createMember() const;
	/// Returns a <tt>EpetraVectorSpaceFactory</tt> object.
	Teuchos::RefCountPtr< const VectorSpaceFactory<Scalar> > smallVecSpcFcty() const;
	/// Returns an allocated <tt>EpetraMultiVector</tt> object. 
	Teuchos::RefCountPtr< MultiVector<Scalar> > createMembers(int numMembers) const;
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > clone() const;

	//@}

	/** @name Overriddend from MPIVectorSpaceBase */
	//@{
	
	///
	MPI_Comm mpiComm() const;
	///
	Index localSubDim() const;

	//@}

private:

#ifdef DOXYGEN_COMPILE
	Epetra_BlockMap                                        *epetra_map;
	EpetraVectorSpaceFactory                               *smallVecSpcFcty;
#else	
	Teuchos::RefCountPtr<const Epetra_BlockMap>            epetra_map_;
	MPI_Comm                                               mpiComm_;
	Index                                                  localSubDim_;
	Teuchos::RefCountPtr<const EpetraVectorSpaceFactory>   smallVecSpcFcty_;
#endif

}; // end class EpetraVectorSpace

// //////////////////////////////////////////////
// Inline members

inline
Teuchos::RefCountPtr<const Epetra_BlockMap>
EpetraVectorSpace::epetra_map() const
{
	return epetra_map_;
}

} // end namespace TSFCore

#endif // TSFCORE_EPETRA_VECTOR_SPACE_HPP
