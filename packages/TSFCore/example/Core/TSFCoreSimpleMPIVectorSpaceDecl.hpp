// ////////////////////////////////////////////////////////////////////////////
// TSFCoreSimpleMPIVectorSpaceDecl.hpp

#ifndef TSFCORE_SIMPLE_MPI_VECTOR_SPACE_DECL_HPP
#define TSFCORE_SIMPLE_MPI_VECTOR_SPACE_DECL_HPP

#include "TSFExtended/src/Core/TSFCoreMPIVectorSpaceBase.hpp"

namespace TSFCore {

///
/** Simple MPI-based vector class.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class SimpleMPIVectorSpace : virtual public MPIVectorSpaceBase<Scalar> {
public:

	///
	/** Constructs to a parallel vector with a fixed number of elements per process.
	 *
	 * @param  mpiComm     [in] The MPI communicator.  This object must be maintained
	 *                     by the client the entire time that <tt>this</tt> is in use.
	 * @pram  localSubDim  [in] The number of elements per processor.  This number must
	 *                     be the same on every processor in this simple implementation.
	 */
	SimpleMPIVectorSpace( MPI_Comm mpiComm, const Index localSubDim );

	/** @name Overridden from VectorSpece */
	//@{

	///
	Index dim() const;
	///
	Teuchos::RefCountPtr<Vector<Scalar> > createMember() const;
	///
	Teuchos::RefCountPtr< const VectorSpace<Scalar> > clone() const;

	//@}

	/** @name Overridden from MPIVectorSpaceBase */
	//@{

	///
	MPI_Comm mpiComm() const;
	///
	Index localOffset() const;
	///
 	Index localSubDim() const;

	//@}

private:

	// //////////////////////////////////////
	// Private data members

	MPI_Comm           mpiComm_;
	Index              localSubDim_;
	Index              globalDim_;
	Index              localOffset_;
	int                numProc_;
	int                procRank_;

	// Not defined and not to be called
	SimpleMPIVectorSpace();
	
}; // end class SimpleMPIVectorSpace

} // end namespace TSFCore

#endif // TSFCORE_SIMPLE_MPI_VECTOR_SPACE_DECL_HPP
