// /////////////////////////////////////////////////////////////////
// TSFCoreSimpleMPIVectorDecl.hpp

#ifndef TSFCORE_SIMPLE_MPI_VECTOR_DECL_HPP
#define TSFCORE_SIMPLE_MPI_VECTOR_DECL_HPP

#include <vector>

#include "TSFExtended/src/Core/TSFCoreMPIVectorBase.hpp"

namespace TSFCore {

///
/** Simple implementation of an MPI-base parallel vector.
 *
 * ToDo: Finish documentation
 */
template<class Scalar>
class SimpleMPIVector : virtual public MPIVectorBase<Scalar> {
public:

	///
	SimpleMPIVector( const MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<Scalar> > &mpiSpace );

	/** @name Overridden from MPIVectorBase */
	//@{

	///
	MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<Scalar> > mpiSpace() const;
	///
	void getLocalData( Scalar** values, ptrdiff_t* stride );

	//@}

private:

	// ///////////////////////////////////////
	// Private data members
	
	MemMngPack::ref_count_ptr<const MPIVectorSpaceBase<Scalar> >   mpiSpace_;
	std::vector<Scalar>                                            localValues_;

	// not defined and not to be called
	SimpleMPIVector();

}; // end class SimpleMPIVector

} // end namespace TSFCore

#endif // TSFCORE_SIMPLE_MPI_VECTOR_DECL_HPP
