// ///////////////////////////////////////////////////////////////
// TSFCoreMultiVectorAllocator.hpp

#ifndef TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP
#define TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP

#include "TSFCore_ConfigDefs.hpp"
#include "TSFCoreVectorSpace.hpp"
#include "ThrowException.hpp"

namespace TSFCore {

///
/** Allocator class to be used with <tt>MemMngPack::AbstractFactoryStd</tt> to create
 * <tt>MultiVector</tt> objects of a given size.
 */
template<class Scalar>
class MultiVectorAllocator {
public:
	///
	MultiVectorAllocator() : numMembers_(0) {}
	///
	typedef MemMngPack::ref_count_ptr<MultiVector<Scalar> >  ptr_t;         // required!
	///
	MultiVectorAllocator( const MemMngPack::ref_count_ptr<const VectorSpace<Scalar> > &vs, int numMembers )
		: vs_(vs), numMembers_(numMembers)
		{
#ifdef _DEBUG
			THROW_EXCEPTION( vs.get()==NULL, std::logic_error, "Error!" );
#endif			
		}
	///
	const ptr_t allocate() const { return vs_->createMembers(numMembers_); }  // required!
private:
	MemMngPack::ref_count_ptr<const VectorSpace<Scalar> >      vs_;
	int                                                        numMembers_;
};

} // namespace TSFCore

#endif // TSFCORE_MULTI_VECTOR_ALLOCATOR_HPP
