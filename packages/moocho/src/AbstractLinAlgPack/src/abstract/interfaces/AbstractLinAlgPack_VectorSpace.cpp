// //////////////////////////////////////////////////////////////////////
// VectorSpace.cpp

#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorSpaceSubSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

VectorSpace::space_ptr_t
VectorSpace::sub_space(const Range1D& rng_in) const
{
	namespace rcp = ReferenceCountingPack;
	const index_type dim = this->dim();
	const Range1D    rng = rng_in.full_range() ? Range1D(1,dim) : rng_in;
#ifdef _DEBUG
	THROW_EXCEPTION(
		rng.ubound() > dim, std::out_of_range
		,"VectorSpace::sub_space(rng): Error, rng = ["<<rng.lbound()<<","<<rng.ubound()<<"] "
		"is not in the range [1,this->dim()] = [1,"<<dim<<"]" );
#endif	
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return space_ptr_t( this, false );
	return rcp::rcp_implicit_cast<space_ptr_t::element_type>(
		rcp::ref_count_ptr<VectorSpaceSubSpace>(
			new VectorSpaceSubSpace(
				space_ptr_t( this, false )
				,rng ) ) );
}

// VectorSpace

VectorSpaceBase::vec_ptr_t  VectorSpace::new_member() const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp_implicit_cast<vec_ptr_t::element_type>(create_member());
}

} // end namespace AbstractLinAlgPack
