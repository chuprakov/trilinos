// //////////////////////////////////////////////////////////////////////
// VectorSpaceSubSpace.cpp

#include <assert.h>

#include "AbstractLinAlgPack/include/VectorSpaceSubSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutableSubView.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

VectorSpaceSubSpace::VectorSpaceSubSpace( const space_ptr_t& full_space, const Range1D& rng )
{
	this->initialize(full_space,rng);
}

void VectorSpaceSubSpace::initialize( const space_ptr_t& full_space, const Range1D& rng )
{
#ifdef _DEBUG
	const index_type n = full_space.get() ? full_space->dim() : 0;
	THROW_EXCEPTION(
		full_space.get() && !rng.full_range() && rng.ubound() > n, std::out_of_range
		,"VectorSpaceSubSpace::initialize(...): Error, "
		"rng = [" << rng.lbound() << "," << rng.ubound() << "] is not in the range "
		"[1,vec->dim()] = [1," << n << "]" );
#endif
	full_space_ = full_space;
	rng_        = full_space.get() && rng.full_range() ? Range1D(1,full_space->dim()) : rng;
}

#ifdef _DEBUG
void VectorSpaceSubSpace::validate_range(const Range1D& rng) const
{
	const index_type n = this->dim();
	THROW_EXCEPTION(
		full_space_.get() == NULL, std::logic_error
		,"VectorSpaceSubSpace::validate_range(rng): Error, Uninitialized" );
	THROW_EXCEPTION(
		full_space_.get() && !rng.full_range() && rng.ubound() > n, std::logic_error
		,"VectorSpaceSubSpace::validate_range(rng): Error, "
		"rng = [" << rng.lbound() << "," << rng.ubound() << "] is not in the range "
		"[1,this->dim] = [1," << n << "]" );
}
#endif

// Overridden from VectorSpaceBase


bool VectorSpaceSubSpace::is_compatible(const VectorSpaceBase& another_space) const
{
	const VectorSpaceSubSpace
		*a_space = dynamic_cast<const VectorSpaceSubSpace*>(&another_space);
	if(!a_space)
		return false;
	return this->rng_ == a_space->rng_ && this->full_space_->is_compatible(*a_space->full_space_);
}

// Overridden form VectorSpace

index_type VectorSpaceSubSpace::dim() const
{
	return full_space_.get() ? rng_.size() : 0;
}

VectorSpace::vec_mut_ptr_t VectorSpaceSubSpace::create_member() const
{
	namespace rcp = ReferenceCountingPack;
	return rcp::rcp_implicit_cast<vec_mut_ptr_t::element_type>(
		rcp::ref_count_ptr<VectorWithOpMutableSubView>(
			new VectorWithOpMutableSubView(
				full_space_->create_member(), rng_ 
				) ) );
}

VectorSpace::space_ptr_t VectorSpaceSubSpace::sub_space(const Range1D& rng_in) const
{
	namespace rcp = ReferenceCountingPack;
	validate_range(rng_in);
	const index_type dim         = this->dim();
	const Range1D    rng         = rng_in.full_range() ? Range1D(1,dim) : rng_in;
	if( rng.lbound() == 1 && rng.ubound() == dim )
		return space_ptr_t( this, false );
	const index_type this_offset = rng_.lbound() - 1;
	return rcp::rcp_implicit_cast<space_ptr_t::element_type>(
		rcp::ref_count_ptr<VectorSpaceSubSpace>(
			new VectorSpaceSubSpace(
				full_space_
				,Range1D( 
					this_offset  + rng.lbound()
					,this_offset + rng.ubound() )
				) ) );
}

} // end namespace AbstractLinAlgPack
