// /////////////////////////////////////////////////////////////////////////
// VectorWithOpSubView.cpp

#include <assert.h>

#include <stdexcept>

#include "AbstractLinAlgPack/include/VectorWithOpSubView.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

VectorWithOpSubView::VectorWithOpSubView( const vec_ptr_t& vec, const Range1D& rng )
	: space_(NULL,rng)
{
	initialize(vec,rng);
}

void VectorWithOpSubView::initialize( const vec_ptr_t& vec, const Range1D& rng )
{
	typedef VectorSpace::space_ptr_t   space_ptr_t;
	space_.initialize(
		vec.get() ? space_ptr_t(&vec->space(),false) : space_ptr_t( NULL )
		,rng
		);
	full_vec_ = vec;
	this->has_changed();
}

// Overridden from VectorWithOp

const VectorSpace& VectorWithOpSubView::space() const
{
	return space_;
}

index_type VectorWithOpSubView::dim() const
{
	return space_.dim();
}

void VectorWithOpSubView::apply_reduction(
	const RTOpPack::RTOp& op
	,const size_t num_vecs, const VectorWithOp** vecs
	,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type first_ele_in, const index_type sub_dim_in, const index_type global_offset_in
	) const
{
#ifdef _DEBUG
	const index_type this_dim = this->dim();
	THROW_EXCEPTION(
		sub_dim_in < 0
		|| !(1 <= first_ele_in && first_ele_in <= this_dim)
		|| ( sub_dim_in > 0 && (sub_dim_in - (first_ele_in - 1) > this_dim) )
		, std::logic_error
		,"VectorWithOpSubView::apply_reduction(...): Error, first_ele_in = "
		<< first_ele_in << ", global_offset_in = " << global_offset_in
		<< ", sub_dim_in = " << sub_dim_in << " and this->dim() = this_dim  = "
		<< this_dim << " are not compatible." );
#endif
	const index_type this_offset = space_impl().rng().lbound() - 1;
	const index_type
		this_sub_dim = ( sub_dim_in 
						 ? sub_dim_in
						 : space_impl().rng().size() - (first_ele_in - 1)
			           );
	full_vec_->apply_reduction(
		op, num_vecs, vecs, num_targ_vecs, targ_vecs, reduct_obj
		,this_offset + first_ele_in     // first_ele
		,this_sub_dim                   // sub_dim
		,global_offset_in               // global_dim
		);
}

value_type VectorWithOpSubView::get_ele(index_type i) const
{
	space_.validate_range(Range1D(i,i));
	return full_vec_->get_ele( space_.rng().lbound() + i - 1 );
}

VectorWithOp::vec_ptr_t
VectorWithOpSubView::sub_view( const Range1D& rng ) const
{
	namespace rcp = ReferenceCountingPack;
	space_.validate_range(rng);
	const index_type this_offset = space_.rng().lbound() - 1;
	return rcp::rcp_implicit_cast<vec_ptr_t::element_type>(
		rcp::ref_count_ptr<VectorWithOpSubView>(
			new VectorWithOpSubView(
				full_vec_
				,Range1D( 
					this_offset  + rng.lbound()
					,this_offset + rng.ubound() )
				) ) );
}

void VectorWithOpSubView::get_sub_vector(
	const Range1D& rng, ESparseOrDense sparse_or_dense, RTOp_SubVector* sub_vec ) const
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		!sub_vec, std::logic_error
		,"VectorWithOpSubView::get_sub_vector(...): Error!" ) ;
#endif
	space_.validate_range(rng);
	const index_type this_offset = space_.rng().lbound() - 1;
	full_vec_->get_sub_vector( rng + this_offset, sparse_or_dense, sub_vec );
	sub_vec->global_offset -= this_offset;
}

void VectorWithOpSubView::free_sub_vector( RTOp_SubVector* sub_vec ) const
{
#ifdef _DEBUG
	THROW_EXCEPTION(
		!sub_vec, std::logic_error
		,"VectorWithOpSubView::free_sub_vector(...): Error!" ) ;
#endif
	const index_type this_offset = space_.rng().lbound() - 1;
	sub_vec->global_offset += this_offset;
	full_vec_->free_sub_vector( sub_vec );
}

} // end namespace AbstractLinAlgPack
