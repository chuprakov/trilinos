// ////////////////////////////////////////////////////////////////
// VectorWithOpMutableSubView.cpp

#include "AbstractLinAlgPack/include/VectorWithOpMutableSubView.h"
#include "ThrowException.h"

namespace AbstractLinAlgPack {

VectorWithOpMutableSubView::VectorWithOpMutableSubView( const vec_mut_ptr_t& vec, const Range1D& rng )
	: VectorWithOpSubView(NULL,rng)
{
	this->initialize(vec,rng);
}

void VectorWithOpMutableSubView::initialize( const vec_mut_ptr_t& vec, const Range1D& rng )
{
	namespace rcp = ReferenceCountingPack;
	VectorWithOpSubView::initialize(
		rcp::rcp_implicit_cast<vec_ptr_t::element_type>(vec)
		,rng );
	full_vec_ = vec;
	this->has_changed();
}

// Overriddend form VectorWithOp

VectorWithOp::vec_ptr_t VectorWithOpMutableSubView::sub_view( const Range1D& rng ) const
{
	return VectorWithOpSubView::sub_view(rng); // Had to override to resolve conflicit!
}

// Overridden from VectorWithOpMutable

void VectorWithOpMutableSubView::apply_transformation(
	const RTOpPack::RTOp& op
	,const size_t num_vecs, const VectorWithOp** vecs
	,const size_t num_targ_vecs, VectorWithOpMutable** targ_vecs
	,RTOp_ReductTarget reduct_obj
	,const index_type global_offset_in, const index_type sub_dim_in
	)
{
#ifdef _DEBUG
	const index_type this_dim = this->dim();
	THROW_EXCEPTION(
		sub_dim_in < 0
		|| ( global_offset_in < 0 && -global_offset_in > this_dim )
		|| ( ( global_offset_in < 0 && sub_dim_in > 0 ) && -global_offset_in + sub_dim_in > this_dim )
		|| ( ( global_offset_in >= 0 && sub_dim_in > 0  ) && sub_dim_in > this_dim )
		, std::logic_error
		,"VectorWithOpMutableSubView::apply_transformation(...): Error, global_offset_in = "
		<< global_offset_in << ", sub_dim_in = " << sub_dim_in << " and this->dim() = this_dim  = "
		<< this_dim << " are not compatible." );
#endif
	const index_type this_offset = space_impl().rng().lbound() - 1;
	const index_type
		this_sub_dim = ( sub_dim_in 
						 ? sub_dim_in
						 : ( space_impl().rng().size()
							 + ( global_offset_in >= 0
								 ? 0 : global_offset_in )
							 ) );
	full_vec_->apply_transformation(
		op, num_vecs, vecs, num_targ_vecs, targ_vecs, reduct_obj
		, global_offset_in - this_offset, this_sub_dim );
}

void VectorWithOpMutableSubView::set_ele( index_type i, value_type val )
{
	space_impl().validate_range(Range1D(i,i));
	return full_vec_->set_ele( space_impl().rng().lbound() + i - 1, val );
}

VectorWithOpMutable::vec_mut_ptr_t
VectorWithOpMutableSubView::sub_view( const Range1D& rng )
{
	namespace rcp = ReferenceCountingPack;
	space_impl().validate_range(rng);
	const index_type this_offset = space_impl().rng().lbound() - 1;
	return rcp::rcp_implicit_cast<vec_mut_ptr_t::element_type>(
		rcp::ref_count_ptr<VectorWithOpMutableSubView>(
			new VectorWithOpMutableSubView(
				full_vec_
				,Range1D( 
					this_offset  + rng.lbound()
					,this_offset + rng.ubound() )
				) ) );
}

void VectorWithOpMutableSubView::set_sub_vector( const RTOp_SubVector& sub_vec_in )
{
	const index_type  this_offset = space_impl().rng().lbound() - 1;
	RTOp_SubVector    sub_vec = sub_vec_in;
	sub_vec.global_offset += this_offset;
	full_vec_->set_sub_vector( sub_vec );
}

} // end namespace AbstractLinAlgPack
