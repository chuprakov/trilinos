// ///////////////////////////////////////////////////////////////////
// BasisSystemCompositeStd.cpp
//
// Copyright (C) 2001 Roscoe Ainsworth Bartlett
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the "Artistic License" (see the web site
//   http://www.opensource.org/licenses/artistic-license.html).
// This license is spelled out in the file COPYING.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// above mentioned "Artistic License" for more details.

#include <assert.h>

#include "AbstractLinAlgPack/src/BasisSystemCompositeStd.hpp"
#include "AbstractLinAlgPack/src/MatrixWithOpNonsingular.hpp"
#include "AbstractLinAlgPack/src/MatrixCompositeStd.hpp"
#include "AbstractLinAlgPack/src/MultiVectorMutable.hpp"
#include "AbstractLinAlgPack/src/VectorSpaceCompositeStd.hpp"
#include "AbstractLinAlgPack/src/LinAlgOpPack.hpp"
#include "ReleaseResource_ref_count_ptr.hpp"
#include "AbstractFactoryStd.hpp"
#include "dynamic_cast_verbose.hpp"
#include "ThrowException.hpp"

namespace {

// Allocates a MultiVectorMutable object given a vector space
// object and the number of columns (num_vecs).

class AllocatorMultiVectorMutable {
public:
	AllocatorMultiVectorMutable(
		const AbstractLinAlgPack::VectorSpace::space_ptr_t&  vec_space
		,AbstractLinAlgPack::size_type                       num_vecs
		)
		:vec_space_(vec_space)
		 ,num_vecs_(num_vecs)
	{}
	typedef MemMngPack::ref_count_ptr<
		AbstractLinAlgPack::MultiVectorMutable>               ptr_t;
	ptr_t allocate() const
	{
		return vec_space_->create_members(num_vecs_);
	}
private:
	AbstractLinAlgPack::VectorSpace::space_ptr_t  vec_space_;
	AbstractLinAlgPack::size_type                 num_vecs_;
}; // end class AllocatorMultiVectorMutable

} // end namespace

namespace AbstractLinAlgPack {

// Static member functions

void BasisSystemCompositeStd::initialize_space_x(
	const VectorSpace::space_ptr_t    &space_xD
	,const VectorSpace::space_ptr_t   &space_xI
	,Range1D                          *var_dep
	,Range1D                          *var_indep
	,VectorSpace::space_ptr_t         *space_x
	)
{
	namespace rcp = MemMngPack;
#ifdef _DEBUG
	THROW_EXCEPTION(
		space_xD.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_space_x(...): Error!" );
	THROW_EXCEPTION(
		space_xI.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_space_x(...): Error!" );
	THROW_EXCEPTION(
	    var_dep == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_space_x(...): Error!" );
	THROW_EXCEPTION(
	    var_indep == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_space_x(...): Error!" );
	THROW_EXCEPTION(
		space_x == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_space_x(...): Error!" );
#endif
	*var_dep   = Range1D(1,space_xD->dim());
	*var_indep = Range1D(var_dep->ubound()+1,var_dep->ubound()+space_xI->dim());
	const VectorSpace::space_ptr_t
		vec_spaces[2] = { space_xD, space_xI };
	*space_x   = rcp::rcp(new VectorSpaceCompositeStd(vec_spaces,2));
}

const BasisSystemCompositeStd::fcty_Gc_ptr_t
BasisSystemCompositeStd::factory_Gc()
{
	namespace rcp = MemMngPack;
	return rcp::rcp( new MemMngPack::AbstractFactoryStd<MatrixWithOp,MatrixCompositeStd>() );
}

void BasisSystemCompositeStd::initialize_Gc(
	const VectorSpace::space_ptr_t    &space_x
	,const Range1D                    &var_dep
	,const Range1D                    &var_indep
	,const VectorSpace::space_ptr_t   &space_c
	,const C_ptr_t                    &C
	,const N_ptr_t                    &N
	,MatrixWithOp                     *Gc
	)
{
	namespace rcp = MemMngPack;
	namespace rmp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	THROW_EXCEPTION(
		space_x.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_Gc(...): Error!" );
	THROW_EXCEPTION(
		space_c.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_Gc(...): Error!" );
	THROW_EXCEPTION(
		C.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_Gc(...): Error!" );
	THROW_EXCEPTION(
		N.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_Gc(...): Error!" );
	THROW_EXCEPTION(
		Gc == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize_Gc(...): Error!" );
#endif

	const size_type
		n = space_x->dim(),
		m = space_c->dim();

	MatrixCompositeStd
		&Gc_comp = dyn_cast<MatrixCompositeStd>(*Gc);
	
	//
	// Gc = [ C'; N' ]
	//

	Gc_comp.reinitialize(n,m);
	// Add the C matrix object
	typedef rcp::ref_count_ptr<rmp::ReleaseResource_ref_count_ptr<MatrixWithOpNonsingular> > C_rr_ptr_ptr_t;
	C_rr_ptr_ptr_t
		C_rr_ptr_ptr = rcp::rcp(new rmp::ReleaseResource_ref_count_ptr<MatrixWithOpNonsingular>(C));
	Gc_comp.add_matrix(
		var_dep.lbound()-1, 0    // row_offset, col_offset
		,1.0                     // alpha
		,C_rr_ptr_ptr->ptr.get() // A
		,C_rr_ptr_ptr            // A_release
		,BLAS_Cpp::trans         // A_trans
		);
	// Add the N matrix object
	typedef rcp::ref_count_ptr<rmp::ReleaseResource_ref_count_ptr<MatrixWithOp> > N_rr_ptr_ptr_t;
	N_rr_ptr_ptr_t
		N_rr_ptr_ptr = rcp::rcp(new rmp::ReleaseResource_ref_count_ptr<MatrixWithOp>(N));
	Gc_comp.add_matrix(
		var_indep.lbound()-1, 0  // row_offset, col_offset
		,1.0                     // alpha
		,N_rr_ptr_ptr->ptr.get() // A
		,N_rr_ptr_ptr            // A_release
		,BLAS_Cpp::trans         // A_trans
		);
	// Finish construction
	Gc_comp.finish_construction( space_x, space_c );
}

void BasisSystemCompositeStd::get_C_N(
	MatrixWithOp               *Gc
	,MatrixWithOpNonsingular   **C
	,MatrixWithOp              **N
	)
{
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	THROW_EXCEPTION(
		Gc == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::get_C_N(...): Error!" );
	THROW_EXCEPTION(
		C == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::get_C_N(...): Error!" );
	THROW_EXCEPTION(
		N == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::get_C_N(...): Error!" );
#endif
	
	// Get reference to concrete Gc matrix subclass
	MatrixCompositeStd
		&Gc_comp = dyn_cast<MatrixCompositeStd>(*Gc);
	// Get referencs to the aggregate C and N matrtices
	MatrixCompositeStd::matrix_list_t::const_iterator
		mat_itr = Gc_comp.matrices_begin(),
		mat_end = Gc_comp.matrices_end();
	if( mat_itr != mat_end ) {
		assert(mat_itr != mat_end);
		*C = &dyn_cast<MatrixWithOpNonsingular>(
			const_cast<MatrixWithOp&>(*(mat_itr++)->A_) );
		assert(mat_itr != mat_end);
		*N = &const_cast<MatrixWithOp&>(*(mat_itr++)->A_);
		assert(mat_itr == mat_end);
	}
	else {
		*C = NULL;
		*N = NULL;
	}
}

void BasisSystemCompositeStd::get_C_N(
	const MatrixWithOp               &Gc
	,const MatrixWithOpNonsingular   **C
	,const MatrixWithOp              **N
	)
{
	using DynamicCastHelperPack::dyn_cast;
#ifdef _DEBUG
	THROW_EXCEPTION(
		C == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::get_C_N(...): Error!" );
	THROW_EXCEPTION(
		N == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::get_C_N(...): Error!" );
#endif
	
	// Get reference to concrete Gc matrix subclass
	const AbstractLinAlgPack::MatrixCompositeStd
		&Gc_comp = dyn_cast<const AbstractLinAlgPack::MatrixCompositeStd>(Gc);
	// Get referencs to the aggregate C and N matrtices
	MatrixCompositeStd::matrix_list_t::const_iterator
		mat_itr = Gc_comp.matrices_begin(),
		mat_end = Gc_comp.matrices_end();
	if( mat_itr != mat_end ) {
		assert(mat_itr != mat_end);
		*C = &dyn_cast<const MatrixWithOpNonsingular>(*(mat_itr++)->A_);
		assert(mat_itr != mat_end);
		*N = &dyn_cast<const MatrixWithOp>(*(mat_itr++)->A_);
		assert(mat_itr == mat_end);
	}
	else {
		THROW_EXCEPTION(
			true, std::invalid_argument
			,"BasisSystemCompositeStd::get_C_N(...): Error, "
			"The Gc matrix object has not been initialized with C and N!" );
	}
}

// Constructors / initializers

BasisSystemCompositeStd::BasisSystemCompositeStd()
	:BasisSystem(MemMngPack::null,MemMngPack::null)
{}

BasisSystemCompositeStd::BasisSystemCompositeStd(
	const VectorSpace::space_ptr_t       &space_x
	,const VectorSpace::space_ptr_t      &space_c
	,const mat_nonsing_fcty_ptr_t        &factory_C
	,const mat_sym_fcty_ptr_t            &factory_transDtD
	,const mat_sym_nonsing_fcty_ptr_t    &factory_S
	)
	:BasisSystem(MemMngPack::null,MemMngPack::null)
{
	namespace mmp = MemMngPack;
	this->initialize(
		space_x,Range1D(1,space_c->dim()),Range1D(space_c->dim()+1,space_x->dim())
		,space_c,factory_C,factory_transDtD,factory_S
		);
}

BasisSystemCompositeStd::BasisSystemCompositeStd(
	const VectorSpace::space_ptr_t       &space_x
	,const Range1D                       &var_dep
	,const Range1D                       &var_indep
	,const VectorSpace::space_ptr_t      &space_c
	,const mat_nonsing_fcty_ptr_t        &factory_C
	,const mat_sym_fcty_ptr_t            &factory_transDtD
	,const mat_sym_nonsing_fcty_ptr_t    &factory_S
	,const mat_fcty_ptr_t                &factory_D
	,const VectorSpace::space_ptr_t      &space_h
	,const mat_fcty_ptr_t                &factory_GhUP
	)
	:BasisSystem(MemMngPack::null,MemMngPack::null)
{
	this->initialize(
		space_x,var_dep,var_indep,space_c,factory_C,factory_transDtD,factory_S
		,factory_D,space_h,factory_GhUP
		);
}
	
void BasisSystemCompositeStd::initialize(
	const VectorSpace::space_ptr_t       &space_x
	,const Range1D                       &var_dep
	,const Range1D                       &var_indep
	,const VectorSpace::space_ptr_t      &space_c
	,const mat_nonsing_fcty_ptr_t        &factory_C
	,const mat_sym_fcty_ptr_t            &factory_transDtD
	,const mat_sym_nonsing_fcty_ptr_t    &factory_S
	,const mat_fcty_ptr_t                &factory_D
	,const VectorSpace::space_ptr_t      &space_h
	,const mat_fcty_ptr_t                &factory_GhUP
	)
{
	namespace rcp = MemMngPack;
	namespace afp = MemMngPack;

	THROW_EXCEPTION(
		space_x.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize(...): Error!" );
	THROW_EXCEPTION(
		var_dep.size() + var_indep.size() != space_x->dim(), std::invalid_argument
		,"BasisSystemCompositeStd::initialize(...): Error!" );
	THROW_EXCEPTION(
		var_dep.lbound() < var_indep.lbound() && (var_dep.lbound() != 1 || var_dep.ubound()+1 != var_indep.lbound())
		, std::invalid_argument
		,"BasisSystemCompositeStd::initialize(...): Error!" );
	THROW_EXCEPTION(
		var_dep.lbound() >= var_indep.lbound() && (var_indep.lbound() != 1 || var_indep.ubound()+1 != var_dep.lbound())
		, std::invalid_argument
		,"BasisSystemCompositeStd::initialize(...): Error!" );
	THROW_EXCEPTION(
		factory_C.get() == NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize(...): Error!" );
	THROW_EXCEPTION(
		space_h.get() == NULL && factory_GhUP.get() != NULL, std::invalid_argument
		,"BasisSystemCompositeStd::initialize(...): Error!" );

	space_x_         = space_x;
	var_dep_         = var_dep;
	var_indep_       = var_indep;
	space_c_         = space_c;
	factory_C_       = factory_C;
	factory_D_       = factory_D;
	space_h_         = space_h;
	inequ_undecomp_  = ( space_h.get() != NULL ? Range1D(1,space_h->dim()) : Range1D::Invalid );
	factory_GhUP_    = factory_GhUP;
	if( factory_D_.get() == NULL ) {
		factory_D_ = afp::abstract_factory_std_alloc<MatrixWithOp,MultiVectorMutable>(
			AllocatorMultiVectorMutable(space_x_->sub_space(var_dep),var_indep.size() ) );
	}
	if( space_h_.get() != NULL && factory_GhUP_.get() == NULL ) {
		factory_GhUP_ = afp::abstract_factory_std_alloc<MatrixWithOp,MultiVectorMutable>(
			AllocatorMultiVectorMutable(space_h_,var_indep.size() ) );
	}
	BasisSystem::initialize(factory_transDtD,factory_S);
}

void BasisSystemCompositeStd::set_uninitialized()
{
	namespace rcp = MemMngPack;

	space_x_         = rcp::null;
	var_dep_         = Range1D::Invalid;
	var_indep_       = Range1D::Invalid;
	factory_C_       = rcp::null;
	factory_D_       = rcp::null;
	space_h_         = rcp::null;
	inequ_undecomp_  = Range1D::Invalid;
	factory_GhUP_    = rcp::null;
}

const VectorSpace::space_ptr_t&
BasisSystemCompositeStd::space_x() const
{
	return space_x_;
}

const VectorSpace::space_ptr_t&
BasisSystemCompositeStd::space_c() const
{
	return space_c_;
}

const VectorSpace::space_ptr_t&
BasisSystemCompositeStd::space_h() const
{
	return space_h_;
}

// To be overridden by subclasses

void BasisSystemCompositeStd::update_D(
	const MatrixWithOpNonsingular&  C
	,const MatrixWithOp&            N
	,MatrixWithOp*                  D
	,EMatRelations                  mat_rel
	) const
{
	using LinAlgOpPack::M_StInvMtM;
	M_StInvMtM( D, -1.0, C, BLAS_Cpp::no_trans, N, BLAS_Cpp::no_trans ); // D = -inv(C)*N
}

void BasisSystemCompositeStd::update_GhUP(
	const MatrixWithOpNonsingular&  C
	,const MatrixWithOp&            N
	,const MatrixWithOp*            D
	,MatrixWithOp*                  GhUP
	,EMatRelations                  mat_rel
	) const
{
	assert(0); // ToDo: Implement when needed!
}

// Overridden from BasisSytem

const BasisSystem::mat_nonsing_fcty_ptr_t
BasisSystemCompositeStd::factory_C() const
{
	return factory_C_;
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemCompositeStd::factory_D() const
{
	return factory_D_;
}

const BasisSystem::mat_fcty_ptr_t
BasisSystemCompositeStd::factory_GhUP() const
{
	return factory_GhUP_;
}

Range1D BasisSystemCompositeStd::var_dep() const
{
	return var_dep_;
}

Range1D BasisSystemCompositeStd::var_indep() const
{
	return var_indep_;
}

Range1D BasisSystemCompositeStd::inequ_undecomp() const
{
	return inequ_undecomp_;
}

void BasisSystemCompositeStd::update_basis(
	const MatrixWithOp*         Gc
	,const MatrixWithOp*        Gh
	,MatrixWithOpNonsingular*   C
	,MatrixWithOp*              D
	,MatrixWithOp*              GcUP
	,MatrixWithOp*              GhUP
	,EMatRelations              mat_rel
	,std::ostream               *out
	) const
{
	using DynamicCastHelperPack::dyn_cast;
	const index_type
		n  = var_dep_.size() + var_indep_.size(),
		m  = var_dep_.size(),
		mI = inequ_undecomp_.size();
	THROW_EXCEPTION(
		n == 0, std::logic_error
		,"BasisSystemCompositeStd::update_basis(...): Error, this must be initialized first!" );
	THROW_EXCEPTION(
		Gc == NULL, std::logic_error
		,"BasisSystemCompositeStd::update_basis(...): Error, Gc can not be NULL!" );
	THROW_EXCEPTION(
		mI == 0 && Gh != NULL, std::logic_error
		,"BasisSystemCompositeStd::update_basis(...): Error, Gh must be NULL!" );
	THROW_EXCEPTION(
		GcUP, std::logic_error
		,"BasisSystemCompositeStd::update_basis(...): Error, GcUP must be NULL!" );
	THROW_EXCEPTION(
		mI == 0 && GhUP != NULL, std::logic_error
		,"BasisSystemCompositeStd::update_basis(...): Error, GhUP must be NULL!" );
	THROW_EXCEPTION(
		C == NULL && D == NULL, std::logic_error
		,"BasisSystemCompositeStd::update_basis(...): Error, C or D must be non-NULL!" );
	// Get references to the aggregate C and N matrices
	const MatrixWithOpNonsingular
		*C_aggr = NULL;
	const MatrixWithOp
		*N_aggr = NULL;
	get_C_N( *Gc, &C_aggr, &N_aggr );
	// Setup C
	if( C ) {
		*C = *C_aggr;
	}
	// Compute D
	if( D ) {
		update_D(*C_aggr,*N_aggr,D,mat_rel);
	}
	// Compute GhUP
	if( GhUP ) {
		update_GhUP(*C_aggr,*N_aggr,D,GhUP,mat_rel);
	}
}

} // end namespace AbstractLinAlgPack
