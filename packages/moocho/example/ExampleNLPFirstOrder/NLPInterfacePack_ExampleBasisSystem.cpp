// ///////////////////////////////////////////////////////////////////
// ExampleBasisSystem.cpp
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

#include "ExampleBasisSystem.h"
#include "AbstractLinAlgPack/include/MatrixSpaceStd.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "RTOpPack/include/RTOpCppC.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace NLPInterfacePack {
 
ExampleBasisSystem::ExampleBasisSystem( const VectorSpace::space_ptr_t& space_x_DI )
{
	this->initialize(space_x_DI);
}
	
void ExampleBasisSystem::initialize( const VectorSpace::space_ptr_t& space_x_DI )
{
	namespace rcp = ReferenceCountingPack;
	space_x_DI_ = space_x_DI;
	if( space_x_DI_.get() ) {
		const size_type
			m = space_x_DI_->dim(),
			n = 2*m;
		assert( n%2 == 0 );
		var_dep_   = Range1D(1,m);
		var_indep_ = Range1D(m+1,n);
		space_C_ = rcp::rcp_implicit_cast<BasisSystem::mat_nonsing_space_ptr_t::element_type>(
			rcp::ref_count_ptr<MatrixSpaceStd<MatrixWithOpNonsingular,MatrixSymDiagonalStd> >(
				new MatrixSpaceStd<MatrixWithOpNonsingular,MatrixSymDiagonalStd>(
					space_x_DI, space_x_DI
					) ) );
		space_D_ = rcp::rcp_implicit_cast<BasisSystem::mat_space_ptr_t::element_type>(
			rcp::ref_count_ptr<MatrixSpaceStd<MatrixWithOp,MatrixSymDiagonalStd> >(
				new MatrixSpaceStd<MatrixWithOp,MatrixSymDiagonalStd>(
					space_x_DI, space_x_DI
					) ) );
	}
	else {
		var_dep_    = Range1D::Invalid;
		var_indep_  = Range1D::Invalid;
		space_C_    = NULL;
		space_D_    = NULL;
	}
}

const BasisSystem::mat_nonsing_space_ptr_t&
ExampleBasisSystem::space_C() const
{
	return space_C_;
}

const BasisSystem::mat_space_ptr_t&
ExampleBasisSystem::space_D() const
{
	return space_D_;
}

Range1D ExampleBasisSystem::var_dep() const
{
	return var_dep_;
}

Range1D ExampleBasisSystem::var_indep() const
{
	return var_indep_;
}

void ExampleBasisSystem::update_basis(
	const MatrixWithOp*         Gc
	,const MatrixWithOp*        Gh
	,MatrixWithOpNonsingular*   C
	,MatrixWithOp*              D
	,MatrixWithOp*              GcUP
	,MatrixWithOp*              GhUP
	)
{
	using DynamicCastHelperPack::dyn_cast;

	const index_type
		n = var_dep_.size() + var_indep_.size(),
		m = var_dep_.size();

	THROW_EXCEPTION(
		!n, std::logic_error
		,"ExampleBasisSystem::update_basis(...): Error, this must be initialized first!" );
	THROW_EXCEPTION(
		!Gc, std::logic_error
		,"ExampleBasisSystem::update_basis(...): Error, Gc can not be NULL!" );
	THROW_EXCEPTION(
		Gh, std::logic_error
		,"ExampleBasisSystem::update_basis(...): Error, Gh must be NULL!" );
	THROW_EXCEPTION(
		GcUP, std::logic_error
		,"ExampleBasisSystem::update_basis(...): Error, GcUP must be NULL!" );
	THROW_EXCEPTION(
		GhUP, std::logic_error
		,"ExampleBasisSystem::update_basis(...): Error, GhUP must be NULL!" );
	THROW_EXCEPTION(
		!C && !D, std::logic_error
		,"ExampleBasisSystem::update_basis(...): Error, C or D must be non-NULL!" );

	// Get reference to concrete Gc matrix subclass
	const AbstractLinAlgPack::MatrixCompositeStd
		&Gc_comp = dyn_cast<const AbstractLinAlgPack::MatrixCompositeStd>(*Gc);
	// Get referencs to the aggregate C and N matrtices
	MatrixCompositeStd::matrix_list_t::const_iterator
		mat_itr = Gc_comp.matrices_begin(),
		mat_end = Gc_comp.matrices_end();
	THROW_EXCEPTION(
		Gc_comp.rows() != n || Gc_comp.cols() != m || mat_itr == mat_end, std::logic_error
		,"ExampleBasisSystem::update_basis(...): Error, Gc is  not intialized properly with "
		"Gc_comp.rows() = " << Gc_comp.rows() << ", Gc_comp.cols() = " << Gc_comp.cols()
		<< ", n = " << n << " and m = " << m << "!" );
	assert(mat_itr != mat_end);
	const MatrixSymDiagonalStd
		&C_aggr = dyn_cast<const MatrixSymDiagonalStd>(*mat_itr->A_);
	assert(mat_itr != mat_end);
	const MatrixSymDiagonalStd
	 	&N_aggr = dyn_cast<const MatrixSymDiagonalStd>(*(++mat_itr)->A_);
	// Setup C
	if( C ) {
		MatrixSymDiagonalStd
			&C_sym_diag = dyn_cast<MatrixSymDiagonalStd>(*C);
		if( C_sym_diag.rows() != m )
			C_sym_diag.initialize( space_x_DI_->create_member() );
		C_sym_diag.diag() = C_aggr.diag();
	}
	// Compute D
	if( D ) {
		MatrixSymDiagonalStd
		 	&D_sym_diag = dyn_cast<MatrixSymDiagonalStd>(*D);
		if( D_sym_diag.rows() != m )
		 	D_sym_diag.initialize( space_x_DI_->create_member() );
		AbstractLinAlgPack::ele_wise_divide(                          // D_diag = - N_diag ./ C_diag
			-1.0, N_aggr.diag(), C_aggr.diag(), &D_sym_diag.diag() );  // ...
	}
}

} // end namespace NLPInterfacePack
