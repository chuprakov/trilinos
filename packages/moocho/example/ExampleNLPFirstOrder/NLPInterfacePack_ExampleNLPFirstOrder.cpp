// //////////////////////////////////////////
// ExampleNLPFirstOrderInfo.cpp
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

#include <stdexcept>

#include "ExampleNLPFirstOrderInfo.h"
#include "AbstractLinAlgPack/include/MatrixSpaceStd.h"
#include "AbstractLinAlgPack/include/MatrixSymDiagonalStd.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "Range1D.h"
#include "ReleaseResource_ref_count_ptr.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace NLPInterfacePack {

ExampleNLPFirstOrderInfo::ExampleNLPFirstOrderInfo(
	const VectorSpace::space_ptr_t&  vec_space
	,value_type                      xo
	,bool                            has_bounds
	,bool                            dep_bounded
	)
	:ExampleNLPFirstOrderDirect(vec_space,xo,has_bounds,dep_bounded)
	,initialized_(false)
{}

// Overridden public members from NLPFirstOrderInfo

void ExampleNLPFirstOrderInfo::set_Gc(MatrixWithOp* Gc)
{
	if(Gc) // Throw an exception if this matrix is the wrong type!
		DynamicCastHelperPack::dyn_cast<AbstractLinAlgPack::MatrixCompositeStd>(*Gc);
	NLPFirstOrderInfo::set_Gc(Gc);
}

const NLPFirstOrderInfo::mat_space_ptr_t&
ExampleNLPFirstOrderInfo::space_Gc() const
{
	return space_Gc_;
}

const NLPFirstOrderInfo::mat_space_ptr_t&
ExampleNLPFirstOrderInfo::space_Gh() const
{
	return ExampleNLPFirstOrderDirect::space_Gh();
}

// Overridden public members from NLP

void ExampleNLPFirstOrderInfo::initialize()
{
	namespace rcp = ReferenceCountingPack;

	ExampleNLPFirstOrderDirect::initialize();
	NLPFirstOrderInfo::initialize();

	space_Gc_ = rcp::rcp_implicit_cast<NLPFirstOrderInfo::mat_space_ptr_t::element_type>(
		rcp::ref_count_ptr<MatrixSpaceStd<MatrixWithOp,MatrixCompositeStd> >(
			new MatrixSpaceStd<MatrixWithOp,MatrixCompositeStd>(
				this->space_c(), this->space_x()
				) ) );
	
	initialized_ = true;
}

bool ExampleNLPFirstOrderInfo::is_initialized() const
{
	return initialized_;
}

// Overridden protected members from NLPFirstOrderInfo

void ExampleNLPFirstOrderInfo::imp_calc_Gc(
	const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const
{
	namespace rcp = ReferenceCountingPack;
	namespace rmp = ResourceManagementPack;
	using DynamicCastHelperPack::dyn_cast;

	const index_type
		n = this->n(),
		m = this->m();

	assert(first_order_info.Gc);  // Should not be NULL if this member is called but check to be sure!
	// Get reference to concrete Gc matrix subclass
	AbstractLinAlgPack::MatrixCompositeStd
		&Gc_comp = DynamicCastHelperPack::dyn_cast<AbstractLinAlgPack::MatrixCompositeStd>(
			*first_order_info.Gc );
	//
	// Initialize Gc with C and N matrices this has not already been done.
	//
	if( Gc_comp.rows() != n || Gc_comp.cols() != m || Gc_comp.matrices_begin() == Gc_comp.matrices_end() ) {
		Gc_comp.reinitalize(n,m);
		// Gc = [ C'; N' ]
		VectorSpace::space_ptr_t
			space_x_DI = this->space_c();
		index_type row_offset = 0;
		for(int k = 0; k < 2; ++k, row_offset += this->m() ) { // k == 0 -> add C',  k == 1 -> add N'
			typedef rcp::ref_count_ptr<rmp::ReleaseResource_ref_count_ptr<MatrixSymDiagonalStd> >
				rr_ptr_ptr_t;
			// this is a mess of a data structure but it is correct
			rr_ptr_ptr_t
				rr_ptr_ptr = new rmp::ReleaseResource_ref_count_ptr<MatrixSymDiagonalStd>(
					new MatrixSymDiagonalStd( space_x_DI->create_member() ) );
			// Add the matrix object
			Gc_comp.add_matrix(
				row_offset, 0         // row_offset, col_offset
				,1.0                  // alpha
				,rr_ptr_ptr->ptr.get() // A
				,rcp::rcp_implicit_cast<MatrixCompositeStd::release_resource_ptr_t::element_type>(
					rr_ptr_ptr )      // A_release
				,BLAS_Cpp::trans      // A_trans
				);
		}
		// Create the composite vector spaces	
		Gc_comp.finish_construction( this->space_c(), this->space_x() );
	}
	//
	// Set the diagonals of the C and N matrices for the current point.
	//
	assert( Gc_comp.matrices_begin() != Gc_comp.matrices_end() );
	// Get references to C and N
	MatrixCompositeStd::matrix_list_t::iterator
		mat_itr = Gc_comp.matrices_begin(),
		mat_end = Gc_comp.matrices_end();
	MatrixSymDiagonalStd
		&C = dyn_cast<MatrixSymDiagonalStd>(*const_cast<MatrixWithOp*>(mat_itr->A_));
	assert(mat_itr != mat_end);
	MatrixSymDiagonalStd
		&N = dyn_cast<MatrixSymDiagonalStd>(*const_cast<MatrixWithOp*>((++mat_itr)->A_));
	// Get x = [ x_D' x_I ]
	VectorWithOp::vec_ptr_t
		x_D = x.sub_view(this->var_dep()),
		x_I = x.sub_view(this->var_indep());
	// Set the diagonals of C and N (this is the only computation going on here)
	C.diag() = *x_I;          // C.diag = x_I - 1.0
	Vp_S( &C.diag(),  -1.0 ); // ...
	N.diag() = *x_D;          // N.diag = x_D - 10.0
	Vp_S( &N.diag(), -10.0 ); // ...

}

void ExampleNLPFirstOrderInfo::imp_calc_Gh(
	const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const
{
	assert(0); // This should never be called!
}

}	// end namespace NLPInterfacePack
