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
#include "AbstractLinAlgPack/include/BasisSystemCompositeStd.h"
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
	namespace rcp = MemMngPack;

	ExampleNLPFirstOrderDirect::initialize();
	NLPFirstOrderInfo::initialize();

	space_Gc_ = BasisSystemCompositeStd::space_Gc(
		this->space_x(), this->space_c() );
	
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
	namespace rcp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	using AbstractLinAlgPack::Vp_S; // Should not have to do this!

	const index_type
		n = this->n(),
		m = this->m();
	const Range1D
		var_dep   = this->var_dep(),
		var_indep = this->var_indep();

	// Get references to aggregate C and N matrices (if allocated)
	MatrixWithOpNonsingular
		*C_aggr = NULL;
	MatrixWithOp
		*N_aggr = NULL;
	BasisSystemCompositeStd::get_C_N(
		first_order_info.Gc, &C_aggr, &N_aggr ); // Will return NULLs if Gc is not initialized

	// Allocate C and N matrix objects if not done yet!
	rcp::ref_count_ptr<MatrixWithOpNonsingular>
		C_ptr = rcp::null;
	rcp::ref_count_ptr<MatrixWithOp>
		N_ptr = rcp::null;
	if( C_aggr == NULL ) {
		const VectorSpace::space_ptr_t
			space_x  = this->space_x(),
			space_xD = space_x->sub_space(var_dep);
		C_ptr  = rcp::rcp(new MatrixSymDiagonalStd(space_xD->create_member()));
		N_ptr  = rcp::rcp(new MatrixSymDiagonalStd(space_xD->create_member()));
		C_aggr = C_ptr.get();
		N_aggr = N_ptr.get();
	}

	// Get references to concreate C and N matrices
	MatrixSymDiagonalStd
		&C = dyn_cast<MatrixSymDiagonalStd>(*C_aggr);
	MatrixSymDiagonalStd
		&N = dyn_cast<MatrixSymDiagonalStd>(*N_aggr);
	// Get x = [ x_D' x_I ]
	VectorWithOp::vec_ptr_t
		x_D = x.sub_view(var_dep),
		x_I = x.sub_view(var_indep);
	// Set the diagonals of C and N (this is the only computation going on here)
	C.diag() = *x_I;          // C.diag = x_I - 1.0
	Vp_S( &C.diag(),  -1.0 ); // ...
	N.diag() = *x_D;          // N.diag = x_D - 10.0
	Vp_S( &N.diag(), -10.0 ); // ...
	// Initialize the matrix object Gc if not done so yet
	if( C_ptr.get() != NULL ) {
		BasisSystemCompositeStd::initialize_Gc(
			this->space_x(), var_dep, var_indep
			,this->space_c()
			,C_ptr, N_ptr
			,first_order_info.Gc
			);
	}
}

void ExampleNLPFirstOrderInfo::imp_calc_Gh(
	const VectorWithOp& x, bool newx, const FirstOrderInfo& first_order_info) const
{
	assert(0); // This should never be called!
}

}	// end namespace NLPInterfacePack
