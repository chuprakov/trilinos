// ///////////////////////////////////////////////////////////////////////
// NLPSerialPreprocessExplJac.cpp
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

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include <typeinfo>
#include <algorithm>

#include "NLPInterfacePack/include/NLPSerialPreprocessExplJac.h"
#include "SparseLinAlgPack/include/MatrixSparseCOORSerial.h"
#include "SparseLinAlgPack/include/PermutationSerial.h"
#include "SparseLinAlgPack/include/VectorDenseEncap.h"
#include "AbstractLinAlgPack/include/MatrixPermAggr.h"
#include "AbstractLinAlgPack/include/MatrixCompositeStd.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/IVector.h"
#include "LinAlgPack/include/PermVecMat.h"
#include "ThrowException.h"
#include "dynamic_cast_verbose.h"
#include "AbstractFactoryStd.h"

namespace NLPInterfacePack {

// NLPSerialPreprocessExplJac

// Constructors / initializers

NLPSerialPreprocessExplJac::NLPSerialPreprocessExplJac(
	const factory_mat_ptr_t     &factory_Gc_orig
	,const factory_mat_ptr_t    &factory_Gh_orig
	)
	:initialized_(false)
{
	this->set_mat_factories(factory_Gc_orig,factory_Gh_orig);
}

void NLPSerialPreprocessExplJac::set_mat_factories(
	const factory_mat_ptr_t     &factory_Gc_orig
	,const factory_mat_ptr_t    &factory_Gh_orig
	)
{
	namespace rcp = MemMngPack;
	namespace afp = MemMngPack;
	if(factory_Gc_orig.get())
		factory_Gc_orig_ = factory_Gc_orig;
	else 
		factory_Gc_orig_ = rcp::rcp(
			new afp::AbstractFactoryStd<MatrixWithOp,MatrixSparseCOORSerial>() );
	if(factory_Gh_orig.get())
		factory_Gh_orig_ = factory_Gh_orig;
	else 
		factory_Gc_orig_ = rcp::rcp(
			new afp::AbstractFactoryStd<MatrixWithOp,MatrixSparseCOORSerial>() );
	factory_Gc_ = rcp::rcp( new afp::AbstractFactoryStd<MatrixWithOp,MatrixPermAggr>() );
	factory_Gh_ = rcp::rcp( new afp::AbstractFactoryStd<MatrixWithOp,MatrixPermAggr>() );
}

// Overridden public members from NLP

void NLPSerialPreprocessExplJac::initialize()
{
	namespace rcp = MemMngPack;

	if( initialized_  && !imp_nlp_has_changed() ) {
		// The subclass NLP has not changed so we can just
		// slip this preprocessing.
		NLPFirstOrderInfo::initialize();
		NLPSerialPreprocess::initialize();  // Some duplication but who cares!
		return;
	}

	// Initialize the base object first
	NLPFirstOrderInfo::initialize();
	NLPSerialPreprocess::initialize();  // Some duplication but who cares!

	const NLP::vec_space_ptr_t
		space_x = this->space_x(),
		space_c = this->space_c(),
		space_h = this->space_h();

	// Initialize the storage for the intermediate quanities
	Gc_nz_full_ = imp_Gc_nz_full();			// Get the estimated number of nonzeros in Gc
	Gc_val_full_.resize(Gc_nz_full_);
	Gc_ivect_full_.resize(Gc_nz_full_);
	Gc_jvect_full_.resize(Gc_nz_full_);
	Gc_perm_new_basis_updated_ = false;
	Gh_nz_full_ = imp_Gh_nz_full();			// Get the estimated number of nonzeros in Gh
	Gh_val_full_.resize(Gh_nz_full_);
	Gh_ivect_full_.resize(Gh_nz_full_);
	Gh_jvect_full_.resize(Gh_nz_full_);
	Gh_perm_new_basis_updated_ = false;

	// If you get here then the initialization went Ok.
	initialized_ = true;
}

bool NLPSerialPreprocessExplJac::is_initialized() const {
	return initialized_;
}

// Overridden public members from NLPFirstOrderInfo

const NLPFirstOrderInfo::mat_fcty_ptr_t
NLPSerialPreprocessExplJac::factory_Gc() const
{
	return factory_Gc_;
}

const NLPFirstOrderInfo::mat_fcty_ptr_t
NLPSerialPreprocessExplJac::factory_Gh() const
{
	return factory_Gh_;
}

void NLPSerialPreprocessExplJac::set_Gc(MatrixWithOp* Gc)
{
	using DynamicCastHelperPack::dyn_cast;
	assert_initialized();
	if( Gc != NULL ) {
		dyn_cast<MatrixPermAggr>(*Gc); // With throw exception if not correct type!
	}
	NLPFirstOrderInfo::set_Gc(Gc);
}

void NLPSerialPreprocessExplJac::set_Gh(MatrixWithOp* Gh)
{
	using DynamicCastHelperPack::dyn_cast;
	assert_initialized();
	if( Gh != NULL ) {
		dyn_cast<MatrixPermAggr>(*Gh); // With throw exception if not correct type!
	}
	NLPFirstOrderInfo::set_Gh(Gh);
}

// Overridden public members from NLPVarReductPerm

bool NLPSerialPreprocessExplJac::get_next_basis(
	Permutation*  P_var,   Range1D* var_dep
	,Permutation* P_equ,   Range1D* equ_decomp
	,Permutation* P_inequ, Range1D* inequ_decomp
	)
{
	const bool new_basis = NLPSerialPreprocess::get_next_basis(
		P_var,var_dep,P_equ,equ_decomp,P_inequ,inequ_decomp
		);
	if( new_basis ) {
		Gc_perm_new_basis_updated_ = false;
		Gh_perm_new_basis_updated_ = false;
	}
	return new_basis;
}

void NLPSerialPreprocessExplJac::set_basis(
	const Permutation   &P_var,   const Range1D  &var_dep
	,const Permutation  *P_equ,   const Range1D  *equ_decomp
	,const Permutation  *P_inequ, const Range1D  *inequ_decomp
	)
{
	NLPSerialPreprocess::set_basis(
		P_var,var_dep,P_equ,equ_decomp,P_inequ,inequ_decomp
		);
	Gc_perm_new_basis_updated_ = false;
	Gh_perm_new_basis_updated_ = false;
}

// Overridden protected members from NLPFirstOrderInfo

void NLPSerialPreprocessExplJac::imp_calc_Gc(
	const VectorWithOp& x, bool newx
	,const FirstOrderInfo& first_order_info
	) const
{
	imp_calc_Gc_or_Gh(
		true        // Calcuate Gc
		,x,newx
		,first_order_info
		);
}

void NLPSerialPreprocessExplJac::imp_calc_Gh(
	const VectorWithOp& x, bool newx
	,const FirstOrderInfo& first_order_info
	) const
{
	imp_calc_Gc_or_Gh(
		false        // Calcuate Gh
		,x,newx
		,first_order_info
		);
}

// protected members

void NLPSerialPreprocessExplJac::assert_initialized() const
{
	THROW_EXCEPTION(
		!initialized_, UnInitialized
		,"NLPSerialPreprocessExplJac : The nlp has not been initialized yet" );
}

// private members

void NLPSerialPreprocessExplJac::imp_calc_Gc_or_Gh(
	bool calc_Gc
	,const VectorWithOp& x, bool newx
	,const FirstOrderInfo& first_order_info
	) const
{
	namespace rcp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;

	assert_initialized();

	const Range1D
		var_dep      = this->var_dep(),
		equ_decomp   = this->equ_decomp(),
		inequ_decomp = this->equ_decomp();
	const size_type
		n       = this->n(),
		n_full  = this->imp_n_full(),
		m_full  = this->imp_m_full(),
		mI_full = this->imp_mI_full(),
		num_cols = (calc_Gc ? m_full : mI_full);


	//
	// Get references to the constituent objects
	//

	// Get the concrete type for the Jacobian matrix
	MatrixPermAggr
		&G_aggr = dyn_cast<MatrixPermAggr>( calc_Gc
											? *first_order_info.Gc
											: *first_order_info.Gh );
	// Get smart pointers to the constituent members
	rcp::ref_count_ptr<MatrixWithOp>
		G_orig = rcp::rcp_const_cast<MatrixWithOp>( G_aggr.mat_orig() );
	rcp::ref_count_ptr<PermutationSerial>
		P_row = rcp::rcp_dynamic_cast<PermutationSerial>(
			rcp::rcp_const_cast<Permutation>( G_aggr.row_perm() ) );  // variable permutation
	rcp::ref_count_ptr<PermutationSerial>
		P_col = rcp::rcp_dynamic_cast<PermutationSerial>(
			rcp::rcp_const_cast<Permutation>( G_aggr.col_perm() ) );  // constraint permutation
	rcp::ref_count_ptr<const MatrixWithOp>
		G_perm = G_aggr.mat_perm();
	// Remove references to G_orig, G_perm, P_row and P_col.
	G_aggr.set_uninitialized();
	// Allocate the original matrix object if not done so yet
	if( G_orig.get() == NULL || G_orig.count() > 1 )
		G_orig = (calc_Gc ? factory_Gc_orig_ : factory_Gh_orig_)->create();
	// Get reference to the MatrixLoadSparseElements interface
	MatrixLoadSparseElements
		&G_lse = dyn_cast<MatrixLoadSparseElements>(*G_orig);

	//
	// Calcuate the full explicit Jacobian
	//

	set_x_full( VectorDenseEncap(x)(), newx, &x_full() );
	if(calc_Gc)
		imp_calc_Gc_full( x_full(), newx, first_order_expl_info() );
	else
		imp_calc_Gh_full( x_full(), newx, first_order_expl_info() );

	// Now get the actual number of nonzeros
	const size_type nz_full = (calc_Gc ? Gc_nz_full_ : Gh_nz_full_);

	// Determine if we need to set the structure and the nonzeros or just the nonzeros
	const bool load_struct = (G_lse.nz() == 0);

	size_type G_nz_previous;
	if( load_struct ) {
		G_lse.reinitialize(n,num_cols,nz_full); // The actual number of nonzeros will be minus the fixed variables
	}
	else {
		G_nz_previous = G_lse.nz();
		G_lse.reset_to_load_values();           // Use row and column indexes already set (better be same insert order!)
	}
		
	// Get pointers to the full COOR matrix just updated
	value_type		*val_full		= &Gc_val_full_[0],
					*val_full_end	= val_full + nz_full;
 	index_type		*ivect_full		= &Gc_ivect_full_[0],
					*jvect_full		= &Gc_jvect_full_[0];

	// Get pointers to buffers to fill with nonzero entries
	value_type			*val    = NULL;
 	index_type			*ivect  = NULL,
		                *jvect  = NULL;
	G_lse.get_load_nonzeros_buffers(
		nz_full // We may actually load less
		,&val
		,load_struct ? &ivect : NULL
		,load_struct ? &jvect : NULL
		);

	// Remove the variables fixed by bounds and adjust the row indices accordingly
	// Note: only update the row and column indices the first time since these
	// will never change and are independent of basis choose.

	value_type  *val_itr    = val;
 	index_type  *ivect_itr  = ivect,
		        *jvect_itr  = jvect,
		        nz          = 0;
	const IVector& var_full_to_remove_fixed = this->var_full_to_remove_fixed();
	if( load_struct ) {
		// Fill values and i and j
		for( ; val_full != val_full_end ; ++val_full, ++ivect_full, ++jvect_full) {
#ifdef _DEBUG
			assert( 0 <= *ivect_full && *ivect_full <= n_full );
#endif
			size_type var_idx = var_full_to_remove_fixed(*ivect_full);
#ifdef _DEBUG
			assert( 0 < var_idx && var_idx <= n_full );
#endif
			if(var_idx <= n) {
				// This is not a fixed variable
				*val_itr++ = *val_full;
				// Also fill the row and column indices
				*ivect_itr++ = var_idx;
				*jvect_itr++ = *jvect_full;
				++nz;
			}
		}
	}
	else {
		// Just fill values
		for( ; val_full != val_full_end ; ++val_full, ++ivect_full) {
#ifdef _DEBUG
			assert( 0 <= *ivect_full && *ivect_full <= n_full );
#endif
			size_type var_idx = var_full_to_remove_fixed(*ivect_full);
#ifdef _DEBUG
			assert( 0 < var_idx && var_idx <= n_full );
#endif
			if(var_idx <= n) {
				// This is not a fixed variable
				*val_itr++ = *val_full;
				++nz;
			}
		}
		// Check that the number of nonzeros added matches the number of nonzeros in G
		THROW_EXCEPTION(
			G_nz_previous != nz, std::runtime_error
			,"NLPSerialPreprocessExplJac::imp_calc_Gc_or_Gh(...): Error, "
			"The number of added nonzeros does not match the number of nonzeros "
			"in the previous matrix load!." );
	}

	// Commit the nonzeros
	G_lse.commit_load_nonzeros_buffers(
		nz  // The actual number of nonzeros to set
		,&val
		,load_struct ? &ivect : NULL
		,load_struct ? &jvect : NULL
		);
	G_lse.finish_construction();

	//
	// Setup permuted view
	//

	// Setup row (variable) permutation
	if( P_row.get() == NULL || P_col.count() > 1 )
	    P_row = rcp::rcp(new PermutationSerial());
	rcp::ref_count_ptr<IVector>        var_perm;
	if( P_row->perm().get() == NULL )  var_perm = rcp::rcp(new IVector(n_full));
	else                               var_perm = rcp::rcp_const_cast<IVector>(P_row->perm());
	*var_perm = this->var_perm();
	P_row->initialize(var_perm,rcp::null);
	// Setup column (constraint) permutation
	if( P_col.get() == NULL || P_col.count() > 1 )
	    P_col = rcp::rcp(new PermutationSerial());
	rcp::ref_count_ptr<IVector>        con_perm;
	if( P_col->perm().get() == NULL )  con_perm = rcp::rcp(new IVector(calc_Gc?m_full:mI_full));
	else                               con_perm = rcp::rcp_const_cast<IVector>(P_col->perm());
	if( calc_Gc )                      *con_perm = this->equ_perm();
	else                               LinAlgPack::identity_perm(con_perm.get());
	P_col->initialize(con_perm,rcp::null);
	// Setup G_perm
	int num_row_part, num_col_part;
	index_type row_part[3], col_part[3];
	if(var_dep.size()) {
		num_row_part = 2;
		row_part[0] = 1;
		row_part[1] = (var_dep.lbound() == 1 ? var_dep.ubound()+1 : var_dep.lbound());
		row_part[2] = n+1;
	}
	else {
		num_row_part = 1;
		row_part[0] = 1;
		row_part[1] = n+1;
	}
	if(calc_Gc) {
		if(equ_decomp.size()) {
			num_col_part = 2;
			col_part[0] = 1;
			col_part[1] = (equ_decomp.lbound() == 1 ? equ_decomp.ubound()+1 : equ_decomp.lbound());
			col_part[2] = m_full+1;
		}
		else {
			num_col_part = 1;
			col_part[0] = 1;
			col_part[1] = m_full+1;
		}
	}
	else {
		if(inequ_decomp.size()) {
			num_col_part = 2;
			col_part[0] = 1;
			col_part[1] = (inequ_decomp.lbound() == 1 ? inequ_decomp.ubound()+1 : inequ_decomp.lbound());
			col_part[2] = mI_full+1;
		}
		else {
			num_col_part = 1;
			col_part[0] = 1;
			col_part[1] = mI_full+1;
		}
	}
	if( G_perm.get() == NULL || (calc_Gc ? !Gc_perm_new_basis_updated_ : !Gh_perm_new_basis_updated_) ) {
		G_perm = G_orig->perm_view(
			P_row.get(),row_part,num_row_part
			,P_col.get(),col_part,num_col_part
			);
	}
	else {
		G_perm = G_orig->perm_view_update(
			P_row.get(),row_part,num_row_part
			,P_col.get(),col_part,num_col_part
			,G_perm
			);
	}
	(calc_Gc ? Gc_perm_new_basis_updated_ : Gh_perm_new_basis_updated_) = true;

	//
	// Reinitialize the aggregate matrix object.
	//

	G_aggr.initialize(G_orig,P_row,P_col,G_perm);

}

} // end namespace NLPInterfacePack
