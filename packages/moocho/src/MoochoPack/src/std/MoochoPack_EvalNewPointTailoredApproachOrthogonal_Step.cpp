// ////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachOrthogonal_Step.cpp
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

#include "ReducedSpaceSQPPack/src/std/EvalNewPointTailoredApproachOrthogonal_Step.h"
#include "ConstrainedOptimizationPack/src/MatrixIdentConcatStd.h"
#include "NLPInterfacePack/src/NLPFirstOrderDirect.h"
#include "AbstractLinAlgPack/src/MatrixCompositeStd.h"
#include "AbstractLinAlgPack/src/MatrixSymWithOpNonsingular.h"
#include "AbstractLinAlgPack/src/MatrixSymInitDiagonal.h"
#include "AbstractLinAlgPack/src/VectorSpace.h"
#include "AbstractLinAlgPack/src/VectorStdOps.h"
#include "AbstractLinAlgPack/src/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/src/AbstractLinAlgPackAssertOp.h"
#include "AbstractLinAlgPack/src/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace ReducedSpaceSQPPack {

EvalNewPointTailoredApproachOrthogonal_Step::EvalNewPointTailoredApproachOrthogonal_Step(
	const deriv_tester_ptr_t                &deriv_tester
	,const bounds_tester_ptr_t              &bounds_tester
	,EFDDerivTesting                        fd_deriv_testing
	)
	:EvalNewPointTailoredApproach_Step(deriv_tester,bounds_tester,fd_deriv_testing)
{}

// protected

void EvalNewPointTailoredApproachOrthogonal_Step::uninitialize_Y_Uv_Uy(
	MatrixWithOp         *Y
	,MatrixWithOp        *Uy
	,MatrixWithOp        *Vy
	)
{
	using DynamicCastHelperPack::dyn_cast;

	MatrixIdentConcatStd
		*Y_orth = Y ? &dyn_cast<MatrixIdentConcatStd>(*Y)  : NULL;
	MatrixCompositeStd
		*Uy_cpst = Uy ? &dyn_cast<MatrixCompositeStd>(*Uy) : NULL;			
	MatrixCompositeStd
		*Vy_cpst = Vy ? &dyn_cast<MatrixCompositeStd>(*Vy) : NULL;

	if(Y_orth)
		Y_orth->set_uninitialized();
	assert(Uy_cpst == NULL); // ToDo: Implement for undecomposed equalities
	assert(Vy_cpst == NULL); // ToDo: Implement for general inequalities
}

void EvalNewPointTailoredApproachOrthogonal_Step::calc_py_Y_Uy_Vy(
	const NLPFirstOrderDirect   &nlp
	,const D_ptr_t              &D
	,VectorWithOpMutable        *py
	,MatrixWithOp               *Y
	,MatrixWithOp               *Uy
	,MatrixWithOp               *Vy
	,EJournalOutputLevel        olevel
	,std::ostream               &out
	)
{
	namespace rcp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	using LinAlgOpPack::syrk;

	const size_type
		n = nlp.n(),
		m = nlp.m(),
		r = nlp.r();
	const Range1D
		var_dep(1,r),
		var_indep(r+1,n),
		con_decomp   = nlp.con_decomp(),
		con_undecomp = nlp.con_undecomp();

	//
	// Get pointers to concreate matrices
	//
	
	MatrixIdentConcatStd
		*Y_orth = Y ? &dyn_cast<MatrixIdentConcatStd>(*Y)  : NULL;
	MatrixCompositeStd
		*Uy_cpst = Uy ? &dyn_cast<MatrixCompositeStd>(*Uy) : NULL;			
	MatrixCompositeStd
		*Vy_cpst = Vy ? &dyn_cast<MatrixCompositeStd>(*Vy) : NULL;

	//
	// Initialize the matrices
	//

	// Y
	if(Y_orth) {
		D_ptr_t  D_ptr = D;
//		if(mat_rel == MATRICES_INDEP_IMPS) {
//			D_ptr = D->clone();
//			THROW_EXCEPTION(
//				D_ptr.get() == NULL, std::logic_error
//				,"DecompositionSystemOrthogonal::update_decomp(...) : Error, "
//				"The matrix class used for the direct sensitivity matrix D = inv(C)*N of type \'"
//				<< typeid(*D).name() << "\' must return return.get() != NULL from the clone() method "
//				"since mat_rel == MATRICES_INDEP_IMPS!" );
//		}
		Y_orth->initialize(
			nlp.space_x()                                     // space_cols
			,nlp.space_x()->sub_space(var_dep)->clone()       // space_rows
			,MatrixIdentConcatStd::BOTTOM                     // top_or_bottom
			,-1.0                                             // alpha
			,D_ptr                                            // D_ptr
			,BLAS_Cpp::no_trans                               // D_trans
			);
	}

	// S
	if(S_ptr_.get() == NULL) {
		S_ptr_ = nlp.factory_S()->create();
	}
	// S = I + (D)'*(D')'
	dyn_cast<MatrixSymInitDiagonal>(*S_ptr_).init_identity(D->space_rows());
	syrk(*D,BLAS_Cpp::trans,1.0,1.0,S_ptr_.get());

	assert(Uy_cpst == NULL); // ToDo: Implement for undecomposed equalities
	assert(Vy_cpst == NULL); // ToDo: Implement for general inequalities

	recalc_py(*D,py,olevel,out);

}

void EvalNewPointTailoredApproachOrthogonal_Step::recalc_py(
	const MatrixWithOp       &D
	,VectorWithOpMutable     *py
	,EJournalOutputLevel     olevel
	,std::ostream            &out
	)
{
	using BLAS_Cpp::no_trans;
	using BLAS_Cpp::trans;
	using AbstractLinAlgPack::Vp_StMtV;
	using AbstractLinAlgPack::V_InvMtV;
	using LinAlgOpPack::V_MtV;

	const MatrixSymWithOpNonsingular   &S = *S_ptr_;

	VectorSpace::vec_mut_ptr_t               // ToDo: make workspace!
		tIa = D.space_rows().create_member(),
		tIb = D.space_rows().create_member();
	//
	// y = (I - D*inv(S)*D')*inv(C)*x
	//   = (I - D*inv(S)*D')*py
	//   = py - D*inv(S)*D'*py
	//
	// =>
	//
	// tIa  = D'*py
	// tIb  = inv(S)*tIa
	// py   += -D*tIb
	//
	V_MtV( tIa.get(), D, trans, *py );        // tIa  = D'*py
	V_InvMtV( tIb.get(), S, no_trans, *tIa ); // tIb  = inv(S)*tIa
	Vp_StMtV( py, -1.0, D, no_trans, *tIb );  // y   += -D*tIb

}

void EvalNewPointTailoredApproachOrthogonal_Step::print_calc_py_Y_Uy_Vy(
	std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Orthogonal decomposition\n"
		<< L << "py = inv(I + D*D') * py <: space_range\n"
		<< L << "Y = [ I ; -D' ] <: space_x|space_range\n"
		<< L << "Uy = ???\n"
		<< L << "Vy = ???\n"
		;
}

}	// end namespace ReducedSpaceSQPPack 
