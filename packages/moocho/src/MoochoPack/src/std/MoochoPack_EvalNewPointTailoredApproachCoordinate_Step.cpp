// /////////////////////////////////////////////////////////////////////////////////
// EvalNewPointTailoredApproachCoordinate_Step.cpp
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

#include "ReducedSpaceSQPPack/include/std/EvalNewPointTailoredApproachCoordinate_Step.h"
#include "ConstrainedOptimizationPack/include/MatrixIdentConcatStd.h"
#include "NLPInterfacePack/include/NLPFirstOrderDirect.h"
#include "AbstractLinAlgPack/include/MatrixWithOp.h"
#include "AbstractLinAlgPack/include/MatrixZero.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

EvalNewPointTailoredApproachCoordinate_Step::EvalNewPointTailoredApproachCoordinate_Step(
	const deriv_tester_ptr_t& 	  deriv_tester
	,const bounds_tester_ptr_t&	  bounds_tester
	,EFDDerivTesting              fd_deriv_testing
	)
	:EvalNewPointTailoredApproach_Step(deriv_tester,bounds_tester,fd_deriv_testing)
{}

// protected

void EvalNewPointTailoredApproachCoordinate_Step::uninitialize_Y_Uv_Uy(
	MatrixWithOp         *Y
	,MatrixWithOp        *Uy
	,MatrixWithOp        *Vy
	)
{
	// Nothing to free
}

void EvalNewPointTailoredApproachCoordinate_Step::calc_py_Y_Uy_Vy(
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
	namespace rcp = ReferenceCountingPack;
	using DynamicCastHelperPack::dyn_cast;

	MatrixIdentConcatStd
		cY = dyn_cast<MatrixIdentConcatStd>(*Y);
	//
	// Y = [      I     ] space_xD  
	//     [    Zero    ] space_xI
	//        space_xD
	//
	VectorSpace::space_ptr_t
		space_x  = nlp.space_x(),
		space_xD = space_x->sub_space(nlp.var_dep())->clone(),
		space_xI = space_x->sub_space(nlp.var_indep())->clone();
	cY.initialize(
		space_x                                                // space_cols
		,space_xD                                              // space_rows
		,MatrixIdentConcatStd::BOTTOM                          // top_or_bottom
		,1.0                                                   // alpha
		,rcp::rcp(
			new MatrixZero(
				space_xI    // space_cols
				,space_xD   // space_rows
				) )                                            // D_ptr
		,BLAS_Cpp::no_trans                                    // D_trans
		);
	// py is not altered here!
}

void EvalNewPointTailoredApproachCoordinate_Step::recalc_py(
	const MatrixWithOp       &D
	,VectorWithOpMutable     *py
	,EJournalOutputLevel     olevel
	,std::ostream            &out
	)
{
	// py is not altered here!
}

void EvalNewPointTailoredApproachCoordinate_Step::print_calc_py_Y_Uy_Vy(
	std::ostream& out, const std::string& L
	) const
{
	out
		<< L << "*** Coordinate decomposition\n"
		<< L << "py_k = py_k\n"
		<< L << "Y = [ I ; 0 ] <: R^(n x m) [0 represented using MatrixZero]\n"
		<< L << "Uy = Gc(var_dep,con_undecomp)\'\n"
		<< L << "Vy = Gh(var_dep,:)\'\n"
		;
}

}	// end namespace ReducedSpaceSQPPack 
