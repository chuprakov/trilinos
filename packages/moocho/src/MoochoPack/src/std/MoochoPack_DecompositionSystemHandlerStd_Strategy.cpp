// ////////////////////////////////////////////////////////////////////////////
// DecompositionSystemHandlerStd_Strategy.cpp
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

#include <ostream>
#include <typeinfo>

#include "ReducedSpaceSQPPack/src/std/DecompositionSystemHandlerStd_Strategy.h"
#include "ReducedSpaceSQPPack/src/ReducedSpaceSQPPackExceptions.h"
#include "ReducedSpaceSQPPack/src/rsqp_algo_conversion.h"
#include "GeneralIterationPack/src/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/src/DecompositionSystem.h"
#include "NLPInterfacePack/src/NLPFirstOrderInfo.h"
#include "AbstractLinAlgPack/src/MatrixWithOpNonsingular.h"
#include "AbstractLinAlgPack/src/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/src/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/src/VectorStdOps.h"
#include "AbstractLinAlgPack/src/VectorWithOpOut.h"
#include "AbstractLinAlgPack/src/assert_print_nan_inf.h"
#include "AbstractLinAlgPack/src/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace ReducedSpaceSQPPack {

// Constructors / initializers

DecompositionSystemHandlerStd_Strategy::DecompositionSystemHandlerStd_Strategy()
{}

// Overridden from DecompositionSystemHandler_Strategy

bool DecompositionSystemHandlerStd_Strategy::update_decomposition(
	rSQPAlgo                                &algo
	,rSQPState                              &s
	,NLPFirstOrderInfo                      &nlp
	,EDecompSysTesting                      decomp_sys_testing
	,EDecompSysPrintLevel                   decomp_sys_testing_print_level
	,bool                                   *new_decomp_selected
	)
{
	using DynamicCastHelperPack::dyn_cast;

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	const size_type
		n  = nlp.n(),
		nb = nlp.num_bounded_x(),
		m  = nlp.m(),
		mI = nlp.mI();
	size_type
		r  = s.decomp_sys().equ_decomp().size();

	// Get the iteration quantity container objects
	IterQuantityAccess<index_type>
		&num_basis_iq = s.num_basis();
	IterQuantityAccess<VectorWithOpMutable>
		&x_iq   = s.x(),
		&nu_iq  = s.nu();
	IterQuantityAccess<MatrixWithOp>
		*Gc_iq  = m  > 0                  ? &s.Gc() : NULL,
		*Gh_iq  = mI > 0                  ? &s.Gh() : NULL,
		*Z_iq   = ( n > m && r > 0 )      ? &s.Z()  : NULL,
		*Y_iq   = ( r > 0 )               ? &s.Y()  : NULL,
		*Uz_iq  = ( m  > 0 && m  > r )    ? &s.Uz() : NULL,
		*Uy_iq  = ( m  > 0 && m  > r )    ? &s.Uy() : NULL,
		*Vz_iq  = ( mI > 0 ) && ( m > 0 ) ? &s.Vz() : NULL,
		*Vy_iq  = ( mI > 0 ) && ( m > 0 ) ? &s.Vy() : NULL;
	IterQuantityAccess<MatrixWithOpNonsingular>
		*R_iq   = ( m > 0 )               ? &s.R()  : NULL;
	
	//
	// Update range/null decomposition
	//

	// Determine if we will test the decomp_sys or not
	const DecompositionSystem::ERunTests
		ds_test_what = ( ( decomp_sys_testing == DST_TEST
						   || ( decomp_sys_testing == DST_DEFAULT
								&& algo.algo_cntr().check_results() ) )
						 ? DecompositionSystem::RUN_TESTS
						 : DecompositionSystem::NO_TESTS );
		
	// Determine the output level for decomp_sys				
	DecompositionSystem::EOutputLevel ds_olevel;
	switch(olevel) {
		case PRINT_NOTHING:
		case PRINT_BASIC_ALGORITHM_INFO:
			ds_olevel = DecompositionSystem::PRINT_NONE;
			break;
		case PRINT_ALGORITHM_STEPS:
		case PRINT_ACTIVE_SET:
			ds_olevel = DecompositionSystem::PRINT_BASIC_INFO;
			break;
		case PRINT_VECTORS:
			ds_olevel = DecompositionSystem::PRINT_VECTORS;
			break;
		case PRINT_ITERATION_QUANTITIES:
			ds_olevel = DecompositionSystem::PRINT_EVERY_THING;
			break;
		default:
			assert(0); // Should not get here!
	};

	// Form the decomposition of Gc and Gh and update the decomposition system matrices
	if( olevel >= PRINT_ALGORITHM_STEPS ) {
		out << "\nUpdating the range/null decompostion matrices ...\n";
	}
	s.decomp_sys().update_decomp(
		&out                               // out
		,ds_olevel                         // olevel
		,ds_test_what                      // test_what
		,Gc_iq->get_k(0)                   // Gc
		,Gh_iq ? &Gh_iq->get_k(0) : NULL   // Gh
		,&Z_iq->set_k(0)                   // Z
		,&Y_iq->set_k(0)                   // Y
		,&R_iq->set_k(0)                   // R
		,Uz_iq ? &Uz_iq->set_k(0) : NULL   // Uz
		,Uy_iq ? &Uy_iq->set_k(0) : NULL   // Uy
		,Vz_iq ? &Vz_iq->set_k(0) : NULL   // Vz
		,Vy_iq ? &Vy_iq->set_k(0) : NULL   // Vy
		,DecompositionSystem::MATRICES_ALLOW_DEP_IMPS // ToDo: Change this!
		);
	s.equ_decomp(   s.decomp_sys().equ_decomp()   );
	s.equ_undecomp( s.decomp_sys().equ_undecomp() );

	*new_decomp_selected = false;

	return true;
}

void DecompositionSystemHandlerStd_Strategy::print_update_decomposition(
	const rSQPAlgo                          &algo
	,const rSQPState                        &s
	,std::ostream                           &out
	,const std::string                      &L
	) const
{
	out
		<< L << "*** Updating the range/null decomposition.\n"
		<< L << "begin update decomposition\n"
		<< L << "(class = \'" << typeid(s.decomp_sys()).name() << "\')\n"
		;
	s.decomp_sys().print_update_decomp( out, L + "  " );
	out
		<< L << "end update decomposition\n"
		;
}

} // end namespace ReducedSpaceSQPPack
