// ///////////////////////////////////////////////////////
// BFGSUpdate_Strategy.cpp
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

#include "ReducedSpaceSQPPack/include/std/BFGSUpdate_Strategy.h"
#include "ReducedSpaceSQPPack/include/ReducedSpaceSQPPackExceptions.h"
#include "AbstractLinAlgPack/include/TestMatrixSymSecantUpdate.h"
#include "AbstractLinAlgPack/include/MatrixSymSecantUpdateable.h"
#include "AbstractLinAlgPack/include/MatrixSymWithOp.h"
#include "AbstractLinAlgPack/include/MatrixWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorSpace.h"
#include "AbstractLinAlgPack/include/VectorWithOpMutable.h"
#include "AbstractLinAlgPack/include/VectorWithOpOut.h"
#include "AbstractLinAlgPack/include/VectorStdOps.h"
#include "AbstractLinAlgPack/include/LinAlgOpPack.h"
#include "dynamic_cast_verbose.h"
#include "ThrowException.h"

namespace ReducedSpaceSQPPack {

BFGSUpdate_Strategy::BFGSUpdate_Strategy(
	bool               rescale_init_identity
	,bool              use_dampening
	,ESecantTesting    secant_testing
	,value_type        secant_warning_tol
	,value_type        secant_error_tol
	)
	:rescale_init_identity_(rescale_init_identity)
	,use_dampening_(use_dampening)
	,secant_testing_(secant_testing)
	,secant_warning_tol_(secant_warning_tol)
	,secant_error_tol_(secant_error_tol)
{}

void BFGSUpdate_Strategy::perform_update(
	VectorWithOpMutable      *s_bfgs
	,VectorWithOpMutable     *y_bfgs
	,bool                    first_update
	,std::ostream            &out
	,EJournalOutputLevel     olevel
	,bool                    check_results
	,MatrixSymWithOp         *B
	,QuasiNewtonStats        * quasi_newton_stats 
	)
{
	namespace rcp = MemMngPack;
	using DynamicCastHelperPack::dyn_cast;
	using AbstractLinAlgPack::dot;
	using AbstractLinAlgPack::Vt_S;
	using AbstractLinAlgPack::Vp_StV;
	using LinAlgOpPack::Vp_V;
	using LinAlgOpPack::V_StV;
	using LinAlgOpPack::V_VpV;
	using LinAlgOpPack::V_VmV;
	using LinAlgOpPack::V_MtV;

	const size_type
		nind = B->rows();
	const value_type
		NOT_CALCULATED = -1.0/3.0;
	value_type
		sTy = NOT_CALCULATED,
		yTy = NOT_CALCULATED;
	bool used_dampening = false;

	MatrixSymSecantUpdateable
	    &B_updatable = dyn_cast<MatrixSymSecantUpdateable>(*B);

	// /////////////////////////////////////////////////////////////
	// Consider rescaling the initial identity hessian before
	// the update.
	// 
	// This was taken from Nocedal & Wright, 1999, p. 200 
	// 
	// Bo = (y'*y)/(y'*s) * I
	// 
	if( first_update && rescale_init_identity() ) {
		if( sTy == NOT_CALCULATED )
			sTy = dot( *s_bfgs, *y_bfgs );
		if( yTy == NOT_CALCULATED )
			yTy = dot( *y_bfgs, *y_bfgs );
		const value_type
			Iscale = yTy/sTy; 
		const value_type
			Iscale_too_small = 1e-5;  // ToDo: Make this adjustable
		if( Iscale >= Iscale_too_small ) {	
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nRescaling the initial identity matrix before the update as:\n"
					<< "Iscale = (y'*y)/(y'*s) = ("<<yTy<<")/("<<sTy<<") = "<<Iscale<<"\n"
					<< "B =  Iscale * eye(n-r) ...\n";
			}
			B_updatable.init_identity( y_bfgs->space(), Iscale );
			if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
				out << "\nB after rescaling = \n" << *B;
			}
		}
		else {
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
				out << "\nSkipping rescaling of the initial identity matrix since:\n"
					<< "Iscale = (y'*y)/(y'*s) = ("<<yTy<<")/("<<sTy<<") = "<<Iscale
					<< " < Iscale_too_small = " << Iscale_too_small << std::endl;
			}
		}
	}

	// ////////////////////////////////////////////////////
	// Modify the s_bfgs and y_bfgs vectors for dampening?
	VectorSpace::vec_mut_ptr_t
		Bs = rcp::null;
	if( use_dampening() ) {
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
		{
			out
				<< "\nConsidering the dampened update ...\n";
		}
		// Bs = Bm1 * s_bfgs
		Bs = y_bfgs->space().create_member();
		V_MtV( Bs.get(), *B, BLAS_Cpp::no_trans, *s_bfgs );
		// sTy = s_bfgs' * y_bfgs
		if( sTy == NOT_CALCULATED )
			sTy = dot( *s_bfgs, *y_bfgs );
		// sTBs = s_bfgs' * Bs
		const value_type sTBs = dot( *s_bfgs, *Bs );
		// Select dampening parameter theta
		const value_type
			theta = ( sTy >= 0.2 * sTBs )
			? 1.0
			: (0.8 * sTBs ) / ( sTBs - sTy );
		if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
		{
			out
				<< "\ns_bfgs'*y_bfgs = " << sTy
				<< ( theta == 1.0 ? " >= " : " < " )
				<< " s_bfgs' * B * s_bfgs = " << sTBs << std::endl;
		}
		if( theta == 1.0 ) {
			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
			{
				out
					<< "Perform the undamped update ...\n";
			}
		}
		else {
			// y_bfgs = theta*y_bfgs + (1-theta)*B*s_bfgs
			Vt_S( y_bfgs, theta );
			Vp_StV( y_bfgs, (1.0-theta), *Bs );

			if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS )
			{
				out
					<< "Dampen the update ...\n"
					<< "theta = " << theta << std::endl
					<< "y_bfgs = theta*y_bfgs + (1-theta)*B*s_bfgs ...\n"
					<< "||y_bfgs||inf = " << y_bfgs->norm_inf() << std::endl;
			}

			if( (int)olevel >= (int)PRINT_VECTORS )
			{
				out
					<< "y_bfgs =\n" << *y_bfgs;
			}

			used_dampening = true;
		}
	}

	// Perform the update if it is defined (s_bfgs' * y_bfgs > 0.0)
				
	VectorSpace::vec_mut_ptr_t
		s_bfgs_save = rcp::null,
		y_bfgs_save = rcp::null;
	if( check_results ) {
		// Save s and y since they may be overwritten in the update.
		s_bfgs_save = s_bfgs->clone();
		y_bfgs_save = y_bfgs->clone();
	}
	try {
		B_updatable.secant_update(
			s_bfgs
			,y_bfgs
			,Bs.get()
			);
	}
	catch( const MatrixSymSecantUpdateable::UpdateSkippedException& excpt ) {
		if( (int)olevel >= (int)PRINT_BASIC_ALGORITHM_INFO ) {
			out << excpt.what() << std::endl
				<< "\nSkipping BFGS update.  B = B ...\n";
		}
		quasi_newton_stats->set_updated_stats(
			QuasiNewtonStats::INDEF_SKIPED );
		return;
	}					
		
	if( (int)olevel >= (int)PRINT_ITERATION_QUANTITIES ) {
		out << "\nB after the BFGS update = \n" << *B;
	}
	
	if( ( check_results && secant_testing() == SECANT_TEST_DEFAULT )
		|| secant_testing() == SECANT_TEST_ALWAYS )
	{
		const bool result =
			AbstractLinAlgPack::TestMatrixSymSecantUpdate(
				*B, *s_bfgs_save, *y_bfgs_save
				, secant_warning_tol(), secant_error_tol()
				, (int)olevel >= (int)PRINT_VECTORS
				, (int)olevel >  (int)PRINT_NOTHING ? &out : NULL
				, (int)olevel >= (int)PRINT_ALGORITHM_STEPS
				);
		if( !result ) {
			const char
				msg[] =	"Error, the secant property for the BFGS update failed\n"
				"Stopping the algorithm ...\n";
			out << msg;
			THROW_EXCEPTION(
				true, TestFailed
				," BFGSUpdate_Strategy::perform_update(...) : " << msg );
		}
	}
	
	quasi_newton_stats->set_updated_stats(
		used_dampening
		? QuasiNewtonStats::DAMPENED_UPDATED
		: QuasiNewtonStats::UPDATED );
}

void BFGSUpdate_Strategy::print_step( std::ostream& out, const std::string& L ) const
{
	out
		<< L << "if use_dampening == true then\n"
		<< L << "    if s'*y >= 0.2 * s'*B*s then\n"
		<< L << "        theta = 1.0\n"
		<< L << "    else\n"
		<< L << "        theta = 0.8*s'*B*s / (s'*B*s - s'*y)\n"
		<< L << "    end\n"
		<< L << "    y = theta*y + (1-theta)*B*s\n"
		<< L << "end\n"
		<< L << "if first_update && rescale_init_identity and y'*s is sufficently positive then\n"
		<< L << "    B = |(y'*y)/(y'*s)| * eye(size(s))\n"
		<< L << "end\n"
		<< L << "if s'*y is sufficently positive then\n"
		<< L << "    *** Peform BFGS update\n"
		<< L << "    (B, s, y ) -> B (through MatrixSymSecantUpdate interface)\n"
		<< L << "    if ( check_results && secant_testing == SECANT_TEST_DEFAULT )\n"
		<< L << "    or ( secant_testing == SECANT_TEST_ALWAYS ) then\n"
		<< L << "        if B*s != y then\n"
		<< L << "            *** The secant condition does not check out\n"
		<< L << "            Throw TestFailed!\n"
		<< L << "        end\n"
		<< L << "    end\n"
		<< L << "end\n"
		;
}

}  // end namespace ReducedSpaceSQPPack
