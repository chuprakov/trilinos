// ////////////////////////////////////////////////////////////////////////////
// LineSearchDirect_Step.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <typeinfo>

#include "../../include/std/LineSearchDirect_Step.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/MeritFuncCalc1DQuadratic.h"
#include "ConstrainedOptimizationPack/include/MeritFuncCalcNLP.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"

bool ReducedSpaceSQPPack::LineSearchDirect_Step::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_inf;
	using LinAlgPack::V_VpV;
	using LinAlgPack::Vp_StV;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLPReduced	&nlp	= algo.nlp();

	EJournalOutputLevel olevel = algo.algo_cntr().journal_output_level();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// /////////////////////////////////////////
	// Set references to iteration quantities
	//
	// Set k+1 first then go back to get k to ensure
	// we have backward storage.
	
	Vector
		&x_kp1 = s.x().set_k(+1).v();
	value_type
		&f_kp1 = s.f().set_k(+1);
	Vector
		&c_kp1 = s.c().set_k(+1).v();

	const value_type
		&f_k = s.f().get_k(0);
	const Vector
		&c_k = s.c().get_k(0).v();
	const Vector
		&x_k = s.x().get_k(0).v();
	const Vector
		&d_k = s.d().get_k(0).v();
	value_type
		&alpha_k = s.alpha().get_k(0);

	// /////////////////////////////////////
	// Compute Dphi_k, phi_kp1 and phi_k

	// Dphi_k
	const value_type
		Dphi_k = merit_func().deriv();
	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "\nDphi_k = "	<< Dphi_k << std::endl;
	}

	if( Dphi_k >= 0 ) {
		throw LineSearchFailure( "LineSearch2ndOrderCorrect_Step::do_step(...) : " 
			"Error, d_k is not a descent direction for the merit function " );
	}

	// ph_kp1
	value_type
		&phi_kp1 = s.phi().set_k(+1) = merit_func().value( f_kp1, c_kp1 );

	// Must compute phi(x) at the base point x_k since the penalty parameter may have changed.
	const value_type
		&phi_k = s.phi().set_k(0) = merit_func().value( f_k, c_k );

	// //////////////////////////////////////
	// Setup the calculation merit function

	// Here f_kp1, and c_kp1 are updated at the same time the
	// line search is being performed.
	nlp.set_f( &f_kp1 );
	nlp.set_c( &c_kp1 );
	MeritFuncCalcNLP
		phi_calc( &merit_func(), &nlp );

	// //////////////////////
	// Do the line search
	
	const VectorSlice xd[2] = { x_k(), d_k() };
	MeritFuncCalc1DQuadratic
		phi_calc_1d( phi_calc, 1, xd, &x_kp1() );

	if( !direct_line_search().do_line_search( phi_calc_1d, phi_k, &alpha_k, &phi_kp1
		, static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ?
			&out : static_cast<std::ostream*>(0)	)		)
	{
		// If the line search failed but the value of the merit function is less than
		// the base point value then just accept it and move on.  This many be a
		// bad thing to do.

		const value_type
			scaled_ared		= (phi_k - phi_kp1)/phi_k,
			keep_on_frac	= 1.0e-10;	// Make adjustable?
		bool keep_on = scaled_ared < keep_on_frac;

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) )
		{
			out
				<< "\nThe maximum number of linesearch iterations has been exceeded "
				<< "(k = " << algo.state().k() << ")\n"
				<< "(phi_k - phi_kp1)/phi_k = " << scaled_ared;
//			if(keep_on) {
//				out
//					<< " < " << keep_on_frac
//					<< "\nso we will accept to step and move on.\n";
//			}
//			else {
				out
//					<< " > " << keep_on_frac
					<< "\nso we will reject the step and declare a line search failure.\n";
//			}
		}
//
//		if( keep_on ) return true;
		
		throw LineSearchFailure("LineSearchDirect_Step::do_step(): Line search failure");
	}

	if( (int)olevel >= (int)PRINT_ALGORITHM_STEPS ) {
		out	<< "\nalpha_k      = "	<< alpha_k				<< std::endl
			<< "\n||x_kp1||inf = "	<< norm_inf( x_kp1 )	<< std::endl
			<< "\nf_kp1        = "	<< f_kp1				<< std::endl
			<< "\n||c_kp1||inf = "	<< norm_inf(c_kp1)		<< std::endl
			<< "\nphi_kp1      = "	<< phi_kp1				<< std::endl;
	}

	if( (int)olevel >= (int)PRINT_VECTORS ) {
		out << "\nx_kp1 =\n"	<< x_kp1
			<< "\nc_kp1 =\n"	<< c_kp1;
	}

	return true;
}

void ReducedSpaceSQPPack::LineSearchDirect_Step::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Preform a line search along the full space search direction d_k.\n"
		<< L << "begin definition of NLP merit function phi.value(f(x),c(x)):\n";

	merit_func().print_merit_func( out, L + "    " );
	
	out	<< L << "end definition\n"
		<< L << "Dphi_k = phi.deriv()\n"
		<< L << "if Dphi_k >= 0 then\n"
		<< L << "    throw line_search_failure\n"
		<< L << "end\n"
		<< L << "phi_kp1 = phi_k.value(f_kp1,c_kp1)\n"
		<< L << "phi_k = phi.value(f_k,c_k)\n"
		<< L << "begin direct line search : \"" << typeid(direct_line_search()).name() << "\"\n";

	direct_line_search().print_algorithm( out, L + "    " );

	out
		<< L << "end direct line search\n"
		<< L << "if maximum number of linesearch iterations are exceeded then\n"
		<< L << "    throw line_search_failure\n"
		<< L << "end\n";
}