// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <ostream>
#include <limits>
#include <sstream>

#include "../../include/std/CheckConvergenceStd_AddedStep.h"
#include "../../include/rSQPAlgoContainer.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/VectorOut.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"

namespace {

// Assert that we have a finite number.  It it is not then throw
// an std::runtime_error exception.
ReducedSpaceSQPPack::value_type
assert_is_number( ReducedSpaceSQPPack::value_type val, const char name[] )
{
	typedef std::numeric_limits<ReducedSpaceSQPPack::value_type> nl_t;
	if( val >= nl_t::max() || val != val ) {
		std::ostringstream omsg;
		omsg
			<<	"CheckConvergenceStd_AddedStep::do_step(...) : Error, "
			<< name << " = " << val << " is not a finite number.";
		throw std::runtime_error( omsg.str() );
	}
	return val;
}

}	// end namespace 

bool ReducedSpaceSQPPack::CheckConvergenceStd_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_inf;

	rSQPAlgo	&algo	= rsqp_algo(_algo);
	rSQPState	&s		= algo.rsqp_state();
	NLPReduced	&nlp	= algo.nlp();

	EIterationInfoOutput olevel = s.iteration_info_output();
	std::ostream& out = algo.track().journal_out();

	// print step header.
	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		using GeneralIterationPack::print_algorithm_step;
		print_algorithm_step( algo, step_poss, type, assoc_step_poss, out );
	}

	// kkt_error = max( ||rGL||inf, ||c||inf )
	const value_type
		kkt_error = ( assert_is_number( s.norm_inf_rGL().set_k(0)
						 = norm_inf(s.rGL().get_k(0)()) , "||rGL_k||inf" ) )
					/ std::_MAX( 1.0
						, assert_is_number( norm_inf(s.Gf().get_k(0)), "||Gf_k||inf" ) );
	// feas_error
	const value_type
		feas_error = assert_is_number( s.norm_inf_c().set_k(0) = norm_inf(s.c().get_k(0)())
						, "||c_k||inf" );

	// step_error
	value_type
		step_error = 0.0;
	if( s.d().updated_k(0) ) {
		const Vector
			&d = s.d().get_k(0),
			&x = s.x().get_k(0);
		Vector::const_iterator
			d_itr = d.begin(),
			x_itr = x.begin();
		while( d_itr != d.end() )
			step_error = std::_MAX( step_error, ::fabs(*d_itr++)/(1.0+::fabs(*x_itr++)) );
	}

	assert_is_number( step_error, "max(d(i)/max(1,x(i)),i=1...n)" );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nkkt_error   = " << kkt_error
			<< "\nfeas_error  = " << feas_error
			<< "\nstep_error  = " << step_error << std::endl;
	}

	const value_type
		kkt_tol		= algo.algo_cntr().kkt_tol(),
		feas_tol	= algo.algo_cntr().feas_tol(),
		step_tol	= algo.algo_cntr().step_tol();

	if( kkt_error < kkt_tol && feas_error < feas_tol && step_error < step_tol ) {

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nFound the solution!!!!!! (k = " << algo.state().k() << ")"
				<< "\nkkt_error   = " << kkt_error	<< " < kkt_tol = "	<< kkt_tol
				<< "\nfeas_error  = " << feas_error	<< " < feas_tol = "	<< feas_tol
				<< "\nstep_error  = " << step_error	<< " < step_tol = "	<< step_tol
					<< std::endl;
		}
		nlp.report_final_solution(
			  s.x().get_k(0)()
			, s.lambda().updated_k(0)	? &s.lambda().get_k(0)()	: 0
			, s.nu().updated_k(0)		? &s.nu().get_k(0)()		: 0
			, true
			);
		algo.terminate(true);	// found min
		return false; // skip the other steps and terminate
	}

	// We are not at the solution so keep going
	return true;
}

void ReducedSpaceSQPPack::CheckConvergenceStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Check to see if the KKT error is small enough for convergence ***\n"
		<< L << "norm_inf_rGL_k = norm(rGL_k,inf)\n"
		<< L << "norm_inf_c_k = norm(c_k,inf)\n"
		<< L << "kkt_err = norm_inf_rGL_k / max(1.0,norm_inf(Gf_k))\n"
		<< L << "feas_err = norm_inf_c_k\n"
		<< L << "if d_k is updated then\n"
		<< L << "    step_err = max( |d_k(i)|/(1+|x_k(i)|), i=1..n )\n"
		<< L << "else\n"
		<< L << "    step_err = 0\n"
		<< L << "end\n"
		<< L << "if kkt_err < kkt_tol and feas_err < feas_tol and step_err < step_tol then\n"
		<< L << "   report optimal x_k, lambda_k and nu_k to the nlp\n"
		<< L << "   terminate, the solution has beed found!\n"
		<< L << "end\n";
}
