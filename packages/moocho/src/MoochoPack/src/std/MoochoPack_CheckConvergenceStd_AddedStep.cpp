// ////////////////////////////////////////////////////////////////////////////
// CheckConvergenceStd_AddedStep.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>

#include <ostream>
#include <limits>
#include <sstream>

#include "../../include/std/CheckConvergenceStd_AddedStep.h"
#include "../../include/rSQPAlgoContainer.h"
#include "../../include/rsqp_algo_conversion.h"
#include "GeneralIterationPack/include/print_algorithm_step.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
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

namespace ReducedSpaceSQPPack {

CheckConvergenceStd_AddedStep::CheckConvergenceStd_AddedStep(
		  EOptErrorCheck opt_error_check
		, EScaleKKTErrorBy scale_kkt_error_by	)
	 : opt_error_check_(opt_error_check), scale_kkt_error_by_(scale_kkt_error_by)
{}

bool CheckConvergenceStd_AddedStep::do_step(Algorithm& _algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss)
{
	using LinAlgPack::norm_inf;
	using LinAlgPack::norm_2;

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

	// scale_kkt_factor
	value_type
		scale_kkt_factor;
	switch(scale_kkt_error_by()) {
		case SCALE_BY_ONE:
			scale_kkt_factor = 1.0;
			break;
		case SCALE_BY_NORM_2_X:
			scale_kkt_factor = 1.0 + s.x().get_k(0).norm_2();
			break;
		case SCALE_BY_NOMR_INF_X:
			scale_kkt_factor = 1.0 + s.x().get_k(0).norm_inf();
			break;
		default:
			assert(0);	// Should never be called
	}


	// opt_err = (||rGL||inf or ||GL||) / (||Gf|| + scale_kkt_factor)
	const value_type
		opt_err =
			( opt_error_check() == OPT_ERROR_REDUCED_GRADIENT_LAGR ?
				assert_is_number( s.rGL().get_k(0).norm_inf() , "||rGL_k||inf" )
				: assert_is_number( s.GL().get_k(0).norm_inf() , "||GL_k||inf" )
			) / (1.0 + assert_is_number( s.Gf().get_k(0).norm_inf(), "||Gf_k||inf" ));

	// feas_err
	const value_type
		feas_err = assert_is_number( s.c().get_k(0).norm_inf(), "||c_k||inf" );

	// kkt_err
	s.kkt_err().set_k(0) = opt_err + feas_err;

	// step_err
	value_type
		step_err = 0.0;
	if( s.d().updated_k(0) ) {
		const Vector
			&d = s.d().get_k(0).v(),
			&x = s.x().get_k(0).v();
		Vector::const_iterator
			d_itr = d.begin(),
			x_itr = x.begin();
		while( d_itr != d.end() )
			step_err = std::_MAX( step_err, ::fabs(*d_itr++)/(1.0+::fabs(*x_itr++)) );
	}

	assert_is_number( step_err, "max(d(i)/max(1,x(i)),i=1...n)" );

	if( static_cast<int>(olevel) >= static_cast<int>(PRINT_ALGORITHM_STEPS) ) {
		out	<< "\nscale_kkt_factor = " << scale_kkt_factor
			<< "\nopt_err          = " << opt_err
			<< "\nfeas_err         = " << feas_err
			<< "\nkkt_err          = " << s.kkt_err().get_k(0)
			<< "\nstep_err         = " << step_err << std::endl;
	}

	const value_type
		opt_tol		= algo.algo_cntr().opt_tol(),
		feas_tol	= algo.algo_cntr().feas_tol(),
		step_tol	= algo.algo_cntr().step_tol();

	if(    opt_err/scale_kkt_factor < opt_tol
		&& feas_err/scale_kkt_factor < feas_tol
		&& step_err < step_tol )
	{

		if( static_cast<int>(olevel) >= static_cast<int>(PRINT_BASIC_ALGORITHM_INFO) ) {
			out	<< "\nFound the solution!!!!!! (k = " << algo.state().k() << ")"
				<< "\nopt_err   = " << opt_err	<< " < opt_tol = "	<< opt_tol
				<< "\nfeas_err  = " << feas_err	<< " < feas_tol = "	<< feas_tol
				<< "\nstep_err  = " << step_err	<< " < step_tol = "	<< step_tol
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

void CheckConvergenceStd_AddedStep::print_step( const Algorithm& algo
	, poss_type step_poss, GeneralIterationPack::EDoStepType type, poss_type assoc_step_poss
	, std::ostream& out, const std::string& L ) const
{
	out
		<< L << "*** Check to see if the KKT error is small enough for convergence ***\n"
		<< L << "if scale_kkt_error_by == SCALE_BY_ONE then\n"
		<< L << "    scale_kkt_factor = 1.0\n"
		<< L << "else if scale_by == SCALE_BY_NORM_2_X then\n"
		<< L << "    scale_kkt_factor = 1.0 + norm_2(x_k)\n"
		<< L << "else if scale_by == SCALE_BY_NOMR_INF_X then\n"
		<< L << "    scale_kkt_factor = 1.0 + norm_inf(x_k)\n"
		<< L << "else\n"
		<< L << "norm_inf_c_k = norm(c_k,inf)\n"
		<< L << "opt_err = (norm_inf(rGL_k) or norm_inf(GL_k))\n"
		<< L << "    / (1.0 + norm_inf(Gf_k))\n"
		<< L << "feas_err = norm_inf_c_k\n"
		<< L << "kkt_err_k = opt_err + feas_err\n"
		<< L << "if d_k is updated then\n"
		<< L << "    step_err = max( |d_k(i)|/(1+|x_k(i)|), i=1..n )\n"
		<< L << "else\n"
		<< L << "    step_err = 0\n"
		<< L << "end\n"
		<< L << "if opt_err/scale_kkt_factor < opt_tol\n"
		<< L << "       and feas_err/scale_kkt_factor < feas_tol\n"
		<< L << "       and step_err < step_tol then\n"
		<< L << "   report optimal x_k, lambda_k and nu_k to the nlp\n"
		<< L << "   terminate, the solution has beed found!\n"
		<< L << "end\n";
}

}	// end namespace ReducedSpaceSQPPack