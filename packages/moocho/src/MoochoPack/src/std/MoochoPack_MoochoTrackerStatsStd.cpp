// ////////////////////////////////////////////////////////////////////////////
// rSQPTrackStatsStd.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <iomanip>

#include "../../include/std/rSQPTrackStatsStd.h"
#include "../../include/rSQPState.h"
#include "../../include/rsqp_algo_conversion.h"
#include "ConstrainedOptimizationPack/include/VectorWithNorms.h"
#include "NLPInterfacePack/include/NLPFirstOrderInfo.h"
#include "SparseLinAlgPack/include/SpVectorClass.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "LinAlgPack/include/GenMatrixOut.h"
#include "LinAlgPack/include/VectorOut.h"
#include "LinAlgPack/include/VectorOp.h"
#include "LinAlgPack/include/PermOut.h"
#include "Misc/include/dynamic_cast_verbose.h"

namespace ReducedSpaceSQPPack {

using std::endl;
using std::setw;
using std::left;
using std::right;
using std::setprecision;

rSQPTrackStatsStd::rSQPTrackStatsStd(
	std::ostream& o, std::ostream& journal_out
	)
	: rSQPTrack(journal_out)
{
	set_output_stream(o);
}

void rSQPTrackStatsStd::set_output_stream(std::ostream& o)
{
	o_ = &o;
	num_QN_updates_ = 0;
	timer_.reset();
	timer_.start();
}

void rSQPTrackStatsStd::output_iteration(const Algorithm& p_algo) const
{
	const rSQPAlgo &algo = rsqp_algo(p_algo);
	const rSQPState &s =algo.rsqp_state();

	// All we have to do here is to just to count the number of quasi-newton updates
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats_(s).updated_k(0) ? &quasi_newton_stats_(s).get_k(0) : NULL );
	if( quasi_newt_stats ) {
		QuasiNewtonStats::EUpdate updated = quasi_newt_stats->updated();
		if( updated == QuasiNewtonStats::DAMPENED_UPDATED || updated == QuasiNewtonStats::UPDATED )
			num_QN_updates_++;
	}
}

void rSQPTrackStatsStd::output_final( const Algorithm& p_algo
	, EAlgoReturn algo_return ) const
{
	using DynamicCastHelperPack::const_dyn_cast;

	const rSQPAlgo			&algo = rsqp_algo(p_algo);
	const rSQPState			&s =algo.rsqp_state();
#ifdef _WINDOWS
	const NLPFirstOrderInfo
		&nlp = dynamic_cast<const NLPFirstOrderInfo&>(algo.nlp()); 
#else
	const NLPFirstOrderInfo
		&nlp = const_dyn_cast<NLPFirstOrderInfo>(algo.nlp()); 
#endif

	// Stop the timer
	timer_.stop();

	// Formating info
	const int
		p      = 18,
		stat_w = 15,
		val_w  = p + 10;

	// Get a Quasi-Newton statistics.
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats_(s).updated_k(0) ? &quasi_newton_stats_(s).get_k(0) : 0 );
	if( quasi_newt_stats ) {
		QuasiNewtonStats::EUpdate updated = quasi_newt_stats->updated();
		if( updated == QuasiNewtonStats::DAMPENED_UPDATED || updated == QuasiNewtonStats::UPDATED )
			num_QN_updates_++;
	}

	// status
	o() << left << setw(stat_w) << "status" << "= "
		<< right << setw(val_w);
	switch( algo_return ) {
		case GeneralIterationPack::TERMINATE_TRUE:
			o() << "solved";
			break;
		case GeneralIterationPack::TERMINATE_FALSE:
			o() << "except";
			break;
		case GeneralIterationPack::MAX_ITER_EXCEEDED:
			o() << "max_iter";
			break;
		case GeneralIterationPack::MAX_RUN_TIME_EXCEEDED:
			o() << "max_run_time";
			break;
		default:
			assert(0);
	}
	o() << "; # solved, except, max_iter, max_run_time\n";
	// niter
	o() << left << setw(stat_w) << "niter" << "= "
		<< right << setw(val_w) << s.k()
		<< "; # Number of rSQP iterations (plus 1?)\n";
	// nfunc
	o() << left << setw(stat_w) << "nfunc" << "= "
		<< right << setw(val_w) << std::_MAX(nlp.num_f_evals(),nlp.num_c_evals())
		<< "; # max( number f(x) evals, number c(x) evals )\n";
	// ngrad
	o() << left << setw(stat_w) << "ngrad" << "= "
		<< right << setw(val_w) << std::_MAX(nlp.num_Gf_evals(),nlp.num_Gc_evals())
		<< "; # max( number Gf(x) evals, number Gc(x) evals )\n";
	// CPU
	o() << left << setw(stat_w) << "CPU" << "= "
		<< right << setw(val_w) << timer_.read()
		<< "; # Number of CPU seconds total\n";
	// obj_func
	o() << left << setw(stat_w) << "obj_func" << "= "
		<< right << setw(val_w);
	if(s.f().updated_k(0))
		o() << s.f().get_k(0);
	else
		o() << "-";
	o() << "; # Objective function value f(x) at final point\n";
	// feas_kkt_err
	o() << left << setw(stat_w) << "feas_kkt_err" << "= "
		<< right << setw(val_w);
	if(s.feas_kkt_err().updated_k(0))
		o() << s.feas_kkt_err().get_k(0);
	else if(s.c().updated_k(0))
		o() << s.c().get_k(0).norm_inf();
	else
		o() << "-";
	o() << "; # Feasibility error at final point (scaled ||c(x)||inf, feas_err_k)\n";
	// opt_kkt_err
	o() << left << setw(stat_w) << "opt_kkt_err" << "= "
		<< right << setw(val_w);
	if(s.opt_kkt_err().updated_k(0))
		o() << s.opt_kkt_err().get_k(0);
	else if(s.rGL().updated_k(0))
		o() << s.rGL().get_k(0).norm_inf();
	else if(s.rGL().updated_k(-1))
		o() << s.rGL().get_k(-1).norm_inf();
	else
		o() << "-";
	o() << "; # Optimality error at final point (scaled ||rGL||inf, opt_err_k)\n";
	// nact
	o() << left << setw(stat_w) << "nact" << "= "
		<< right << setw(val_w);
	if(s.nu().updated_k(0))
		o() << s.nu().get_k(0).nz();
	else if(s.nu().updated_k(-1))
		o() << s.nu().get_k(-1).nz();
	else
		o() << "-";
	o() << "; # Number of total active constraints at the final point\n";
	// nbasis_change
	o() << left << setw(stat_w) << "nbasis_change" << "= "
		<< right << setw(val_w) << s.num_basis()
		<< "; # Number of basis changes\n";
	// nquasi_newton
	o() << left << setw(stat_w) << "nquasi_newton" << "= "
		<< right << setw(val_w) << num_QN_updates_
		<< "; # Number of quasi-newton updates\n";

}

} // end namespace ReducedSpaceSQPPack
