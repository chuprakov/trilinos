// ////////////////////////////////////////////////////////////////////////////
// rSQPTrackConsoleStd.cpp

// disable VC 5.0 warnings about debugger limitations
#pragma warning(disable : 4786)	

#include <assert.h>
#include <iomanip>

#include "../../include/std/rSQPTrackConsoleStd.h"
#include "../../include/rSQPState.h"
#include "../../include/rsqp_algo_conversion.h"
#include "../../include/std/quasi_newton_stats.h"
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

// Static members
int		rSQPTrackConsoleStd::w_i2_		= 2;
char	rSQPTrackConsoleStd::ul_i2_[]	= "--";
int		rSQPTrackConsoleStd::w_i4_		= 4;
char	rSQPTrackConsoleStd::ul_i4_[]	= "----";
int		rSQPTrackConsoleStd::p2_		= 2;
int		rSQPTrackConsoleStd::w_p2_		= 8;
char	rSQPTrackConsoleStd::ul_p2_[]	= "--------";
int		rSQPTrackConsoleStd::p3_		= 3;
int		rSQPTrackConsoleStd::w_p3_		= 9;
char	rSQPTrackConsoleStd::ul_p3_[]	= "---------";

rSQPTrackConsoleStd::rSQPTrackConsoleStd(std::ostream& o, std::ostream& journal_out)
	: rSQPTrack(journal_out), o_(&o), printed_lines_(NUM_PRINT_LINES)
{
	timer_.start();
}

void rSQPTrackConsoleStd::set_output_stream(std::ostream& o)
{
	o_ = &o;
	timer_.reset();
	timer_.start();
}

void rSQPTrackConsoleStd::output_iteration(const Algorithm& p_algo) const
{
	const rSQPAlgo &algo = rsqp_algo(p_algo);
	const rSQPState &s =algo.rsqp_state();

	if(s.k() == 0) {
		print_top_header(s,algo);
		printed_lines_ = NUM_PRINT_LINES;
	}
	
	// Output the table's header
	if(printed_lines_ == NUM_PRINT_LINES) {
		printed_lines_ = 0;
		print_header(s,algo);
	}

	// ///////////////////////////////
	// Output a row for the iteration
	
	// Get a Quasi-Newton statistics.
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats(s).updated_k(0) ? &quasi_newton_stats(s).get_k(0) : 0 );

	// k
	o() << " " << right << setw(w_i4_) << s.k();
	// f
	if( s.f().updated_k(0) )
		o() << " " << setprecision(p3_) << right << setw(w_p3_) << s.f().get_k(0);
	else
		o() << " " << right << setw(w_p3_) << "-";
	// ||c||s
	if( s.feas_kkt_err().updated_k(0) )
		o() << " " << setprecision(p3_) << right << setw(w_p3_) << s.feas_kkt_err().get_k(0);
	else
		o() << " " << right << setw(w_p3_) << "-";
	// ||rGL||s
	if( s.opt_kkt_err().updated_k(0) )
		o() << " " << setprecision(p3_) << right << setw(w_p3_) << s.opt_kkt_err().get_k(0);
	else
		o() << " " << right << setw(w_p3_) << "-";
	// QN
	if( quasi_newt_stats ) {
		o() << " " << right << setw(w_i2_);
		switch( quasi_newt_stats->updated() ) {
			case QuasiNewtonStats::UNKNOWN:
				o() << "-";
				break;
			case QuasiNewtonStats:: REINITIALIZED:
				o() << "IN";
				break;
			case QuasiNewtonStats::DAMPENED_UPDATED:
				o() << "DU";
				break;
			case QuasiNewtonStats::UPDATED:
				o() << "UP";
				break;
			case QuasiNewtonStats::SKIPED:
				o() << "SK";
				break;
			case QuasiNewtonStats::INDEF_SKIPED:
				o() << "IS";
				break;
			default:
				assert(0);
		}
	}
	else {
		o() << " " << right << setw(w_i2_) << "-";
	}
	// #act
	if( s.nu().updated_k(0) )
		o() << " " << right << setw(w_i4_) << s.nu().get_k(0).nz();
	else
		o()	<< " " << right << setw(w_i4_) << "-";
	// ||Ypy||2
	if( s.Ypy().updated_k(0) )
		o() << " "<< setprecision(p2_)  << right << setw(w_p2_) << s.Ypy().get_k(0).norm_2();
	else
		o() << " " << right << setw(w_p2_) << "-";
	// ||Zpz||2
	if( s.Zpz().updated_k(0) )
		o() << " " << setprecision(p2_) << right << setw(w_p2_) << s.Zpz().get_k(0).norm_2();
	else
		o() << " " << right << setw(w_p2_) << "-";
	// ||d||inf
	if( s.d().updated_k(0) )
		o() << " " << setprecision(p2_) << right << setw(w_p2_) << s.d().get_k(0).norm_inf();
	else
		o() << " " << right << setw(w_p2_) << "-";
	// alpha
	if( s.alpha().updated_k(0) )
		o() << " " << setprecision(p2_) << right << setw(w_p2_) << s.alpha().get_k(0);
	else
		o() << " " << right << setw(w_p2_) << "-";

	o() << std::endl;

	++printed_lines_;
}

void rSQPTrackConsoleStd::output_final( const Algorithm& p_algo
	, EAlgoReturn algo_return ) const
{
	using DynamicCastHelperPack::const_dyn_cast;

	const rSQPAlgo			&algo = rsqp_algo(p_algo);
	const rSQPState			&s =algo.rsqp_state();
	const NLPFirstOrderInfo
#ifdef _WINDOWS
		&nlp = dynamic_cast<const NLPFirstOrderInfo&>(algo.nlp()); 
#else
		&nlp = const_dyn_cast<NLPFirstOrderInfo>(algo.nlp()); 
#endif

	// Output the table's header for the first iteration
	if(s.k() == 0) {
		print_top_header(s,algo);
		print_header(s,algo);
	}
	else {
		o()
			<< " " << right << ul_i4_ 		// "k"
			<< " " << right << ul_p3_		// "f"
			<< " " << right << ul_p3_		// "||c||s"
			<< " " << right << ul_p3_		// "||rGL||s"
			<< " " << right << ul_i2_		// "QN"
			<< " " << right << ul_i4_		// "#act"
			<< endl;
	}

	// //////////////////////////////////////////
	// Output a row for the final iteration
	
	// Get a Quasi-Newton statistics.
	const QuasiNewtonStats	*quasi_newt_stats =
		( quasi_newton_stats(s).updated_k(0) ? &quasi_newton_stats(s).get_k(0) : 0 );

	// k
	o() << " " << right << setw(w_i4_) << s.k();
	// f
	if( s.f().updated_k(0) )
		o() << " " << setprecision(p3_) << right << setw(w_p3_) << s.f().get_k(0);
	else
		o() << " " << right << setw(w_p3_) << "-";
	// ||c||s
	if( s.feas_kkt_err().updated_k(0) )
		o() << " " << setprecision(p3_) << right << setw(w_p3_) << s.feas_kkt_err().get_k(0);
	else
		o() << " " << right << setw(w_p3_) << "-";
	// ||rGL||s
	if( s.opt_kkt_err().updated_k(0) )
		o() << " " << setprecision(p3_) << right << setw(w_p3_) << s.opt_kkt_err().get_k(0);
	else
		o() << " " << right << setw(w_p3_) << "-";
	// QN
	if( quasi_newt_stats ) {
		o() << " " << right << setw(w_i2_);
		switch( quasi_newt_stats->updated() ) {
			case QuasiNewtonStats::UNKNOWN:
				o() << "-";
				break;
			case QuasiNewtonStats:: REINITIALIZED:
				o() << "IN";
				break;
			case QuasiNewtonStats::DAMPENED_UPDATED:
				o() << "DU";
				break;
			case QuasiNewtonStats::UPDATED:
				o() << "UP";
				break;
			case QuasiNewtonStats::SKIPED:
				o() << "SK";
				break;
			case QuasiNewtonStats::INDEF_SKIPED:
				o() << "IS";
				break;
			default:
				assert(0);
		}
	}
	else {
		o() << " " << right << setw(w_i2_) << "-";
	}
	// #act
	if( s.nu().updated_k(0) )
		o() << " " << right << setw(w_i4_) << s.nu().get_k(0).nz();
	else
		o()	<< " " << right << setw(w_i4_) << "-";

	o() << endl;

	// Print total time
	o() << "\nTotal time = " << timer_.read() << " sec\n";

	switch( algo_return ) {
		case GeneralIterationPack::TERMINATE_TRUE:
			o() << "\nJackpot! You have found the solution!!!!!!\n";
			break;
		case GeneralIterationPack::TERMINATE_FALSE:
			o() << "\nOops!  Not the solution.  Some error has occured!\n";
			break;
		case GeneralIterationPack::MAX_ITER_EXCEEDED:
			o() << "\nOops!  Not the solution.  Maximum number of SQP iteration exceeded!\n";
			break;
		case GeneralIterationPack::MAX_RUN_TIME_EXCEEDED:
			o() << "\nOops!  Not the solution.  Maximum runtime exceeded!\n";
			break;
		default:
			assert(0);
	}

	o()	<< "\nNumber of function evaluations:\n"
		<<     "-------------------------------\n"
		<< "f(x)  : " << nlp.num_f_evals() << endl
		<< "c(x)  : " << nlp.num_c_evals() << endl
		<< "Gf(x) : " << nlp.num_Gf_evals() << endl
		<< "Gc(x) : " << nlp.num_Gc_evals() << endl;
}

void rSQPTrackConsoleStd::print_top_header(const rSQPState &s
	, const rSQPAlgo &algo) const
{
	o()	<< "\n\n********************************\n"
		<< "*** Start of rSQP Iterations ***\n"
		<< "n = " << s.x().get_k(0).cv().size()
		<< ", m = " << s.c().get_k(0).cv().size()
		<< ", nz = ";
	if( s.Gc().updated_k(0) )
		o()	<< s.Gc().get_k(0).nz() << endl;
	else
		o()	<< "?\n";
	if( algo.nlp().scale_f() != 1.0 ) {
		o()	<< "f(x) is scaled by : " << algo.nlp().scale_f() << endl;
	}
}

void rSQPTrackConsoleStd::print_header(const rSQPState &s
	, const rSQPAlgo &algo) const
{
	o()
		<< endl
		<< " " << left << setw(w_i4_) << "k"
		<< " " << left << setw(w_p3_) << "f"
		<< " " << left << setw(w_p3_) << "||c||s"
		<< " " << left << setw(w_p3_) << "||rGL||s"
		<< " " << left << setw(w_i2_) << "QN"
		<< " " << left << setw(w_i4_) << "#act"
		<< " " << left << setw(w_p2_) << "||Ypy||2"
		<< " " << left << setw(w_p2_) << "||Zpz||2"
		<< " " << left << setw(w_p2_) << "||d||inf"
		<< " " << left << setw(w_p2_) << "alpha"
		<< endl
		<< " " << right << ul_i4_ 		// "k"
		<< " " << right << ul_p3_		// "f"
		<< " " << right << ul_p3_		// "||c||s"
		<< " " << right << ul_p3_		// "||rGL||s"
		<< " " << right << ul_i2_		// "QN"
		<< " " << right << ul_i4_		// "#act"
		<< " " << right << ul_p2_		// "||Ypy||2"
		<< " " << right << ul_p2_		// "||Zpz||2"
		<< " " << right << ul_p2_		// "||d||inf"
		<< " " << right << ul_p2_		// "alpha"
		<< endl;
}

} // end namespace ReducedSpaceSQPPack
