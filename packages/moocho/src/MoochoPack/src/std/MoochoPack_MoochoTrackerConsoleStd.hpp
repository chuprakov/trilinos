// ////////////////////////////////////////////////////////////////////////////
// rSQPTrackConsoleStd.h

#ifndef RSQP_TRACK_CONSOLE_STD_H
#define RSQP_TRACK_CONSOLE_STD_H

#include "ReducedSpaceSQPPack/include/rSQPTrack.h"
#include "Misc/include/stpwatch.h"

namespace ReducedSpaceSQPPack {

///
/** This rSQP iteration class provides a tablular output suitable for
  * an 80 char wide console.
  * 
  * Specifically, these object produces a table with the following fields:
  * \begin{itemize}
  * \item k : The SQP iteration counter starting at zero generally.
  * \item f : (iteration quantity f_k)
  * 			The value of the objective function at the current iteration.
  * 			This value may be scaled and the scaling factor will be
  * 			printed out before the table is produced.
  * \item ||c||s : (iteration quantity feas_kkt_err_k)
  * 			The scaled value of the constraints norm (usually the infinity
  * 			norm) at the current iteration.  This is the value compared
  * 			to the convergence criteria feas_tol (see the step
  * 			\Ref{CheckConvergenceStd_AddedStep}).
  * \item ||rGL||s : (iteration quantity opt_kkt_err_k)
  * 			The scaled value of the norm (usually the infinity
  * 			norm) of the reduced gradient of the Lagrangian
  * 			at the current iteration.  This is the value compared
  * 			to the convergence criteria opt_tol (see the step
  * 			\Ref{CheckConvergenceStd_AddedStep}).
  * \item QN : (iteration quantity quasi_newton_stats_k)
  * 			Information about the quasi-Newton update of the reduced
  * 			Hessian rHL_k.
  * 			\begin{description}
  * 			\item[IN] : rHL_k was reinitialized (identity?)
  * 			\item[UP] : A standard quasi-Newton update was performed on
  * 				rHL_km1 -> rHL_k
  * 			\item[DU] : A dampened quasi-Newton update (BFGS) was performed
  * 				on rHL_km1 -> rHL_k.
  * 			\item[SK] : The quasi-Newton update was skipped because the
  * 				current iterate was in the wrong region.
  * 			\item[IS] : The quasi-Newton update (BFGS) was skipped because
  * 				it was not positive definite or illdefined.
  * 			\end{description}	
  *	\item #act : (iteration quantity nu_k.nz())
  *				The number of active variable bounds at the current iteration.
  *	\item ||Ypy||2 : (iteration quantity Ypy_k)
  *				The 2 norm of the Range space (feasibility) contribution to the
  *				full step d.
  *	\item ||Zpz||2 : (iteration quantity Zpz_k)
  *				The 2 norm of the Null space (optimality) contribution to the
  *				full step d.
  *	\item ||d||inf : (iteration quantity d_k)
  *				The infinity norm of the step vector for the primal unknows
  *				x_kp1 = x_k + alpha_k * d_k.
  * \end{itemize}
  * 
  * The above quantities can tell you a lot about the progress of the SQP
  * algorithm.
  * 
  * ToDo: Finish discussion.
  * 
  * After the algorithm is finished the total solution time will be printed as
  * well as if the solution was found or not.  Also, the number of function and
  * gradient evaluations will be printed.  Note that the timer is started from
  * the moment this object is created or when set_output_stream(...) is called.
  */
class rSQPTrackConsoleStd : public rSQPTrack {
public:

	/// Construct with an output stream (console presumably)
	rSQPTrackConsoleStd(std::ostream& o, std::ostream& journal_out);

	/// Set the output stream for console outputting and restart the timer.
	void set_output_stream(std::ostream& o);

	// /////////////////////////////////////////////////////////
	// Overridden from AlgorithmTrack

	///
	void output_iteration(const Algorithm& algo) const;

	///
	void output_final(const Algorithm& algo, EAlgoReturn algo_return) const;

protected:

	/// Print the top header to the output
	void print_top_header(const rSQPState &s, const rSQPAlgo& algo) const;

	/// Print the header to the output
	void print_header(const rSQPState &s, const rSQPAlgo& algo) const;

	std::ostream& o() const
	{	return *const_cast<rSQPTrackConsoleStd*>(this)->o_; }

private:

	// ///////////////////////////////////////////
	// Private types

	enum { NUM_PRINT_LINES = 10 };
	
	// ///////////////////////////////////////////
	// Private data members

	std::ostream*						o_;
	mutable StopWatchPack::stopwatch	timer_;
	mutable int							printed_lines_;

	// Static formating info.
	static int		w_i2_;
	static char		ul_i2_[];
	static int		w_i4_;
	static char		ul_i4_[];
	static int		p2_;
	static int		w_p2_;
	static char		ul_p2_[];
	static int		p3_;
	static int		w_p3_;
	static char		ul_p3_[];

	// ///////////////////////////////////////////
	// Private member funcitons

	// Not defined and not to be called
	rSQPTrackConsoleStd();
};	// end class rSQPTrackConsoleStd

}	// end namespace ReducedSpaceSQPPack 

#endif	// RSQP_TRACK_CONSOLE_STD_H