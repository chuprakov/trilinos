// //////////////////////////////////////////////////////////
// QPSolverRelaxedQPOPT.h

#ifndef QP_SOLVER_RELAXED_QPOPT_H
#define QP_SOLVER_RELAXED_QPOPT_H

#include "QPSolverRelaxedQPOPTSOL.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Node base clase for the primal QP solvers QPOPT and QPSOL.
  *
  * In this implementation it is required that G only support the \Ref{MatrixWithOp}
  * interface and is therefore quite flexible in the QPs it can solve.
  */
class QPSolverRelaxedQPOPT : public QPSolverRelaxedQPOPTSOL
{
public:

	///
	typedef QPSolverRelaxedQPOPTSOL inherited;

	///
	/** Set the maximum number of QP iterations as max_itr = max_qp_iter_frac * n.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac )

	///
	QPSolverRelaxedQPOPT(
		  value_type        max_qp_iter_frac	= 10.0
		);

	///
	~QPSolverRelaxedQPOPT();

	// /////////////////////////////////
	// Overridden from QPSolverRelaxed

	///
	void release_memory();

protected:

	// /////////////////////////////////////////////////////////////
	// Overridden from QPSolverRelaxedQPOPTSOL

	///
	f_int liwork(f_int N, f_int NCLIN) const;
	///
	f_int lrwork(f_int N, f_int NCLIN) const;
	///
	EInform call_qp_solver(bool warm_start);

private:

	// ////////////////////////////
	// Private types

	///
	enum EQPOPTInform {
		STRONG_LOCAL_MIN			= 0,
		WEAK_LOCAL_MIN				= 1,
		UNBOUNDED					= 2,
		INFEASIBLE					= 3,
		ITMAX_EXCEEDED				= 4,
		MAX_DOF_TOO_SMALL			= 5,
		INVALID_INPUT				= 6,
		PROB_TYPE_NOT_REGOG			= 7
	};

	// ////////////////////////////
	// Private data members

	// extra QPOPT control and input parameters.

	// control

	f_int						ITMAX_;
	f_dbl_prec					BIGBND_;
	f_dbl_prec					FEATOL_;

	// input/output

	f_int						LDA_;
	f_int						LDH_;
	f_dbl_prec*					H_;
	f_int						INFORM_;

	// ////////////////////////////
	// Private member functions

	// not defined and not to be called.
	QPSolverRelaxedQPOPT(const QPSolverRelaxedQPOPT&);
	QPSolverRelaxedQPOPT& operator=(const QPSolverRelaxedQPOPT&);

};	// end class QPSolverRelaxedQPOPT

}	// end namespace ConstrainedOptimizationPack

#endif // QP_SOLVER_RELAXED_QPOPT_H
