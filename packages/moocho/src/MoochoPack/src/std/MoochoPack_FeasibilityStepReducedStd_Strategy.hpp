// ////////////////////////////////////////////////////////////////////////////////
// FeasibilityStepReducedStd_Strategy.h
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

#ifndef FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H
#define FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H

#include "FeasibilityStep_Strategy.h"
#include "QuasiRangeSpaceStep_Strategy.h"
#include "d_bounds_iter_quant.h"
#include "ConstrainedOptimizationPack/include/QPSolverRelaxed.h"
#include "ConstrainedOptimizationPack/include/QPSolverRelaxedTester.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "Misc/include/StandardCompositionMacros.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ReducedSpaceSQPPack {

///
/** Implements the feasibility step computation for reduced space SQP.
 */
class FeasibilityStepReducedStd_Strategy : public FeasibilityStep_Strategy {
public:

	/// <<std comp>> members for the qp solver
	STANDARD_COMPOSITION_MEMBERS( QuasiRangeSpaceStep_Strategy, quasi_range_space_step )

	typedef ConstrainedOptimizationPack::QPSolverRelaxedTester
		QPSolverRelaxedTester;

	/// <<std comp>> members for the qp solver
	STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxed, qp_solver )

	/// <<std comp>> members for comparision object compatible with Gc
	STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxedTester, qp_tester )
		
	///
	enum EQPObjective {
		OBJ_MIN_FULL_STEP           // min 1/2 * (Y*wy + Z*wz)'*(Y*wy + Z*wz)
		,OBJ_MIN_NULL_SPACE_STEP    // min 1/2 * wz'*wz
		,OBJ_RSQP                   // min qp_grad_k'*wz + 1/2 * wz'*rHL_k*wz
	};

	///
	/** Set what is used for the QP objective.
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPObjective, qp_objective )

	///
	enum EQPTesting { QP_TEST_DEFAULT, QP_TEST, QP_NO_TEST };

	///
	/** Set how and if the QP solution is tested.
	  *
	  * ToDo: Finish documentation.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPTesting, qp_testing )

	/// Construct and initialize
	FeasibilityStepReducedStd_Strategy(
		const quasi_range_space_step_ptr_t   &quasi_range_space_step
		,const qp_solver_ptr_t               &qp_solver
		,const qp_tester_ptr_t               &qp_tester
		,EQPObjective                        qp_objective     = OBJ_MIN_NULL_SPACE_STEP
		,EQPTesting                          qp_testing       = QP_TEST_DEFAULT
		);

	// ////////////////////////////////////////////
	// Overridden from FeasibilityStep_Strategy

	///
	/** Computes a feasibility step by computing simple range and null space components.
	 *
	 * ToDo: Finish documentation!
	 *
	 */
 	bool compute_feasibility_step(
		std::ostream& out, EJournalOutputLevel olevel, rSQPAlgo *algo, rSQPState *s
		,const VectorSlice& xo, const VectorSlice& c_xo, VectorSlice* w
	  	);

	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;

private:
	//
	typedef ReferenceCountingPack::ref_count_ptr<const MatrixWithOp> Hess_ptr_t;
	//
	d_bounds_iq_member			d_bounds_;
	int                         current_k_;
	Hess_ptr_t                  Hess_ptr_;
	Vector                      grad_store_;
	GenMatrix                   Hess_store_;

}; // end class FeasibilityStepReducedStd_Strategy

} // end namespace ReducedSpaceSQPPack

#endif // FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H
