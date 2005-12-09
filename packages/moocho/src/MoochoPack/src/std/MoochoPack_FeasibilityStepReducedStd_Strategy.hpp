// ////////////////////////////////////////////////////////////////////////////////
// MoochoPack_FeasibilityStepReducedStd_Strategy.hpp
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

#include "MoochoPack_FeasibilityStep_Strategy.hpp"
#include "MoochoPack_QuasiRangeSpaceStep_Strategy.hpp"
#include "MoochoPack_d_bounds_iter_quant.hpp"
#include "IterationPack_CastIQMember.hpp"
#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "ConstrainedOptPack_QPSolverRelaxedTester.hpp"
#include "AbstractLinAlgPack_MatrixOp.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "AbstractLinAlgPack_VectorSpace.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace MoochoPack {

///
/** Implements the feasibility step computation for reduced space SQP.
 */
class FeasibilityStepReducedStd_Strategy : public FeasibilityStep_Strategy
{
public:

	/// <<std comp>> members for the qp solver
	STANDARD_COMPOSITION_MEMBERS( QuasiRangeSpaceStep_Strategy, quasi_range_space_step )

	typedef ConstrainedOptPack::QPSolverRelaxedTester
		QPSolverRelaxedTester;

	/// QP solver
	STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxed, qp_solver )

	/// Comparision object compatible with Gc
	STANDARD_COMPOSITION_MEMBERS( QPSolverRelaxedTester, qp_tester )
		
	///
	enum EQPObjective {
		OBJ_MIN_FULL_STEP           ///< min 1/2 * (Y*wy + Z*wz)'*(Y*wy + Z*wz)
		,OBJ_MIN_NULL_SPACE_STEP    ///< min 1/2 * wz'*wz
		,OBJ_RSQP                   ///< min qp_grad_k'*wz + 1/2 * wz'*rHL_k*wz
	};

	///
	/** Set what is used for the QP objective.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EQPObjective, qp_objective )

	///
	enum EQPTesting {
		QP_TEST_DEFAULT     ///< Decide based on olevel input to <tt>compute_feasibility_step(...)</tt>
		,QP_TEST            ///< Perform the tests
		,QP_NO_TEST         ///< Don't perform the tests
	};

	///
	/** Set how and if the QP solution is tested.
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
	/** Computes a feasibility step by computing simple quasi-range and null space components.
	 *
	 * ToDo: Finish documentation!
	 *
	 */
 	bool compute_feasibility_step(
		std::ostream& out, EJournalOutputLevel olevel, NLPAlgo *algo, NLPAlgoState *s
		,const Vector& xo, const Vector& c_xo, VectorMutable* w
	  	);

	///
	void print_step( std::ostream& out, const std::string& leading_str ) const;

private:

	IterationPack::CastIQMember<VectorMutable>  dl_iq_;
	IterationPack::CastIQMember<VectorMutable>  du_iq_;
	int                                                      current_k_;
	Teuchos::RefCountPtr<const MatrixOp>            Hess_ptr_;
	VectorSpace::vec_mut_ptr_t                               grad_store_;
	DMatrix                                                Hess_store_;

}; // end class FeasibilityStepReducedStd_Strategy

} // end namespace MoochoPack

#endif // FEASIBILITY_STEP_REDUCED_STD_STRATEGY_H
