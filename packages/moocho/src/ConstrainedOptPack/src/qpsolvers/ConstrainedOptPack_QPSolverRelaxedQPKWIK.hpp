// ///////////////////////////////////////////////////////////////////////////////////////
// QPSolverRelaxedQPKWIK.h
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

#ifndef QP_SOLVER_RELAXED_QPKWIK_H
#define QP_SOLVER_RELAXED_QPKWIK_H

#include <vector>

#include "QPSolverRelaxed.h"
#include "LinAlgPack/include/VectorClass.h"
#include "LinAlgPack/include/GenMatrixClass.h"
#include "Misc/include/StandardMemberCompositionMacros.h"

namespace ConstrainedOptimizationPack {

///
/** Solves Quadratic Programming (QP) problem using the primal-dual active-set
  * solver QPKWIK.
  *
  * In this implementation it is required that G support the \Ref{MatrixExtractInvCholFactor}
  * interface and is therefore quite restrictive on the likes of QPs it can solve.
  */
class QPSolverRelaxedQPKWIK : public QPSolverRelaxed
{
public:

	///
	/** Set the maximum number of QP iterations as max_itr = max_qp_iter_frac * n.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac )

	///
	/** Set the value of an infinite bound.
	  */
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, infinite_bound )

	///
	QPSolverRelaxedQPKWIK(
		  value_type        max_qp_iter_frac	= 10.0
		  ,value_type       infinite_bound      = 1e+20
		);

	///
	~QPSolverRelaxedQPKWIK();

	// /////////////////////////////////
	// Overridden from QPSolverRelaxed

	///
	QPSolverStats get_qp_stats() const;

	///
	void release_memory();

protected:

	// /////////////////////////////////
	// Overridden from QPSolverRelaxed

	///
	QPSolverStats::ESolutionType imp_solve_qp(
		  std::ostream* out, EOutputLevel olevel, ERunTests test_what
		, const VectorSlice& g, const MatrixWithOp& G
		, value_type etaL
		, const SpVectorSlice& dL, const SpVectorSlice& dU
		, const MatrixWithOp* E, BLAS_Cpp::Transp trans_E, const VectorSlice* b
			, const SpVectorSlice* eL, const SpVectorSlice* eU
		, const MatrixWithOp* F, BLAS_Cpp::Transp trans_F, const VectorSlice* f
		, value_type* obj_d
		, value_type* eta, VectorSlice* d
		, SpVector* nu
		, SpVector* mu, VectorSlice* Ed
		, VectorSlice* lambda, VectorSlice* Fd
		);

private:

	// //////////////////////////////////////////////////////////////
	// Private types

	///
	typedef std::vector<f_int>			IBND_t;
	///
	typedef std::vector<f_int>			IACTSTORE_t;
	///
	typedef std::vector<f_int>			IACT_t;
	///
	typedef std::vector<f_int>			ISTATE_t;

	// //////////////////////////////////////////////////////////////
	// Private Data Members.

	QPSolverStats			qp_stats_;

	// Inverse mapping for IBND_INV(j) == k <-> IBND(k) == j
	IBND_t  IBND_INV_;

	// Parameters to QPKWIK

	///
	f_int		N_;
	///
	f_int		M1_;
	///
	f_int		M2_;
	///
	f_int		M3_;
	///
	Vector		GRAD_;
	///
	GenMatrix	UINV_AUG_;
	///
	f_int		LDUINV_AUG_;
	///
	IBND_t		IBND_;
	///
	Vector		BL_;
	///
	Vector		BU_;
	///
	GenMatrix	A_;
	///
	f_int		LDA_;
	///
	Vector		YPY_;
	///
	f_int		IYPY_;
	///
	f_int		WARM_;
	///
	f_dbl_prec	NUMPARAM_[3];
	///
	f_int       MAX_ITER_;

	// Input / Output

	///
	Vector		X_;
	///
	f_int		NACTSTORE_;
	///
	IACTSTORE_t	IACTSTORE_;
	///
	f_int		INF_;
	
	// Output

	///
	f_int		NACT_;
	///
	IACT_t		IACT_;
	///
	Vector		UR_;
	///
	f_dbl_prec	EXTRA_;
	///
	f_int		ITER_;
	///
	f_int		NUM_ADDS_;
	///
	f_int		NUM_DROPS_;
	
	// Internal state

	///
	ISTATE_t	ISTATE_;
	

	// Workspace

	///
	f_int		LRW_;
	///
	Vector		RW_;

	// ////////////////////////////
	// Private member functions

};	// end class QPSolverRelaxedQPKWIK

}	// end namespace ConstrainedOptimizationPackTypes

#endif	// QP_SOLVER_RELAXED_QPKWIK_H
