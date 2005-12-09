// ///////////////////////////////////////////////////////////////////////////////////////
// ConstrainedOptPack_QPSolverRelaxedQPKWIK.hpp
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

#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

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

	/** @name Initialization */
	//@{

	/// Set the maximum number of QP iterations as max_itr = max_qp_iter_frac * n.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac )

	/// Set the value of an infinite bound.
	STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, infinite_bound )

	///
	QPSolverRelaxedQPKWIK(
		  value_type        max_qp_iter_frac	= 10.0
		  ,value_type       infinite_bound      = 1e+20
		);

	///
	~QPSolverRelaxedQPKWIK();

	//@}

	/** @name Overridden from QPSolverRelaxed */
	//@{

	///
	QPSolverStats get_qp_stats() const;
	///
	void release_memory();

	//@}

protected:

	/** @name Overridden from QPSolverRelaxed */
	//@{

	///
	QPSolverStats::ESolutionType imp_solve_qp(
		std::ostream* out, EOutputLevel olevel, ERunTests test_what
		,const Vector& g, const MatrixSymOp& G
		,value_type etaL
		,const Vector* dL, const Vector* dU
		,const MatrixOp* E, BLAS_Cpp::Transp trans_E, const Vector* b
		,const Vector* eL, const Vector* eU
		,const MatrixOp* F, BLAS_Cpp::Transp trans_F, const Vector* f
		,value_type* obj_d
		,value_type* eta, VectorMutable* d
		,VectorMutable* nu
		,VectorMutable* mu, VectorMutable* Ed
		,VectorMutable* lambda, VectorMutable* Fd
		);

	//@}

private:

	// //////////////////////////////////////////////////////////////
	// Private types

	///
	typedef std::vector<index_type>  IBND_t;
	///
	typedef std::vector<index_type>  IACTSTORE_t;
	///
	typedef std::vector<index_type>  IACT_t;
	///
	typedef std::vector<index_type>  ISTATE_t;

	// //////////////////////////////////////////////////////////////
	// Private Data Members.

	QPSolverStats   qp_stats_;

	// Inverse mapping for IBND_INV(j) == k <-> IBND(k) == j
	IBND_t          IBND_INV_;

	// Parameters to QPKWIK

	///
	index_type      N_;
	///
	index_type      M1_;
	///
	index_type      M2_;
	///
	index_type      M3_;
	///
	DVector          GRAD_;
	///
	DMatrix       UINV_AUG_;
	///
	index_type      LDUINV_AUG_;
	///
	IBND_t          IBND_;
	///
	DVector          BL_;
	///
	DVector          BU_;
	///
	DMatrix       A_;
	///
	index_type		LDA_;
	///
	DVector          YPY_;
	///
	index_type      IYPY_;
	///
	index_type      WARM_;
	///
	value_type      NUMPARAM_[3];
	///
	index_type      MAX_ITER_;

	// Input / Output

	///
	DVector          X_;
	///
	index_type      NACTSTORE_;
	///
	IACTSTORE_t     IACTSTORE_;
	///
	index_type      INF_;
	
	// Output

	///
	index_type      NACT_;
	///
	IACT_t          IACT_;
	///
	DVector          UR_;
	///
	value_type      EXTRA_;
	///
	index_type      ITER_;
	///
	index_type      NUM_ADDS_;
	///
	index_type      NUM_DROPS_;
	
	// Internal state

	///
	ISTATE_t        ISTATE_;

	// Workspace

	///
	index_type      LRW_;
	///
	DVector          RW_;

}; // end class QPSolverRelaxedQPKWIK

} // end namespace ConstrainedOptimizationPackTypes

#endif // QP_SOLVER_RELAXED_QPKWIK_H
