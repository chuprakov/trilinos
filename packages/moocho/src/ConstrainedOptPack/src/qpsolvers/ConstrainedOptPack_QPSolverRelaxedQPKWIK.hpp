// @HEADER
// ***********************************************************************
// 
// Moocho: Multi-functional Object-Oriented arCHitecture for Optimization
//                  Copyright (2003) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef QP_SOLVER_RELAXED_QPKWIK_H
#define QP_SOLVER_RELAXED_QPKWIK_H

#include <vector>

#include "ConstrainedOptPack_QPSolverRelaxed.hpp"
#include "DenseLinAlgPack_DVectorClass.hpp"
#include "DenseLinAlgPack_DMatrixClass.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

namespace ConstrainedOptPack {

/** \brief Solves Quadratic Programming (QP) problem using the primal-dual active-set
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
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, max_qp_iter_frac );

  /// Set the value of an infinite bound.
  STANDARD_MEMBER_COMPOSITION_MEMBERS( value_type, infinite_bound );

  /** \brief . */
  QPSolverRelaxedQPKWIK(
      value_type        max_qp_iter_frac	= 10.0
      ,value_type       infinite_bound      = 1e+20
    );

  /** \brief . */
  ~QPSolverRelaxedQPKWIK();

  //@}

  /** @name Overridden from QPSolverRelaxed */
  //@{

  /** \brief . */
  QPSolverStats get_qp_stats() const;
  /** \brief . */
  void release_memory();

  //@}

protected:

  /** @name Overridden from QPSolverRelaxed */
  //@{

  /** \brief . */
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

  /** \brief . */
  typedef std::vector<index_type>  IBND_t;
  /** \brief . */
  typedef std::vector<index_type>  IACTSTORE_t;
  /** \brief . */
  typedef std::vector<index_type>  IACT_t;
  /** \brief . */
  typedef std::vector<index_type>  ISTATE_t;

  // //////////////////////////////////////////////////////////////
  // Private Data Members.

  QPSolverStats   qp_stats_;

  // Inverse mapping for IBND_INV(j) == k <-> IBND(k) == j
  IBND_t          IBND_INV_;

  // Parameters to QPKWIK

  /** \brief . */
  index_type      N_;
  /** \brief . */
  index_type      M1_;
  /** \brief . */
  index_type      M2_;
  /** \brief . */
  index_type      M3_;
  /** \brief . */
  DVector          GRAD_;
  /** \brief . */
  DMatrix       UINV_AUG_;
  /** \brief . */
  index_type      LDUINV_AUG_;
  /** \brief . */
  IBND_t          IBND_;
  /** \brief . */
  DVector          BL_;
  /** \brief . */
  DVector          BU_;
  /** \brief . */
  DMatrix       A_;
  /** \brief . */
  index_type		LDA_;
  /** \brief . */
  DVector          YPY_;
  /** \brief . */
  index_type      IYPY_;
  /** \brief . */
  index_type      WARM_;
  /** \brief . */
  value_type      NUMPARAM_[3];
  /** \brief . */
  index_type      MAX_ITER_;

  // Input / Output

  /** \brief . */
  DVector          X_;
  /** \brief . */
  index_type      NACTSTORE_;
  /** \brief . */
  IACTSTORE_t     IACTSTORE_;
  /** \brief . */
  index_type      INF_;
  
  // Output

  /** \brief . */
  index_type      NACT_;
  /** \brief . */
  IACT_t          IACT_;
  /** \brief . */
  DVector          UR_;
  /** \brief . */
  value_type      EXTRA_;
  /** \brief . */
  index_type      ITER_;
  /** \brief . */
  index_type      NUM_ADDS_;
  /** \brief . */
  index_type      NUM_DROPS_;
  
  // Internal state

  /** \brief . */
  ISTATE_t        ISTATE_;

  // Workspace

  /** \brief . */
  index_type      LRW_;
  /** \brief . */
  DVector          RW_;

}; // end class QPSolverRelaxedQPKWIK

} // end namespace ConstrainedOptimizationPackTypes

#endif // QP_SOLVER_RELAXED_QPKWIK_H
