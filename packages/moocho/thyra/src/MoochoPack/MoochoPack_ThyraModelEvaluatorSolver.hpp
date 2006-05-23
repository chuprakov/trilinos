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

#ifndef MOOCHOPACK_THYRA_MODEL_EVALUATOR_SOLVER_HPP
#define MOOCHOPACK_THYRA_MODEL_EVALUATOR_SOLVER_HPP

#include "MoochoPack_MoochoSolver.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Thyra_DefaultNominalBoundsOverrideModelEvaluator.hpp"
#include "Thyra_DefaultFinalPointCaptureModelEvaluator.hpp"

namespace MoochoPack {

/** \brief NLP Solver class for models represented through
 * <tt>Thyra::ModelEvaluator</tt>.
 *
 * ToDo: Finish documetation!
 */
class ThyraModelEvaluatorSolver {
public:
  
  /** \brief . */
  ThyraModelEvaluatorSolver();
  
  /** \brief . */
  void setupCLP(
    Teuchos::CommandLineProcessor *clp
    );
  
  /** \brief . */
  MoochoSolver& getSolver();
  
  /** \brief . */
  const MoochoSolver& getSolver() const;
  
  /** \brief . */
  void setModel(
    const Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> > &model
    ,const int                                                     p_idx  = 0
    ,const int                                                     g_idx  = 0
    );
    
  /** \brief . */
  void readInitialGuess(
    std::ostream *out = NULL
    );
  
  /** \brief . */
  MoochoSolver::ESolutionStatus	solve();
  
  /** \brief . */
  void writeFinalSolution(
    std::ostream *out = NULL
    ) const;


private:

  MoochoSolver                                               solver_;
  
  Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> >   origModel_;
  int                                                        p_idx_;
  int                                                        g_idx_;
  
  Teuchos::RefCountPtr<Thyra::DefaultNominalBoundsOverrideModelEvaluator<value_type> > nominalModel_;
  
  Teuchos::RefCountPtr<Thyra::DefaultFinalPointCaptureModelEvaluator<value_type> > finalPointModel_;

  Teuchos::RefCountPtr<Thyra::ModelEvaluator<value_type> > outerModel_;
  
  bool                do_sim_;
  bool                use_direct_;
  bool                use_black_box_;
  bool                use_finite_diff_;
  Thyra::DirectionalFiniteDiffCalculatorTypes::EFDMethodType fd_method_type_;
  double              fd_step_len_;
  double              fwd_newton_tol_;
  int                 fwd_newton_max_iters_;
  std::string         stateGuessFileBase_;
  double              scaleStateGuess_;
  std::string         paramGuessFileBase_;
  double              scaleParamGuess_;
  std::string         stateSoluFileBase_;
  std::string         paramSoluFileBase_;

};

} // namespace MoochoPack

#endif	// MOOCHOPACK_THYRA_MODEL_EVALUATOR_SOLVER_HPP
