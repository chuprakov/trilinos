// @HEADER
// ***********************************************************************
// 
//               TSFCore: Trilinos Solver Framework Core
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

// ///////////////////////////////////////////////////
// Epetra_LinearSystemScaler.hpp

#ifndef EPETRA_LINEAR_SYSTEM_SCALER_HPP
#define EPETRA_LINEAR_SYSTEM_SCALER_HPP

#include "Epetra_ProductOperator.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"

class Epetra_Operator;
class Epetra_RowMatrix;
class Epetra_Vector;
class Epetra_MultiVector;

namespace Epetra {

///
/** A utility class for scaling Epetra linear systems.
 *
 * The purpose of this class is to make the use of a number of options
 * for forward and/or adjoint solves as transparent as possible.  This
 * class only depends Epetra and Teuchos and is therefore very
 * general.  Clients can select left and/or right scaling for the
 * forward and/or adjoint solves and can automatically scale
 * Epetra_Operator objects for the matrix operator and (the optional)
 * preconditoner and RHS and solution Epetra_MultiVector objects.
 *
 * This scaling class will not alter any of the Epetra_Operator
 * objects that are passed to it.  Only Epetra_MultiVector objects
 * will be modified in calls to <tt>transformRhs()</tt> and
 * <tt>transformSol()</tt>.
 *
 * Clients can also consider this scaling class to be stateless with
 * respect to any particular scaling or linear system.  This is
 * important when multiple linear systems must be maintained
 * simultaneously.  Therefore, this class can be considered to be a
 * type of factory object.
 *
 * ToDo: Describe uses cases this class supports!
 *
 * Note, the default copy constructor and assignment operators are
 * allowed since they have the correct semantics.  Copying (or
 * assigning) a <tt>LinearSystemScaler</tt> object means copying its
 * options.
 */
class LinearSystemScaler {
public:

	/** @name Public types */
	//@{

	///
	enum EFwdSolveLeftScaling { FWD_SOLVE_LEFT_SCALING_NONE, FWD_SOLVE_LEFT_SCALING_ROW_SUM };
	///
	enum EFwdSolveRightScaling { FWD_SOLVE_RIGHT_SCALING_NONE, FWD_SOLVE_RIGHT_SCALING_COL_SUM };
	///
	enum EAdjSolveLeftScaling { ADJ_SOLVE_LEFT_SCALING_NONE, ADJ_SOLVE_LEFT_SCALING_COL_SUM };
	///
	enum EAdjSolveRightScaling { ADJ_SOLVE_RIGHT_SCALING_NONE, ADJ_SOLVE_RIGHT_SCALING_ROW_SUM };
	///
	enum EFwdSolvePrec { FWD_SOLVE_PREC_DEFAULT, FWD_SOLVE_PREC_LEFT, FWD_SOLVE_PREC_RIGHT };
	///
	enum EAdjSolvePrec { ADJ_SOLVE_PREC_DEFAULT, ADJ_SOLVE_PREC_LEFT, ADJ_SOLVE_PREC_RIGHT };
	///
	class Scaling {
	public:
		///
		Scaling()
			{}
		///
		Scaling(
			const Teuchos::RefCountPtr<const Epetra_Vector>    &leftScaling
			,const Teuchos::RefCountPtr<const Epetra_Vector>   &rightScaling
			)
			:leftScaling_(leftScaling)
			,rightScaling_(rightScaling)
			{}
		///
		Teuchos::RefCountPtr<const Epetra_Vector> leftScaling() const { return leftScaling_; }
		///
		Teuchos::RefCountPtr<const Epetra_Vector> rightScaling() const { return rightScaling_; }
	private:
		Teuchos::RefCountPtr<const Epetra_Vector>  leftScaling_;
		Teuchos::RefCountPtr<const Epetra_Vector>  rightScaling_;
	}; // class Scaling

	//@}

	/** @name Constructors / initializers / accessors */
	//@{

  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EFwdSolveLeftScaling, fwdSolveLeftScaling );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EFwdSolveRightScaling, fwdSolveRightScaling );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EAdjSolveLeftScaling, adjSolveLeftScaling );
  ///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EAdjSolveRightScaling, adjSolveRightScaling );
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EFwdSolvePrec, fwdSolvePrec );
	///
	STANDARD_MEMBER_COMPOSITION_MEMBERS( EAdjSolvePrec, adjSolvePrec );

	///
	/** Construct uninitialized but with default option values.
	 */
 	LinearSystemScaler(
		const EFwdSolveLeftScaling    fwdSolveLeftScaling   = FWD_SOLVE_LEFT_SCALING_NONE
		,const EFwdSolveRightScaling  fwdSolveRightScaling  = FWD_SOLVE_RIGHT_SCALING_NONE
		,const EAdjSolveLeftScaling   adjSolveLeftScaling   = ADJ_SOLVE_LEFT_SCALING_NONE
		,const EAdjSolveRightScaling  adjSolveRightScaling  = ADJ_SOLVE_RIGHT_SCALING_NONE
		,const EFwdSolvePrec          fwdSolvePrec          = FWD_SOLVE_PREC_DEFAULT
		,const EAdjSolvePrec          adjSolvePrec          = ADJ_SOLVE_PREC_DEFAULT
		);

	//@}

	/** @name Scaling methods */
	//@{

	// ToDo: Add a method to set row and/or column scalings!

	///
	/** Compute automatic scaling factors given an Epetra_RowMatrix.
	 *
	 * This function computes row and column scalings vectors and stores
	 * them in the returned <tt>Scaling</tt> object.
	 */
	Teuchos::RefCountPtr<const Scaling> computeScaling(
		const Epetra_RowMatrix    &Op
		,const Teuchos::ETransp   Op_trans = Teuchos::NO_TRANS
		) const;

	///
	/** Apply scaling factors to an operator and its preconditioner for
	 * a forward or an adjoint solve.
	 *
	 * @param  scaling
	 *             [in] Row and column scalings.
	 * @param  scaling_trans
	 *             [in] Determine how the row and column scalings in 
	 *             <tt>scaling</tt> are to be viewed.
	 * @param  Op  [in] Operator to scale.  Note, this object's reference
	 *             will be embedded in <tt>Op</tt> but the only
	 *             way that <tt>*Op</tt> will be affected is in calls
	 *             to <tt>Epetra_Operator::SetUseTranspose()</tt>.
	 * @param  Op_trans
	 *             [in] Determines if <tt>Op</tt> is viewed as the transpose
	 *             or not.
	 * @param  Prec
	 *             [in] Preconditioner to scale.  Note, this object's reference
	 *             will be embedded in <tt>Prec</tt> but the only
	 *             way that <tt>*Op</tt> will be affected is in calls
	 *             to <tt>Epetra_Operator::SetUseTranspose()</tt>.
	 *             Note, it is allowed for <tt>Prec.get()==NULL</tt> if there
	 *             is no preconditioner.
	 * @param  Prec_trans
	 *             [in] Determines if <tt>Prec</tt> is viewed as the transpose
	 *             or not.
	 * @param  Prec_inverse
	 *             [in] Determines if the application of the preconditioner is
	 *             applied using <tt>Prec->Apply(...)</tt>
	 *             (<tt>Prec_inverse=ProductOperator::APPLY_MODE_APPLY</tt>)
	 *             or using <tt>Prec->ApplyInverse(...)</tt>
	 *             (<tt>Prec_inverse=ProductOperator::APPLY_MODE_APPLY_INVERSE</tt>).
	 * @param  Solve_trans
	 *             [in] Determines if operators are to be returned for the
	 *             forward (Teuchos::NO_TRANS) or the adjoint (!Teuchos::NO_TRANS) solve.
	 * @param  Op_solve
	 *             [out] A perhaps new <tt>Epetra_Operator</tt> object that
	 *             represents the scaled operator and perhaps the aggregate
	 *             preconditioner.
	 * @param  Prec_solve_side
	 *             [in] Determines if the external preconditioner is to be applied
	 *             on the left or on the right.
	 * @param  Prec_solve_inverse
	 *             [in] Determines if the application of the returned preconditioner is
	 *             applied using <tt>Prec_solve->->Apply(...)</tt>
	 *             (<tt>Prec_inverse=ProductOperator::APPLY_MODE_APPLY</tt>)
	 *             or using <tt>Prec_solve->->ApplyInverse(...)</tt>
	 *             (<tt>Prec_inverse=ProductOperator::APPLY_MODE_APPLY_INVERSE</tt>).
	 * @param  Prec_solve
	 *             [out] A perhaps new <tt>Epetra_Operator</tt> object that
	 *             represents the scaled preconditioner.
	 *
	 * Preconditions:<ul>
	 * <li><tt>scaling.get()!=NULL</tt>
	 * <li><tt>Op.get()!=NULL</tt>
	 * <li><tt>Op_solve!=NULL</tt>
	 * <li><tt>[Prec.get()!=NULL] Prec_solve!=NULL</tt>
	 * </ul>
	 *
	 * Postconditions:<ul>
	 * <li><tt>Op_solve->get()!=NULL</tt>
	 * <li><tt>[Solve_trans==Teuchos::NO_TRANS && Prec.get()!=NULL && fwdSolvePrec()!=FWD_SOLVE_PREC_DEFAULT]</tt>
	 *         <tt>Prec</tt> will be aggregated with <tt>Op</tt> and <tt>Prec_solve->get()==NULL</tt>
	 * <li><tt>[Solve_trans!=Teuchos::NO_TRANS && Prec.get()!=NULL && adjSolvePrec()!=ADJ_SOLVE_PREC_DEFAULT]</tt>
	 *         <tt>Prec</tt> will be aggregated with <tt>Op</tt> and <tt>Prec_solve->get()==NULL</tt>
	 * <li><tt>[Solve_trans==Teuchos::NO_TRANS && Prec.get()!=NULL && fwdSolvePrec()==FWD_SOLVE_PREC_DEFAULT]</tt>
	 *         <tt>Prec</tt> will be scaled independently and will be imbedded in <tt>Prec_solve->get()!=NULL</tt>
	 * <li><tt>[Solve_trans!=Teuchos::NO_TRANS && Prec.get()!=NULL && adjSolvePrec()==FWD_SOLVE_PREC_DEFAULT]</tt>
	 *         <tt>Prec</tt> will be scaled independently and will be imbedded in <tt>Prec_solve->get()!=NULL</tt>
	 * </ul>
	 *
	 * This function returns new objects for the operator and the
	 * preconditioner for a forward or an adjoint solve.  Note, this
	 * function should be called separately to get opeators for the
	 * forward and adjoint solves since <tt>this</tt> may contain
	 * different options for these two systems.  Depending on the
	 * options of <tt>this</tt> set, the preconditioner with either be
	 * aggregated with the operation in the returned
	 * <tt>Op_solve</tt> object or not (see the above
	 * postconditions).
	 */
	void generateSolveOps(
		const Scaling                                 &scaling
		,const Teuchos::ETransp                       scaling_trans
		,const Teuchos::RefCountPtr<Epetra_Operator>  &Op
		,const Teuchos::ETransp                       Op_trans
		,const Teuchos::RefCountPtr<Epetra_Operator>  &Prec
		,const Teuchos::ETransp                       Prec_trans
		,const ProductOperator::EApplyMode            Prec_inverse
		,const Teuchos::ETransp                       Solve_trans
		,Teuchos::RefCountPtr<Epetra_Operator>        *Op_solve
		,const Teuchos::ESide                         Prec_solve_side
		,const ProductOperator::EApplyMode            Prec_solve_inverse
		,Teuchos::RefCountPtr<Epetra_Operator>        *Prec_solve
		,std::ostream                                 *out                         = NULL
		) const;

	///
	/** Transform the RHS before a linear system is solved.
	 *
	 * @param  Op  [in] Contains the scalings and perhaps the
	 *             aggregate preconditioner needed to scale
	 *             the origninal RHS.
	 * @param  Solve_trans
	 *             [in] Determines if the forward (Teuchos::NO_TRANS) or
	 *             adjoint (!Teuchos::NO_TRANS) system will be solved or not.
	 * @param  Rhs [in/out] On input, is the original RHS to the
	 *             linear system to be solved.  On output, the
	 *             scaled and/or preconditioned RHS to be directly
	 *             sent to the linear solver.
	 *
	 * Warning! The client is responsible for the effect of an external
	 * preconditioner that may may be used in the linear solve that is
	 * not embedded in <tt>Op</tt>.  This will be significant if 
	 * a left external preconditioner is used for example.
	 */
	void preSolveTransformRhs(
		Epetra_Operator               &Op_solve
		,const Teuchos::ETransp       Solve_trans
		,Epetra_MultiVector           *Rhs
		) const;

	///
	/** Transform the solution after a linear system is solved.
	 *
	 * @param  Op  [in] Contains the scalings and perhaps the
	 *             aggregate preconditioner needed to scale
	 *             the solved for .
	 * @param  Solve_trans
	 *             [in] Determines if the forward (Teuchos::NO_TRANS) or
	 *             adjoint (!Teuchos::NO_TRANS) system will be solved or not.
	 * @param  Solu
	 *             [in/out] On input, the solved solution from the
	 *             linear solver.  On output, the scaled and/or
	 *             unpreconditioned solution to the original
	 *             linear system.
	 *
	 * Warning! The client is responsible for the effect of an external
	 * preconditioner that may have been used in the linear solve that
	 * is not embedded in <tt>Op</tt>.  This will be significant if a
	 * right external preconditioner is used for example.
	 */
	void postSolveTransformSolu(
		Epetra_Operator               &Op_solve
		,const Teuchos::ETransp       Solve_trans
		,Epetra_MultiVector           *Solu
		) const;

	//@}

};

} // namespace Epetra

#endif // EPETRA_LINEAR_SYSTEM_SCALER_HPP
