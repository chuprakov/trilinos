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
// Epetra_LinearSystemScaler.cpp

#include "Epetra_LinearSystemScaler.hpp"
#include "Epetra_DiagonalOperator.hpp"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"
#include "Teuchos_dyn_cast.hpp"

namespace {

enum EOpType { PLAIN_OP, LEFT_PREC_OP, RIGHT_PREC_OP, LEFT_RIGHT_PREC_OP };

inline
Teuchos::ETransp trans_trans( Teuchos::ETransp trans1, Teuchos::ETransp trans2 )
{
	if( trans1 == trans2 ) return Teuchos::NO_TRANS;
	return Teuchos::TRANS;
}

typedef Epetra::LinearSystemScaler ELSS;

inline
const char* toString( ELSS::EFwdSolveLeftScaling o )
{ return ( o==ELSS::FWD_SOLVE_LEFT_SCALING_NONE ? "FWD_SOLVE_LEFT_SCALING_NONE" : "FWD_SOLVE_LEFT_SCALING_ROW_SUM" ); }

inline
const char* toString( ELSS::EFwdSolveRightScaling o )
{ return ( o==ELSS::FWD_SOLVE_RIGHT_SCALING_NONE ? "FWD_SOLVE_RIGHT_SCALING_NONE" : "FWD_SOLVE_RIGHT_SCALING_COL_SUM" ); }

inline
const char* toString( ELSS::EAdjSolveLeftScaling o )
{ return ( o==ELSS::ADJ_SOLVE_LEFT_SCALING_NONE ? "ADJ_SOLVE_LEFT_SCALING_NONE" : "ADJ_SOLVE_LEFT_SCALING_COL_SUM" ); }

inline
const char* toString( ELSS::EAdjSolveRightScaling o )
{ return ( o==ELSS::ADJ_SOLVE_RIGHT_SCALING_NONE ? "ADJ_SOLVE_RIGHT_SCALING_NONE" : "ADJ_SOLVE_RIGHT_SCALING_ROW_SUM" ); }

inline
const char* toString( ELSS::EFwdSolvePrec o )
{ return ( o==ELSS::FWD_SOLVE_PREC_DEFAULT ? "FWD_SOLVE_PREC_DEFAULT" : ( o==ELSS::FWD_SOLVE_PREC_LEFT ? "FWD_SOLVE_PREC_LEFT" : "FWD_SOLVE_PREC_RIGHT" ) ); }

inline
const char* toString( ELSS::EAdjSolvePrec o )
{ return ( o==ELSS::ADJ_SOLVE_PREC_DEFAULT ? "ADJ_SOLVE_PREC_DEFAULT" : ( o==ELSS::ADJ_SOLVE_PREC_LEFT ? "ADJ_SOLVE_PREC_LEFT" : "ADJ_SOLVE_PREC_RIGHT" ) ); }

const std::string fwd_op_type_name = "fwd_op_type";
const std::string adj_op_type_name = "adj_op_type";

} // namespace

namespace Epetra {

LinearSystemScaler::LinearSystemScaler (
	const EFwdSolveLeftScaling     fwdSolveLeftScaling
	,const EFwdSolveRightScaling   fwdSolveRightScaling
	,const EAdjSolveLeftScaling    adjSolveLeftScaling
	,const EAdjSolveRightScaling   adjSolveRightScaling
	,const EFwdSolvePrec           fwdSolvePrec
	,const EAdjSolvePrec           adjSolvePrec
	)
	:fwdSolveLeftScaling_(fwdSolveLeftScaling)
	,fwdSolveRightScaling_(fwdSolveRightScaling)
	,adjSolveLeftScaling_(adjSolveLeftScaling)
	,adjSolveRightScaling_(adjSolveRightScaling)
	,fwdSolvePrec_(fwdSolvePrec)
	,adjSolvePrec_(adjSolvePrec)
{}

// Scaling methods

Teuchos::RefCountPtr<const LinearSystemScaler::Scaling>
LinearSystemScaler::computeScaling(
	const Epetra_RowMatrix   &Op
	,const Teuchos::ETransp  Op_trans
	) const
{
	Teuchos::RefCountPtr<Epetra_Vector>
		leftScaling  = Teuchos::rcp(new Epetra_Vector(Op.OperatorRangeMap())),
		rightScaling = Teuchos::rcp(new Epetra_Vector(Op.OperatorDomainMap()));
	TEST_FOR_EXCEPTION(0!=Op.InvRowSums(*leftScaling),std::runtime_error,"Error!");
	TEST_FOR_EXCEPTION(0!=Op.InvColSums(*rightScaling),std::runtime_error,"Error!");
	return Teuchos::rcp(
		new Scaling(
			Op_trans == Teuchos::NO_TRANS  ? leftScaling  : rightScaling
			,Op_trans == Teuchos::NO_TRANS ? rightScaling : leftScaling
			)
		);
}

void LinearSystemScaler::generateSolveOps(
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
	,std::ostream                                 *out
	) const
{
#ifdef _DEBUG
	TEST_FOR_EXCEPTION(
		Op.get()==NULL || Op_solve==NULL || ( Prec.get()!=NULL && Prec_solve==NULL )
		,std::invalid_argument
		,"LinearSystemScaler::applySolveScalings(...): Error!"
		);
	// ToDo: Validate maps of operator and preconditioner!
#endif // _DEBUG
	//
	// Sort through the options and decide what you are going
	// to do.
	//
	if(out) {
		*out << "\n*** Entering Epetra::LinearSystemScaler::generateSolveOps(...) ...\n"
				 << "\nGenerating scaled and/or aggregate operator and preconditioner for "
				 << (Solve_trans==Teuchos::NO_TRANS?"forward":"adjont") << " solve with options:";
		if(Solve_trans==Teuchos::NO_TRANS) {
			*out << "\n  fwdSolveLeftScaling  = " << toString(fwdSolveLeftScaling())
					 << "\n  fwdSolveRightScaling = " << toString(fwdSolveRightScaling())
					 << "\n  fwdSolvePrec         = " << toString(fwdSolvePrec());
		}
		else {
			*out << "\n  adjSolveLeftScaling  = " << toString(adjSolveLeftScaling())
					 << "\n  adjSolveRightScaling = " << toString(adjSolveRightScaling())
					 << "\n  adjSolvePrec         = " << toString(adjSolvePrec());
		}
		*out << std::endl;
	}
	
	// ToDo: Figure this out!

	//
	// Setup the aggregate operator
	//
	// Here we will setup the following product opeator:
	//
	//    M = L * A * R
	//
	//        where:
	//                L is the left preconditioner and/or scaling
	//                A is the original operator Op
	//                R is the right precondtiioner and/or scaling
	//                M is the aggregate output operator Op_solve
	//
	Teuchos::RefCountPtr<Epetra_Operator>   Op_array[3];
	Teuchos::ETransp                        Op_trans_array[3];
	ProductOperator::EApplyMode             Op_inverse_array[3];
	int ops_added = 0;
	//
	// Add the aggregate left preconditioner and/or scaling operator
	//
	bool hasLeftPrec = false;
	Teuchos::RefCountPtr<Epetra_Operator> leftScalingOp;
	if(
		( Solve_trans==Teuchos::NO_TRANS
			&& fwdSolveLeftScaling()!=FWD_SOLVE_LEFT_SCALING_NONE
			&& fwdSolvePrec()==FWD_SOLVE_PREC_DEFAULT
			)
		||
		( Solve_trans!=Teuchos::NO_TRANS
			&& adjSolveLeftScaling()!=ADJ_SOLVE_LEFT_SCALING_NONE
			&& adjSolvePrec()==ADJ_SOLVE_PREC_DEFAULT
			)
		)
	{
		if(out)
			*out << "\nGenerating a simple scaling operator on the left ...\n";
		leftScalingOp
			= Teuchos::rcp(
				new DiagonalOperator(
					Teuchos::rcp(	Solve_trans==Teuchos::NO_TRANS ? &Op->OperatorRangeMap() : &Op->OperatorDomainMap(), false )
					,(Solve_trans==Teuchos::NO_TRANS)==(scaling_trans) ? scaling.leftScaling() : scaling.rightScaling()
					)
				);
		Op_array[ops_added] = leftScalingOp;
		Op_trans_array[ops_added] = Teuchos::NO_TRANS;
		Op_inverse_array[ops_added] =  ProductOperator::APPLY_MODE_APPLY;
		++ops_added;
		hasLeftPrec = true;
	}
	else if(
		( Solve_trans==Teuchos::NO_TRANS
			&& fwdSolveLeftScaling()==FWD_SOLVE_LEFT_SCALING_NONE
			&& fwdSolvePrec()==FWD_SOLVE_PREC_LEFT
			)
		||
		( Solve_trans!=Teuchos::NO_TRANS
			&& adjSolveLeftScaling()==ADJ_SOLVE_LEFT_SCALING_NONE
			&& adjSolvePrec()==ADJ_SOLVE_PREC_LEFT
			)
		)
	{
		if(out)
			*out << "\nAggregating the preconditioner as a left preconditioner with the operator ...\n";
		Op_array[ops_added] = Prec;
		Op_trans_array[ops_added] = trans_trans(Prec_trans,Solve_trans);
		Op_inverse_array[ops_added] = Prec_inverse;
		++ops_added;
		hasLeftPrec = true;
	}
	else {
		// ToDo: Implement more options!
	}
	//
	// Add the operator
	//
	if(1) {
		Op_array[ops_added] = Op;
		Op_trans_array[ops_added] = trans_trans(Op_trans,Solve_trans);
		Op_inverse_array[ops_added] = ProductOperator::APPLY_MODE_APPLY;
		++ops_added ;
	}
	//
	// Add the aggregate right preconditioner and/or scaling
	//
	bool hasRightPrec = false;
	Teuchos::RefCountPtr<Epetra_Operator> rightScalingOp;
	if(
		( Solve_trans==Teuchos::NO_TRANS
			&& fwdSolveRightScaling()!=FWD_SOLVE_RIGHT_SCALING_NONE
			&& fwdSolvePrec()==FWD_SOLVE_PREC_DEFAULT
			)
		||
		( Solve_trans!=Teuchos::NO_TRANS
			&& adjSolveRightScaling()!=ADJ_SOLVE_RIGHT_SCALING_NONE
			&& adjSolvePrec()==ADJ_SOLVE_PREC_DEFAULT
			)
		)
	{
		if(out)
			*out << "\nGenerating a simple scaling operator on the right ...\n";
		rightScalingOp 
			= Teuchos::rcp(
				new DiagonalOperator(
					Teuchos::rcp(	Solve_trans==Teuchos::NO_TRANS ? &Op->OperatorDomainMap() : &Op->OperatorRangeMap(), false )
					,(Solve_trans==Teuchos::NO_TRANS)==(scaling_trans) ? scaling.rightScaling() : scaling.leftScaling()
					)
				);
		Op_array[ops_added] = rightScalingOp;
		Op_trans_array[ops_added] = Teuchos::NO_TRANS;
		Op_inverse_array[ops_added] = ProductOperator::APPLY_MODE_APPLY;
		++ops_added;
		hasRightPrec = true;
	}
	else if(
		( Solve_trans==Teuchos::NO_TRANS
			&& fwdSolveRightScaling()==FWD_SOLVE_RIGHT_SCALING_NONE
			&& fwdSolvePrec()==FWD_SOLVE_PREC_RIGHT
			)
		||
		( Solve_trans!=Teuchos::NO_TRANS
			&& adjSolveRightScaling()==ADJ_SOLVE_RIGHT_SCALING_NONE
			&& adjSolvePrec()==ADJ_SOLVE_PREC_RIGHT
			)
		)
	{
		if(out)
			*out << "\nAggregating the preconditioner as a right preconditioner with the operator ...\n";
		Op_array[ops_added] = Prec;
		Op_trans_array[ops_added] = trans_trans(Prec_trans,Solve_trans);
		Op_inverse_array[ops_added] = Prec_inverse;
		++ops_added;
		hasRightPrec = true;
	}
	else {
		// ToDo: Implement more options!
	}
	// Remember what type of preconditioned operator this is
	// and store it in the smart point of the first 
	// RefCountPtr array argument.  That way, we can
	// retrieve it later.
	EOpType opType;
	if( hasLeftPrec ) {
		if( hasRightPrec ) {
			opType = LEFT_RIGHT_PREC_OP;
		}
		else {
			opType = LEFT_PREC_OP;
		}
	}
	else if( hasRightPrec ) {
		opType = RIGHT_PREC_OP;
	}
	else {
		opType = PLAIN_OP;
	}
	Teuchos::set_extra_data<EOpType>(
		opType,Solve_trans==Teuchos::NO_TRANS ? fwd_op_type_name : adj_op_type_name, &Op_array[0], Teuchos::POST_DESTROY, false
		);
	// Create the aggregate preconditioner and/or scaled operator
	*Op_solve = Teuchos::rcp(
		new Epetra::ProductOperator( ops_added, Op_array, Op_trans_array, Op_inverse_array )
		);
	//
	// Setup the external preconditioner
	//
	if(
		Prec.get()
		&&
		(
			(Solve_trans==Teuchos::NO_TRANS && fwdSolvePrec()==FWD_SOLVE_PREC_DEFAULT)
			|| 
			(Solve_trans!=Teuchos::NO_TRANS && adjSolvePrec()==ADJ_SOLVE_PREC_DEFAULT)
			)
		)
	{
		if(out)
			*out << "\nGenerating an external preconditioner ...\n";
		Teuchos::RefCountPtr<Epetra_Operator> Op_array[3];
		Teuchos::ETransp Op_trans_array[3];
		ProductOperator::EApplyMode Op_inverse_array[3];
		int ops_added = 0;
		if( Prec_solve_inverse == ProductOperator::APPLY_MODE_APPLY_INVERSE ) {
			//
			// Here the aggregate preconditioner will be applied as its its
			// inverse.  Therefore, we must scale the precondtioner just
			// like we scaled the operator.
			//
			// Add left scaling
			if(
				( Solve_trans==Teuchos::NO_TRANS
					&& fwdSolveLeftScaling()!=FWD_SOLVE_LEFT_SCALING_NONE
					)
				||
				( Solve_trans!=Teuchos::NO_TRANS
					&& adjSolveLeftScaling()!=ADJ_SOLVE_LEFT_SCALING_NONE
					)
				)
			{
				// We scaled the operator by the left so we need to scale the
				// external preconditioner by the inverse scaling on the right
				Op_array[ops_added] = leftScalingOp;
				Op_trans_array[ops_added] = Teuchos::NO_TRANS;
				Op_inverse_array[ops_added] = ProductOperator::APPLY_MODE_APPLY;
				++ops_added;
			}
			// Add the preconditioner itself
			if(1) {
				Op_array[ops_added] = Prec;
				Op_trans_array[ops_added] = trans_trans(Prec_trans,Solve_trans);
				Op_inverse_array[ops_added]
					= (
						Prec_inverse==ProductOperator::APPLY_MODE_APPLY
						?ProductOperator::APPLY_MODE_APPLY_INVERSE
						:ProductOperator::APPLY_MODE_APPLY
						);
				++ops_added;
			}
			*Prec_solve = Teuchos::rcp(
				new Epetra::ProductOperator( ops_added, Op_array, Op_trans_array, Op_inverse_array )
				);
			// Add right scaling
			// ToDo: Implement!
		}
		else {
			assert(0); // ToDo: Implement this case if every needed!
		}
	}
	if(out)
		*out << "\n*** Leaving Epetra::LinearSystemScaler::generateSolveOps(...)\n";
}

void LinearSystemScaler::preSolveTransformRhs(
	Epetra_Operator         &Op_solve
	,const Teuchos::ETransp Solve_trans
	,Epetra_MultiVector     *Rhs
	) const
{
	using Teuchos::dyn_cast;
	Epetra::ProductOperator &prod_Op = dyn_cast<Epetra::ProductOperator>(Op_solve);
	if( prod_Op.num_Op() == 1 ) {
		// There is only the original operator so no transformation is needed!
	}
	else {
		const EOpType opType
			= Teuchos::get_extra_data<EOpType>(
				prod_Op.Op(0), Solve_trans==Teuchos::NO_TRANS ? fwd_op_type_name : adj_op_type_name
				);
		if( opType==LEFT_PREC_OP || opType==LEFT_RIGHT_PREC_OP ) {
			Epetra_MultiVector  scaledRhs(Rhs->Map(),Rhs->NumVectors());  // ToDo: Make cache data member
			prod_Op.applyConstituent(0,Teuchos::NO_TRANS,ProductOperator::APPLY_MODE_APPLY,*Rhs,&scaledRhs);
			*Rhs = scaledRhs;
		}
	}
}

void LinearSystemScaler::postSolveTransformSolu(
	Epetra_Operator         &Op_solve
	,const Teuchos::ETransp Solve_trans
	,Epetra_MultiVector     *Solu
	) const
{
	using Teuchos::dyn_cast;
	Epetra::ProductOperator
		&prod_Op = dyn_cast<Epetra::ProductOperator>(Op_solve);
	if( prod_Op.num_Op() == 1 ) {
		// There is only the original operator so no transformation is needed!
	}
	else {
		const EOpType opType
			= Teuchos::get_extra_data<EOpType>(
				prod_Op.Op(0), Solve_trans==Teuchos::NO_TRANS ? fwd_op_type_name : adj_op_type_name
				);
		if( opType==RIGHT_PREC_OP || opType==LEFT_RIGHT_PREC_OP ) {
			Epetra_MultiVector  scaledSolu(Solu->Map(),Solu->NumVectors());  // ToDo: Make cache data member
			prod_Op.applyConstituent(
				opType==RIGHT_PREC_OP?1:2
				,Teuchos::NO_TRANS,ProductOperator::APPLY_MODE_APPLY,*Solu,&scaledSolu
				);
			*Solu = scaledSolu;
		}
	}
}

} // namespace Epetra
