#include "TSFInverseOperator.h"

#include "TSFError.h"
#include "TSFUtils.h"
#include "TSFInverseAdjointOperator.h"
#include "TSFAdjointOperator.h"
#include "TSFLinearSolver.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"

using namespace TSF;


TSFInverseOperator::TSFInverseOperator(const TSFLinearOperator& op,
																			 const TSFLinearSolver& solver)
	: TSFLinearOperatorBase(op.range(), op.domain()), op_(op), solver_(solver)
{}

void TSFInverseOperator::apply(const TSFVector& in, 
															 TSFVector& out) const
{
	try
		{
			if (solver_.isNull())
				{
					out = domain().createMember();
					op_.applyInverse(in, out);
				}
			else
				{
					op_.applyInverse(solver_, in, out);
				}
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseOperator::apply()");
		}
}

void TSFInverseOperator::applyAdjoint(const TSFVector& in, 
																			TSFVector& out) const
{
	try
		{
			if (solver_.isNull())
				{
					op_.applyInverseAdjoint(in, out);
				}
			else
				{
					op_.applyInverseAdjoint(solver_, in, out);
				}
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseOperator::applyAdjoint()");
		}
}

void TSFInverseOperator::applyInverse(const TSFVector& in, 
																			TSFVector& out) const
{
	try
		{
			op_.apply(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseOperator::applyInverse()");
		}
}



void TSFInverseOperator::applyInverseAdjoint(const TSFVector& in, 
																						 TSFVector& out) const
{
	try
		{
			op_.applyAdjoint(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseOperator::applyInverseAdjoint()");
		}
}




void TSFInverseOperator::getInverse(const TSFLinearOperator& /* self */,
																		TSFLinearOperator& inv) const 
{
	inv = op_;
}

void TSFInverseOperator::getInverse(const TSFLinearSolver& /* solver */,
																		const TSFLinearOperator& /* self */,
																		TSFLinearOperator& inv) const 
{
	inv = op_;
}



void TSFInverseOperator::getAdjoint(const TSFLinearOperator& /* self */,
																		TSFLinearOperator& adj) const 
{
	adj = new TSFInverseAdjointOperator(op_, solver_);
}

void TSFInverseOperator::getInverseAdjoint(const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& invAdj) const 
{
	invAdj = new TSFAdjointOperator(op_);
}

void TSFInverseOperator::getInverseAdjoint(const TSFLinearSolver& /* solver */,
																					 const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& invAdj) const 
{
	invAdj = new TSFAdjointOperator(op_);
}






