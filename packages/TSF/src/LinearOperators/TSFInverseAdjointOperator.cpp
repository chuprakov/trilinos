#include "TSFInverseAdjointOperator.h"

#include "TSFError.h"
#include "TSFUtils.h"
#include "TSFInverseOperator.h"
#include "TSFAdjointOperator.h"
#include "TSFLinearSolver.h"
#include "TSFVectorSpaceBase.h"

using namespace TSF;


TSFInverseAdjointOperator::TSFInverseAdjointOperator(const TSFLinearOperator& op,
																			 const TSFLinearSolver& solver)
	: TSFLinearOperatorBase(op.domain(), op.range()), op_(op), solver_(solver)
{}

void TSFInverseAdjointOperator::apply(const TSFVector& in, 
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
			TSFError::trace(e, "in TSFInverseAdjointOperator::apply()");
		}
}

void TSFInverseAdjointOperator::applyAdjoint(const TSFVector& in, 
																						 TSFVector& out) const
{
	try
		{
			if (solver_.isNull())
				{
					op_.applyInverse(in, out);
				}
			else
				{
					op_.applyInverse(solver_, in, out);
				}
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseAdjointOperator::applyAdjoint()");
		}
}

void TSFInverseAdjointOperator::applyInverse(const TSFVector& in, 
																						 TSFVector& out) const
{
	try
		{
			op_.applyAdjoint(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseAdjointOperator::applyInverse()");
		}
}

void TSFInverseAdjointOperator::applyInverse(const TSFLinearSolver& /* solver */,
																			const TSFVector& in, 
																			TSFVector& out) const
{
	try
		{
			op_.apply(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseAdjointOperator::applyInverse()");
		}
}

void TSFInverseAdjointOperator::applyInverseAdjoint(const TSFVector& in, 
																										TSFVector& out) const
{
	try
		{
			op_.apply(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseAdjointOperator::applyInverseAdjoint()");
		}
}


void TSFInverseAdjointOperator::applyInverseAdjoint(const TSFLinearSolver& /*solver*/,
																						 const TSFVector& in, 
																						 TSFVector& out) const
{
	try
		{
			op_.apply(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFInverseAdjointOperator::applyInverseAdjoint()");
		}
}

void TSFInverseAdjointOperator::getInverse(const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& inv) const 
{
	inv = new TSFAdjointOperator(op_);
}

void TSFInverseAdjointOperator::getInverse(const TSFLinearSolver& /* solver */,
																					 const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& inv) const 
{
	inv = new TSFAdjointOperator(op_);
}



void TSFInverseAdjointOperator::getAdjoint(const TSFLinearOperator& /* self */,
																		TSFLinearOperator& adj) const 
{
	adj = new TSFInverseOperator(op_, solver_);
}

void TSFInverseAdjointOperator::getInverseAdjoint(const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& invAdj) const 
{
	invAdj = op_;
}

void TSFInverseAdjointOperator::getInverseAdjoint(const TSFLinearSolver& /* solver */,
																					 const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& invAdj) const 
{
	invAdj = op_;
}






