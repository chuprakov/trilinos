#include "TSFAdjointOperator.h"

#include "TSFError.h"
#include "TSFUtils.h"
#include "TSFInverseOperator.h"
#include "TSFInverseAdjointOperator.h"
#include "TSFLinearSolver.h"
#include "TSFVectorSpaceBase.h"
#include "TSFAdjointMatrixOperator.h"
using namespace TSF;


TSFAdjointOperator::TSFAdjointOperator(const TSFLinearOperator& op)
	: TSFLinearOperatorBase(op.range(), op.domain()), op_(op)
{}

void TSFAdjointOperator::apply(const TSFVector& in, 
															 TSFVector& out) const
{
	try
		{
			op_.applyAdjoint(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFAdjointOperator::apply()");
		}
}

void TSFAdjointOperator::applyAdjoint(const TSFVector& in, 
																			TSFVector& out) const
{
	try
		{
			op_.apply(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFAdjointOperator::applyAdjoint()");
		}
}

void TSFAdjointOperator::applyInverse(const TSFVector& in, 
																			TSFVector& out) const
{
	try
		{
			op_.applyInverseAdjoint(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFAdjointOperator::applyInverse()");
		}
}

void TSFAdjointOperator::applyInverse(const TSFLinearSolver& solver,
																			const TSFVector& in, 
																			TSFVector& out) const
{
	try
		{
			op_.applyInverseAdjoint(solver, in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFAdjointOperator::applyInverse()");
		}
}

void TSFAdjointOperator::applyInverseAdjoint(const TSFVector& in, 
																						 TSFVector& out) const
{
	try
		{
			op_.applyInverse(in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFAdjointOperator::applyInverseAdjoint()");
		}
}


void TSFAdjointOperator::applyInverseAdjoint(const TSFLinearSolver& solver,
																						 const TSFVector& in, 
																						 TSFVector& out) const
{
	try
		{
			op_.applyInverse(solver, in, out);
		}
	catch(exception& e)
		{
			TSFError::trace(e, "in TSFAdjointOperator::applyInverseAdjoint()");
		}
}

void TSFAdjointOperator::getInverse(const TSFLinearOperator& /* self */,
																		TSFLinearOperator& inv) const 
{
	inv = new TSFInverseAdjointOperator(op_);
}

void TSFAdjointOperator::getInverse(const TSFLinearSolver& solver,
																		const TSFLinearOperator& /* self */,
																		TSFLinearOperator& inv) const 
{
	inv = new TSFInverseAdjointOperator(op_, solver);
}


void TSFAdjointOperator::getAdjoint(const TSFLinearOperator& /* self */,
																		TSFLinearOperator& adj) const 
{
	adj = op_;
}

void TSFAdjointOperator::getInverseAdjoint(const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& invAdj) const 
{
	invAdj = new TSFInverseOperator(op_);
}

void TSFAdjointOperator::getInverseAdjoint(const TSFLinearSolver& solver,
																					 const TSFLinearOperator& /* self */,
																					 TSFLinearOperator& invAdj) const 
{
	invAdj = new TSFInverseOperator(op_, solver);
}

bool TSFAdjointOperator::isMatrixOperator() const
{
	return op_.isMatrixOperator();	
}

const TSFSmartPtr<const TSFMatrixOperator>
TSFAdjointOperator::getMatrix() const
{
	TSFSmartPtr<const TSFAdjointMatrixOperator> adj_mat_op;
	if(op_.isMatrixOperator()) {
		TSFSmartPtr<TSFMatrixOperator> mat_op;
		mat_op = tsf_const_cast<TSFMatrixOperator>(op_.getMatrix());
		adj_mat_op = new TSFAdjointMatrixOperator(mat_op);
	}
	return tsf_implicit_cast<const TSFMatrixOperator>(adj_mat_op);
}



void TSFAdjointOperator::getRow(int row, TSFArray<int>& indices, 
                                TSFArray<TSFReal>& values) const
{
  /* very inefficient method of getting row  */
  try
    {
      TSFVector e = op_.domain().createMember();
      e.zero();
      e[row] = 1.0;
      TSFVector r = op_.range().createMember();
      r = op_ * e;

      for (int i = 0; i < r.space().dim(); i++)
        {
          if (r[i] != 0)
            {
              indices.append(i);
              values.append(r[i]);
            }
        }
    }
  catch(exception& e)
    {
      TSFError::trace(e, "in TSFAdjointMatrixOperator::getRow()");
    }

}
