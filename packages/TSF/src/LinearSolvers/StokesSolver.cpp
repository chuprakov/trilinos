#include "StokesSolver.h"
#include "TSFUtils.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFBlockLinearOperator.h"
#include "TSFIdentityOperator.h"
#include "TSFZeroOperator.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

#include "TSFBlas.h"
#include "TSFError.h"
#include "TSFOut.h"
#include "TSFMPI.h"

using namespace TSF;


StokesSolver::StokesSolver(const TSFLinearSolver& innerSolver,
													 const TSFLinearSolver& outerSolver)
	: TSFLinearSolverBase(), innerSolver_(innerSolver), outerSolver_(outerSolver)
{}


bool StokesSolver::solve(const TSFLinearOperator& op,
												 const TSFVector& rhs,
												 TSFVector& soln) const 
{
	try
		{
			if (op.numBlockRows() != 2 || op.numBlockCols()!=2)
				{
					TSFError::raise("StokesSolver::solve expects a 2x2 block operator");
				}
			if (rhs.numBlocks() != 2)
				{
					TSFError::raise("StokesSolver::solve expects a 2-block vector"); 
				}
			
			TSFLinearOperator A = op.getBlock(0,0);
			TSFLinearOperator B = op.getBlock(0,1);			
			TSFLinearOperator C = op.getBlock(1,0);
			TSFLinearOperator D = op.getBlock(1,1);
			
			TSFVector a = rhs.getBlock(0);
			TSFVector b = rhs.getBlock(1);

			TSFLinearOperator Ainv = A.inverse(innerSolver_);

			TSFOut::println("forming modified operator");
			TSFLinearOperator modifiedOp 
				= new TSFBlockLinearOperator(op.domain(), op.range());
			
			modifiedOp.setBlock(0, 0, new TSFIdentityOperator(A.range()));
			modifiedOp.setBlock(0, 1, Ainv*B);
			modifiedOp.setBlock(1, 0, C);
			modifiedOp.setBlock(1, 1, new TSFZeroOperator(D.domain(), D.range()));
			
			TSFOut::println("forming modified RHS");
			TSFVector modifiedRHS = op.domain().createMember();
			modifiedRHS.setBlock(0, Ainv*a);
			modifiedRHS.setBlock(1, b);

			TSFOut::println("starting outer solve");
			bool rtn = outerSolver_.solve(modifiedOp, modifiedRHS, soln);
			if (!rtn)
				{
					TSFOut::println("outer solve failed in StokesSolver::solve()");
				}
			return rtn;
		}
	catch(exception& e)
		{
			TSFOut::println(e.what());
			return false;
		}
}

