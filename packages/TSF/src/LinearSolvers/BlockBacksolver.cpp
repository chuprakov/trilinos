#include "BlockBacksolver.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFError.h"
#include "TSFOut.h"
#include "TSFDeferredLinearCombination.h"

using namespace TSF;

BlockBacksolver::BlockBacksolver(const TSFLinearSolver& s)
	: solvers_(tuple(s))
{}

BlockBacksolver::BlockBacksolver(const TSFLinearSolver& s0, 
																 const TSFLinearSolver& s1)
	: solvers_(tuple(s0, s1))
{}

BlockBacksolver::BlockBacksolver(const TSFLinearSolver& s0, 
																 const TSFLinearSolver& s1,
																 const TSFLinearSolver& s2)
	: solvers_(tuple(s0, s1, s2))
{}

BlockBacksolver::BlockBacksolver(const TSFArray<TSFLinearSolver>& s)
	: solvers_(s)
{}



bool BlockBacksolver::solve(const TSFLinearOperator& op, 
														const TSFVector& rhs, 
														TSFVector& soln) const
{
	int n = op.numBlockRows();
	if (n != op.numBlockCols())
		{
			TSFError::raise("BlockBacksolver::solve given non-square block operator");
		}

	soln = rhs.space().createMember();

	for (int i=n-1; i>=0; i--)
		{
			cerr << "solving block " << i << endl;
			TSFLinearSolver solver;
			if (solvers_.length()==1)
				{
					solver = solvers_[0];
				}
			else
				{
					solver = solvers_[i];
				}
			
			TSFVector bi = rhs.getBlock(i);
			TSFVector tmp = bi;
			for (int j=i+1; j<n; j++)
				{
					if (op.getBlock(i,j).isZeroOperator()) continue;
					tmp = tmp - op.getBlock(i,j)*soln.getBlock(j);
				}
			
			TSFVector x;
			if (!solver.solve(op.getBlock(i,i), tmp, x))
				{
					TSFOut::println("Block backsolve failed at diagonal block "
												+ TSF::toString(i));
					return false;
				}
					
			soln.setBlock(i, x);
		}
	return true;
}
