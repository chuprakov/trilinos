#include "SchurComplementSolver.h"
#include "TSFUtils.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

#include "TSFBlas.h"
#include "TSFError.h"
#include "TSFOut.h"
#include "TSFMPI.h"
#include "TSFDeferredLinearCombination.h"

using namespace TSF;


SchurComplementSolver::SchurComplementSolver(const TSFLinearSolver& innerSolver,
																						 const TSFLinearSolver& outerSolver)
	: TSFLinearSolverBase(), innerSolver_(innerSolver), outerSolver_(outerSolver)
{}

SchurComplementSolver::~SchurComplementSolver(){;}

bool SchurComplementSolver::solve(const TSFLinearOperator& op,
																	const TSFVector& rhs,
																	TSFVector& soln) const 
{
	try
		{
			if (op.numBlockRows() != 2 || op.numBlockCols()!=2)
				{
					TSFError::raise("SchurComplementSolver::solve expects a 2x2 block operator");
				}
			if (rhs.numBlocks() != 2)
				{
					TSFError::raise("SchurComplementSolver::solve expects a 2-block vector"); 
				}
			
			TSFLinearOperator A = op.getBlock(0,0);
			TSFLinearOperator B = op.getBlock(0,1);			
			TSFLinearOperator C = op.getBlock(1,0);
			TSFLinearOperator D = op.getBlock(1,1);
			
			TSFVector a = rhs.getBlock(0);
			TSFVector b = rhs.getBlock(1);
			
			TSFLinearOperator Ainv = A.inverse(innerSolver_);
			
			TSFVector p = b - C*Ainv*a;
			
			TSFVector y = (D - C*Ainv*B).inverse(outerSolver_)*p;
			
			TSFVector x = Ainv*(a - B*y);
			
			soln = rhs.space().createMember();
			soln.setBlock(0, x);
			soln.setBlock(1, y);
			return true;
		}
	catch(exception& e)
		{
			TSFOut::println(e.what());
			return false;
		}
}

