  #include "BlockForwardsolver.h"
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

  BlockForwardsolver::BlockForwardsolver(const TSFLinearSolver& s)
	: solvers_(tuple(s))
  {}

  BlockForwardsolver::BlockForwardsolver(const TSFLinearSolver& s0, const TSFLinearSolver& s1)
	: solvers_(tuple(s0, s1))
  {}

  BlockForwardsolver::BlockForwardsolver(const TSFLinearSolver& s0, const TSFLinearSolver& s1, const TSFLinearSolver& s2)
	: solvers_(tuple(s0, s1, s2))
  {}

  BlockForwardsolver::BlockForwardsolver(const TSFArray<TSFLinearSolver>& s)
	: solvers_(s)
  {}


  bool BlockForwardsolver::solve(const TSFLinearOperator& op, const TSFVector& rhs, TSFVector& soln) const
  {
   int n = op.numBlockRows();
   if (n != op.numBlockCols())
   {			
   TSFError::raise("BlockBacksolver::solve given non-square block operator");
   }

   soln = rhs.space().createMember();
 
   for (int i=0; i<n; i++)
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

   for (int j=i-1; j>=0; j--)
   {
     if (op.getBlock(i,j).isZeroOperator()) continue;
     tmp = tmp - op.getBlock(i,j)*soln.getBlock(j);
   }

   TSFVector x = soln.getBlock(i);

   if (!solver.solve(op.getBlock(i,i), tmp, x))
      {
      TSFOut::println("Block backsolve failed at diagonal block " + TSF::toString(i));
      return false;
      }

      soln.setBlock(i, x);

   }
   return true;
  }
