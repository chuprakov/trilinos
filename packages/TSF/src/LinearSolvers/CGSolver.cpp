#include "CGSolver.h"
#include "TSFUtils.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"

#include "TSFParameterListImplem.h"
#include "TSFBlas.h"
#include "TSFError.h"
#include "TSFOut.h"
#include "TSFMPI.h"
#include "TSFDeferredLinearCombination.h"

using namespace TSF;


CGSolver::CGSolver(const TSFParameterList& params)
	: TSFLinearSolverBase(params)
{;}

TSFParameterList CGSolver::defaultParameters() const 
{
  TSFParameterList rtn("CG Parameters");
  rtn.addParameter(TSFParameter("max iterations", 500));
  rtn.addParameter(TSFParameter("convergence tolerance", 1.0e-10));
  return rtn;
}

bool CGSolver::solve(const TSFLinearOperator& A,
										 const TSFVector& b,
										 TSFVector& soln) const 
{
  
  int maxiters = params_.getParameter("max iterations").getInt();
  double tol = params_.getParameter("convergence tolerance").getInt();

	TSFReal normOfB = sqrt(b.dot(b));

	/* check for trivial case of zero rhs */
	if (normOfB < tol) 
		{
			soln = b.space().createMember();
			soln.zero();
			return true;
		}

	/* check for initial zero residual */
	TSFVector x0 = b.space().createMember();
	TSFVector r = b.copy();
	TSFVector p = r.copy();


	int myRank = TSFMPI::getRank();

	for (int k=1; k<=maxiters; k++)
		{
          //if (myRank==0) TSFOut::printf("iteration %d\n", k);
			TSFVector Ap = A*p;
			TSFReal r20 = r*r;
			TSFReal pAp = p*Ap;
			if (TSFUtils::chop(pAp)==0) 
				{
                  cerr << "CG Iteration " << k << endl;
                  //TSFError::raise("failure mode 1 in CG");
                  pAp = 1.0e-16;
				}
			TSFReal alpha = r20/pAp;
			x0 = x0 + alpha*p;
			r = r - alpha*Ap;
			TSFReal r21 = r*r;

			if (sqrt(r21) < tol*normOfB) 
				{
					soln = x0.copy();
					return true;
				}
			TSFReal beta = r21/r20;
			p = r + beta*p;
		}
	TSFError::raise("CG failed to converge");
	return false;
}
	



