#include "BICGSTABSolver.h"
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
#include "TSFPreconditionerFactoryBase.h"
#include "TSFPreconditionerFactory.h"
#include "TSFDeferredLinearCombination.h"

using namespace TSF;


BICGSTABSolver::BICGSTABSolver(const TSFParameterList& params)
  : TSFLinearSolverBase(params), 
    preconditionerFactory_()
{;}

BICGSTABSolver::BICGSTABSolver(const TSFPreconditionerFactory& pf,
			       const TSFParameterList& params)
  : TSFLinearSolverBase(params),
    preconditionerFactory_(pf)
{;}

BICGSTABSolver::BICGSTABSolver(const double& tol, int maxiters)
  : TSFLinearSolverBase(), 
    preconditionerFactory_()
{
  params_ = defaultParameters();
  params_.setValue("max iterations", maxiters);
  params_.setValue("convergence tolerance", tol);
}

BICGSTABSolver::BICGSTABSolver(const TSFPreconditionerFactory& pf,
                               const double& tol, int maxiters)
  : TSFLinearSolverBase(), 
    preconditionerFactory_(pf)
{
  params_ = defaultParameters();
  params_.setValue("max iterations", maxiters);
  params_.setValue("convergence tolerance", tol);
}

BICGSTABSolver::~BICGSTABSolver(){;}

TSFParameterList BICGSTABSolver::defaultParameters() const 
{
  TSFParameterList rtn("BICGSTAB Parameters");
  rtn.addParameter(TSFParameter("max iterations", 500));
  rtn.addParameter(TSFParameter("convergence tolerance", 1.0e-10));
  return rtn;
}

bool BICGSTABSolver::solve(const TSFLinearOperator& op,
			   const TSFVector& b,
			   TSFVector& soln) const 
{
  TSFPreconditioner p = preconditionerFactory_.createPreconditioner(op);
  if (p.isIdentity())
    {
      return solveUnpreconditioned(op, b, soln);
    }
  else if (!p.hasRight())
    {
      TSFLinearOperator A = p.left()*op;
      TSFVector newRHS = b.space().createMember();
      p.left().apply(b, newRHS);
      return solveUnpreconditioned(A, newRHS, soln);
    }
  else if (!p.hasLeft())
    {
      TSFLinearOperator A = op * p.right();
      TSFVector intermediateSoln = soln.space().createMember();
      bool success = solveUnpreconditioned(A, b, intermediateSoln);
      if (success) p.right().apply(intermediateSoln, soln);
      return success;
    }
  else
    {
      TSFLinearOperator A = p.left() * op * p.right();
      TSFVector newRHS;
      p.left().apply(b, newRHS);
      TSFVector intermediateSoln;
      bool success = solveUnpreconditioned(A, newRHS, intermediateSoln);
      if (success) p.right().apply(intermediateSoln, soln);
      return success;
    }
}


bool BICGSTABSolver::solveUnpreconditioned(const TSFLinearOperator& op,
					   const TSFVector& b,
					   TSFVector& soln) const 
{
  int maxiters = params_.getParameter("max iterations").getInt();
  double tol = params_.getParameter("convergence tolerance").getDouble();

  TSFReal normOfB = sqrt(b.dot(b));

  /* check for trivial case of zero rhs */
  if (normOfB < tol) 
    {
      soln = b.space().createMember();
      soln.zero();
      return true;
    }

  /* check for initial zero residual */
  TSFVector x0 = b.copy();
  TSFVector r0 = b.space().createMember();
  TSFVector tmp = b.space().createMember();


  op.apply(x0, tmp);
  r0.subtract(b, tmp);
	
  if (sqrt(r0.dot(r0)) < tol*normOfB) 
    {
      soln = x0;
      return true;
    }

  TSFVector p0 = r0.copy();
  //    p0.randomize();
  TSFVector r0Hat = r0.copy();
  TSFVector xMid = b.space().createMember();
  TSFVector rMid = b.space().createMember();
  TSFVector ArMid = b.space().createMember();
  TSFVector x = b.space().createMember();
  TSFVector r = b.space().createMember();
  TSFVector s = b.space().createMember();
  TSFVector ap = b.space().createMember();

  int myRank = TSFMPI::getRank();

  for (int k=1; k<=maxiters; k++)
    {
      // ap = A*p0
      op.apply(p0, ap);

      TSFReal den = ap.dot(r0Hat);
      if (TSFUtils::chop(sqrt(fabs(den))/normOfB)==0) 
	{
	  TSFOut::rootPrintf("BICGSTAB solver failure mode 1 on iteration %d of %d: Ap*r0Hat=%g normOfB=%g\n",
			     k, maxiters, den, normOfB);
	  TSFError::raise("BICGSTAB::solve failure mode 1");
	}
			
      TSFReal a0 = r0.dot(r0Hat)/den;
			
      xMid = x0 + a0*p0;
      //xMid.axpy(a0, p0, x0);

      rMid = r0 - a0*ap;
      //rMid.axpy(-a0, ap, r0);

      // check for convergence
      TSFReal resid = rMid.norm2()/normOfB;
      if (resid < tol) 
	{
	  soln = xMid; 
	  if (myRank==0 && verbosity_ > 0)
	    {
	      TSFOut::rootPrintf("BICGSTAB converged to resid=%g after %d iterations\n",
				 resid, k);
	    }
	  return true;
	}

      // ArMid = A*rMid
      op.apply(rMid, ArMid);

      den = ArMid.dot(ArMid);
      if (TSFUtils::chop(sqrt(fabs(den))/normOfB)==0)  
	{
	  TSFError::raise("BICGSTAB::solve failure mode 2");
	}

      TSFReal w = rMid.dot(ArMid)/den;
			
      x = xMid + w*rMid;
      //x.axpy(w, rMid, xMid);
			
      r = rMid - w*ArMid;
      //r.axpy(-w, ArMid, rMid);

      // check for convergence
      resid = sqrt(r.dot(r))/normOfB;
      if (resid < tol) 
	{
	  soln = x;
	  if (myRank==0 && verbosity_ > 0)
	    {
	      TSFOut::rootPrintf("BICGSTAB converged to resid=%g after %d iterations\n",
				 resid, k);
	    }
	  return true;
	}

      den = w*(r0.dot(r0Hat));
      if (TSFUtils::chop(sqrt(fabs(den))/normOfB)==0) 
	{
	  TSFError::raise("BICGSTAB::solve failure mode 3");
	}
      TSFReal beta = a0*(r.dot(r0Hat))/den;

      p0 = r + beta*(p0 - w*ap);
      //s.axpy(-w, ap, p0);
      //p0.axpy(beta, s, r);

      r0 = r.copy();
      x0 = x.copy();

      if (myRank==0 && verbosity_ > 1 ) TSFOut::rootPrintf("iteration %d %g\n", k, resid);
    }
  //TSFError::raise("BICGSTAB failed to converge");
  TSFOut::rootPrintf("BICGSTAB failed to converge");
  return false;
	
}

	



