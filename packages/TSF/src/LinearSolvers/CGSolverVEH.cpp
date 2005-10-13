#include <vector>
#include <math.h>
#include <float.h>
// #include <stdio.h>
// #include <stdlib.h>
// #include <errno.h>
// #include <fcntl.h>
// #include <unistd.h>
// #include <string.h>
// #include <fstream>
#include "CGSolverVEH.h"
#include "TSFUtils.h"
#include "TSFLinearOperator.h"
#include "TSFLinearOperatorBase.h"
#include "TSFVectorSpaceBase.h"
#include "TSFVectorBase.h"
#include "TSFVector.h"
#include "TSFVectorSpace.h"
#include "TSFVectorType.h"

#include "TSFBlas.h"
#include "TSFError.h"
#include "TSFOut.h"
#include "TSFMPI.h"
#include "TSFDeferredLinearCombination.h"
#include "TSFPreconditionerFactoryBase.h"
#include "TSFPreconditionerFactory.h"
#include "DenseSerialVectorSpace.h"
#include "DenseSerialVector.h"
#include "StrUtils.h"
#include "TSFOperatorSourceBase.h"

// #include "PetraVectorSpace.h"
// #include "PetraVector.h"
// #include "PetraVectorType.h"
// #include "PetraMatrix.h"

// #include "Epetra_MpiComm.h"
// #include "Epetra_Comm.h"
// #include "Epetra_Map.h"
// #include "Epetra_MultiVector.h"
// #include "Epetra_Vector.h"
// #include "Epetra_CrsMatrix.h"

// CG: modified to include breakdown checking and full checkpointing


using namespace TSF;
using namespace std;


CGSolverVEH::CGSolverVEH(const TSFReal& tol, int maxIters)
  : TSFLinearSolverBase(),
    tol_(tol), 
    maxIters_(maxIters),
    preconditionerFactory_()
{;}


CGSolverVEH::CGSolverVEH(const TSFPreconditionerFactory& pf, 
                         const TSFReal& tol, 
			 int maxIters) 
  : TSFLinearSolverBase(), 
    tol_(tol), 
    maxIters_(maxIters), 
    preconditionerFactory_(pf)
{;}


CGSolverVEH::~CGSolverVEH(){;}

bool CGSolverVEH::solve(const TSFLinearOperator& A,
                        const TSFVector& b,
                        TSFVector& soln) const 
{

  if (verbosity_ > 0) 
    {
      TSFOut::rootPrintln("In CGSolverVEH:");
    }
  if (verbosity_ > 2) // print initial setup information
    {
      TSFOut::rootPrintln("  CG parameters:");
      TSFOut::rootPrintln("    maxIters = " + toString(maxIters_));
      TSFOut::rootPrintln("    tol      = " + toString(tol_)     );
    }

  // Don't have symmetric preconditioners available yet
  TSFPreconditioner p = preconditionerFactory_.createPreconditioner(A);
  if (p.isIdentity())
    {
      if (verbosity_ > 1) TSFOut::rootPrintf("  preco = I\n");
      bool success = solveUnpreconditioned(A, b, soln);
      cerr << "success? " << success << endl;
      return success;
      //      return solveUnpreconditioned(A, b, soln);
    }
  else if (!p.hasRight())
    {
      if (verbosity_ > 1) TSFOut::rootPrintf("  Using left preco\n");
      TSFLinearOperator M = p.left()*A;
      TSFVector newRHS = b.space().createMember();
      p.left().apply(b, newRHS);
      return solveUnpreconditioned(M, newRHS, soln);
    }
  else if (!p.hasLeft())
    {
      if (verbosity_ > 1) TSFOut::rootPrintln("  Using right preco\n");
      TSFLinearOperator M = A * p.right();
      TSFVector intermediateSoln;
      bool success = solveUnpreconditioned(M, b, intermediateSoln);
      if (success) p.right().apply(intermediateSoln, soln);
      return success;
    }
  else
    {
      if (verbosity_ > 1) TSFOut::rootPrintf("  Using left and right preco\n");
      TSFLinearOperator M = p.left() * A * p.right();
      TSFVector newRHS;
      p.left().apply(b, newRHS);
      TSFVector intermediateSoln;
      bool success = solveUnpreconditioned(M, newRHS, intermediateSoln);
      if (success) p.right().apply(intermediateSoln, soln);
      return success;
    }
}




bool CGSolverVEH::solveUnpreconditioned(const TSFLinearOperator& A,
					const TSFVector& b,
					TSFVector& soln) const 
{


  int myRank = TSFMPI::getRank();
  int numProcs = TSFMPI::getNProc();
	
  TSFReal normOfB = b.norm2();

  TSFReal breakdown_tol = DBL_EPSILON;

  /* Starting with x0 = 0 and r = b */
  TSFVector x0 = b.space().createMember();
  TSFVector r = b.copy();
  TSFVector p = r.copy();
  
  // 0th iteration
  TSFVector Ap = A*p;
  TSFReal r20 = r*r;
  TSFReal pAp = p*Ap;
  if (TSFUtils::chop(pAp)==0) 
    {
      TSFError::raise("failure mode 1 in CG");
    }
  TSFReal alpha = r20/pAp;
  x0 = x0 + alpha*p;
  r = r - alpha*Ap;
  TSFReal r21 = r*r;
  TSFReal beta = r21/r20;
  p = r + beta*p;
  
  TSFReal normOfR = r.norm2();

  if (verbosity_ > 1) 
    TSFOut::rootPrintf("  Initial (unscaled) resid = %g\n", 
                       normOfR);
  
  if (verbosity_ > 1) 
    TSFOut::rootPrintf("  Initial scaled resid = %g\n", 
                       normOfR/normOfB);
  
  for (int iter=0; iter<=maxIters_; iter++)
    {
      normOfR = r.norm2();
      // write residual data into file for Edward plotting code
      //			if (myRank == 0)
      //      	{
      //      	  logMyData(0, realIterationCount_, normOfR);
      //				}
      // write out residuals for plotting at end
      //       if (myRank == 0)
      //         {
      //           char CGresid[100];
      //           sprintf(CGresid, "CGresid.txt");
      //           ofstream ofsresid(CGresid, ios::ate);
      //           if (ofsresid.bad())
      //             {
      //               TSFError::raise("CG failed to open output file");
      //               return false;
      //             }
      //           ofsresid << normOfR << endl;
      //           ofsresid.close();
      //         }

      Ap = A*p;
      pAp = p*Ap;

      //   if (TSFUtils::chop(pAp)==0) 
      //      {
      //	TSFError::raise("failure mode 1 in CG");
      //      }

      if (fabs(pAp) < breakdown_tol)
        {
          // check for breakdown -- following breakdown check in AZTEC
          TSFReal breakdown_check = 100.0 * p.norm2() * Ap.norm2() * DBL_EPSILON;
          if (fabs(pAp) < breakdown_check)
            {
              // CG has hit breakdown
              TSFError::raise("CG breakdown 2");
            }
          else
            {
              // update breakdown_tol and continue
              breakdown_tol = 0.1 * fabs(pAp);
            }
        }
			
      r20 = r*r;
      alpha = r20/pAp;
      x0 = x0 + alpha*p;
      r = r - alpha*Ap;
      normOfR = r.norm2();

      if (myRank==0 && verbosity_ > 3)
	{
	  TSFOut::rootPrintf("  Iter = %d\t Current scaled resid = %g\n", 
			     iter, normOfR/normOfB);
	}
      
      r21 = r*r;

      if (sqrt(r21) < tol_*normOfB) 
        {
          soln = x0.copy();
	  if (verbosity_ > 0)
	    TSFOut::rootPrintf("CG converged (1) in %d iters: final scaled resid = %g\n", 
			       iter+1, sqrt(r21)/normOfB);

          return true;
        }

      beta = r21/r20;
      p = r + beta*p;

			
    }
  cerr << "failed to converge" << endl;
  // TSFError::raise("CG failed to converge");
  return false;
}
	



