#include <vector>
#include <math.h>
#include "GMRESSolver.h"
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

#include "DenseSerialVectorSpace.h"

using namespace TSF;
using namespace std;


GMRESSolver::GMRESSolver(const TSFReal& tol, int maxIters, int m)
  : TSFLinearSolverBase(), 
  tol_(tol), maxIters_(maxIters), m_(m), preconditionerFactory_()
{;}

GMRESSolver::GMRESSolver(const TSFPreconditionerFactory& pf,
			 const TSFReal& tol, int maxIters, int m)
  : TSFLinearSolverBase(), 
  tol_(tol), maxIters_(maxIters), m_(m), preconditionerFactory_(pf)
{;}


bool GMRESSolver::solve(const TSFLinearOperator& A,
			const TSFVector& b,
			TSFVector& soln) const 
{
  if (verbosity_ > 0) 
    {
      TSFOut::rootPrintln("In GMRES:");
    }
  if (verbosity_ > 2)
    {
      TSFOut::rootPrintln("  GMRES parameters:");
      TSFOut::rootPrintln("    maxIters = " + toString(maxIters_));
      TSFOut::rootPrintln("    tol      = " + toString(tol_)     );
      TSFOut::rootPrintln("    restart  = " + toString(m_)       );
    }
  // TSFTimeMonitor t(myTimer_);
  TSFPreconditioner p = preconditionerFactory_.createPreconditioner(A);
  if (p.isIdentity())
    {
      if (verbosity_ > 1) TSFOut::rootPrintf("  preco = I\n");
      return solveUnpreconditioned(A, b, soln);
    }
  else if (!p.hasRight())
    {
      if (verbosity_ > 1) TSFOut::rootPrintf("  Using left preco\n");
      // cout << "Using left preco" << endl;
      TSFLinearOperator M = p.left()*A;
      TSFVector newRHS = b.space().createMember();
      p.left().apply(b, newRHS);
      return solveUnpreconditioned(M, newRHS, soln);
    }
  else if (!p.hasLeft())
    {
      //cout << "Using right preco" << endl;
      if (verbosity_ > 1) TSFOut::rootPrintln("  Using right preco = \n");
      TSFLinearOperator M = A * p.right();
      TSFVector intermediateSoln = soln.space().createMember();
      //cerr << "In GMRES: going to solve\n";
      bool success = solveUnpreconditioned(M, b, intermediateSoln);
      if (success) p.right().apply(intermediateSoln, soln);
      return success;
    }
  else
    {
      if (verbosity_ > 1) TSFOut::rootPrintf("  Using left and right preco = \n");
      // cout << "Using left and right preco" << endl;
      TSFLinearOperator M = p.left() * A * p.right();
      //TSFVector newRHS = b.space().createMember();
      TSFVector newRHS = p.left().range().createMember();
/*       cerr << "left = \n"; */
/*       p.left().describe(); */
/*       cerr << "b = \n"; */
/*       b.describe(); */
/*       if (b.space() != p.left().domain()) */
/*         { */
/*           cerr << "In GMRES: error in l & r preco\n"; */
/*         } */
/*       else */
/*         { */
/*           cerr << "In GMRES: b is ok\n"; */
/*         } */
      //p.left().apply(b, newRHS);
      newRHS = p.left() * b;
      //cerr << "got new RHS\n";
      TSFVector intermediateSoln = soln.space().createMember();
      //cerr << "   going to solve" << endl;
      bool success = solveUnpreconditioned(M, newRHS, intermediateSoln);
      if (success) p.right().apply(intermediateSoln, soln);
      return success;
    }
}


bool GMRESSolver::solveUnpreconditioned(const TSFLinearOperator& A,
					const TSFVector& b,
					TSFVector& soln) const 
{
  // following GMRES from Matlab
  //cerr << "   Starting iterations\n";
  TSFReal normOfB = sqrt(b.dot(b));

	/* check for trivial case of zero rhs */
  if (normOfB < tol_) 
    {
      soln = b.space().createMember();
      soln.zero();
      return true;
    }

  TSFVectorSpace mSpace = new DenseSerialVectorSpace(m_ + 1);

  TSFVector x0 = soln.space().createMember(); // x0 = 0
  TSFVector r0 = b.space().createMember(); 
  TSFVector vh = b.space().createMember(); 
  TSFVector tmp = b.space().createMember(); 
  TSFVector residVec = b.space().createMember(); 
  TSFVector u = b.space().createMember(); 
  TSFVector vrf = b.space().createMember(); 
  //  TSFVector y = b.space().createMember(); 

  TSFVector h = mSpace.createMember();
  TSFVector f = mSpace.createMember();
  TSFVector q = mSpace.createMember();
  TSFVector mtmp = mSpace.createMember();
  TSFVector y = mSpace.createMember(); 

  vector<TSFVector> V(m_ + 1);
  vector<TSFVector> W(m_ + 1);
  vector<TSFVector> QT(m_ + 1);
  vector<TSFVector> R(m_ + 1);

  int myRank = TSFMPI::getRank();


  for (unsigned int k = 0; k < V.size(); k++) 
    {		
      V[k] = soln.space().createMember(); // V = n x (m+1)
      W[k] = soln.space().createMember(); // W = n x (m+1)
      QT[k] = mSpace.createMember(); // QT = (m+1)x(m+1)
      R[k] = mSpace.createMember(); // R = (m+1)x(m+1)
    }


  TSFReal relTol = tol_ * normOfB; // relative tolerance
  TSFReal phibar;
  TSFReal rt;
  TSFReal c;
  TSFReal s;
  TSFReal temp;

  int CONV = 0; // not converged yet

  // r0 =  b - A*x0;
  A.apply(x0, tmp);
  r0.subtract(b, tmp);
  //  tmp = A*x0;
  //  r0 = b - tmp;
  residVec = r0.copy();
  TSFReal normOfResidVec;
  normOfResidVec = residVec.norm2();
  if (verbosity_ > 1) TSFOut::rootPrintf("  Initial resid = %g\n", normOfResidVec);
  //cout << "initial resid = " << normOfResidVec << endl;

	// Loop i = 1 : maxIters unless convergence (or failure)
  for (int i=0; i<maxIters_; i++)
    {
      h.zero();
      f.zero();
      for (unsigned int z = 0; z < V.size(); z++) 
	{
	  V[z].zero();
	  W[z].zero();
	  QT[z].zero();
	  R[z].zero();
	}

      vh = residVec.copy();
      h[0] = vh.norm2();
      //cout << "new resid norm = " << h[0] << endl;
      double newtemp = 1.0 / h[0];
      //      temp = 1.0 / h[0];
      V[0] = newtemp*vh;
      // V[0] = vh;
      QT[0][0] = 1.0;
      phibar = h[0];

      // inner loop from 0:restart
      for(int j=0; j<m_; j++)
	{
	  u = A*V[j];

	  for( int k=0; k<=j; k++)
	    {
	      h[k] = V[k] * u;
	      u = u - h[k] * V[k];
	    }
					
	  h[j+1] = u.norm2();
	  V[j+1] = u * (1.0 / h[j+1]);

	  for (int k=0; k<=j; k++)
	    {
	      for (int z=0; z<=j; z++) 
		R[k][j] = R[k][j] + QT[k][z] * h[z];
	    }

	  rt = R[j][j];

	  // find cos(theta) and sin(theta) of Givens rotation
	  if (h[j+1] == 0)
	    {
	      c = 1.0; // theta = 0 
	      s = 0.0;
	    }
	  else if (fabs(h[j+1]) > fabs(rt))
	    {
	      temp = rt / h[j+1];
	      s = 1.0 / sqrt(1.0 + fabs(temp)*fabs(temp)); // pi/4 < theta < 3pi/4
	      c = - temp * s;
	    }
	  else
	    {
	      temp = h[j+1] / rt;
	      c = 1.0 / sqrt(1.0 + fabs(temp)*fabs(temp)); // -pi/4 <= theta < 0 < theta <= pi/4
	      s = - temp * c;
	    }
					
	  R[j][j] = c * rt - s * h[j+1]; // left out conj on c and s
	  // cout << "c = " << c << "  s = " << s 
	  //		 << endl;

	  for(int k=0; k<=j; k++)
	    q[k] = QT[j][k];
	  for(int k=0; k<=j; k++)
	    {
	      QT[j][k] = c * q[k];
	      QT[j+1][k] = s * q[k];
	    }
	  QT[j][j+1] = -s;
	  QT[j+1][j+1] =  c;
	  f[j] = c * phibar;
	  phibar = s * phibar;

	  if (j < m_-1)
	    {
	      W[j] = V[j].copy();
	      for(int k=0; k<=j-1; k++)
		W[j] = W[j] - R[k][j]*W[k];
	      W[j] = W[j] * (1.0 / R[j][j]);

	      soln = soln + f[j]*W[j];
	    }
	  else
	    {
	      mtmp.zero();
	      // back solve to get tmp vector to form vrf
	      mtmp[j] = f[j] / R[j][j];
	      for(int k=j-1; k>=0; k--)
		{
		  mtmp[k] = f[k];
		  for(int z=k+1; z<=j; z++)
		    mtmp[k] = mtmp[k] - R[k][z]*mtmp[z];
		  mtmp[k] = mtmp[k] / R[k][k];
		}
							
	      vrf.zero();
	      for(int k=0; k<=j; k++)
		vrf = vrf + mtmp[k] * V[k];
							
	      soln = x0 + vrf;
	    }
					
	  // update current resid norm
	  tmp.zero();
	  tmp = A*soln;
	  residVec = b - tmp;
	  normOfResidVec = residVec.norm2();
	  //TSFOut::rootPrintln("norm of resid " + TSF::toString(normOfResidVec));
	  
	  if (myRank==0 && verbosity_ > 3)
	    {
	      TSFOut::rootPrintf("  Iter = %d\t Current resid = %g\n", j, normOfResidVec/normOfB);
	    }

	  // check for convergence
	  if (normOfResidVec < relTol)
	    {
	      if (j < m_-1)
		{
		  // compute more accurate soln to test convergence
		  y.zero();
		  // back solve to get y(0:j) = R(0:j,0:j) \ f(0:j);
		  y[j] = f[j] / R[j][j];
		  for(int k=j-1; k>=0; k--)
		    {
		      y[k] = f[k];
		      for(int z=k+1; z<=j; z++)
			y[k] = y[k] - R[k][z]*y[z];
		      y[k] = y[k] / R[k][k];
		    }
		  // soln = x0 + V(:,0:j) * y(0:j)
		  soln = x0.copy();
		  for(int k=0; k<=j; k++)
		    soln = soln + y[k] * V[k];
									
		  tmp.zero();
		  tmp = A*soln;
		  residVec = b - tmp;
		  normOfResidVec = residVec.norm2();
		  //cout << "possibly converged" << endl;
		  //cout << "More accurate resid = " << normOfResidVec << endl;
		}
	      // test for convergence
	      if (normOfResidVec < relTol)
		{
		  // we're done
		  CONV = 1; // converged
		  if (myRank==0 && verbosity_ > 0)
		    {
		      TSFOut::rootPrintf("  GMRES converged (1)! resid = %g\n", normOfResidVec/normOfB);
		      TSFOut::rootPrintf("  iter = %d\n", j);
		    }
		  return true;
		  // break;
		}
	    }
	}
			
      if (!CONV)
	{
	  // not converged yet; update x0 and resid and restart
	  x0 = soln.copy();
	  tmp.zero();
	  tmp = A*x0;
	  residVec = b - tmp;
	  normOfResidVec = residVec.norm2();
      if (verbosity_ > 0)
        {
          TSFOut::rootPrintf("  GMRES not converged; current scaled resid = %g\n ", 
                             normOfResidVec/normOfB);
          //cout << "not converged; current scaled resid = " << normOfResidVec/normOfB << endl;
        }
	}
      else
	{
	  if (myRank==0 && verbosity_ > 0)
	    {
	      TSFOut::rootPrintf("  GMRES converged (2)! scaled norm of resid = %g\n", normOfResidVec/normOfB);
	    }
	  return true;
	}
      //	break;
			
    }
	
	
	
  return false;
}
	



